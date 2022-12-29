//
// Created by 胡栋月 on 6/12/22.
//
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>

#include "head/geometry.h"
#include "head/light.h"
#include "head/specialShape.h"
#include "head/triangleMesh.h"
#include "head/tools.h"

#if !defined METHOD_GEOMETRY || !defined CULLING || !FACE_NORMAL || !SMOOTH_SHADING
//#define SMOOTH_SHADING
//#define METHOD_GEOMETRY
//#define CULLING
//#define FACE_NORMAL

/**
 * compute the direction for reflection
 */
Vec3f reflection(const Vec3f &incidentDir, const Vec3f &Nhit) {
    return incidentDir - 2 * incidentDir.dotProduct(Nhit) * Nhit;
}

/**
 *  compute the refraction
 *
 */

Vec3f refraction(const Vec3f &incidentDir, const Vec3f &Nhit, const float &ior) {
    Vec3f normal = Nhit;
    float etaIncident = 1, etaRefract = ior;
    float NdotNormal = clamp(-1, 1, normal.dotProduct(incidentDir)) ;
    float cosIncident;
    if(NdotNormal < 0) { // from air to another medium
        cosIncident = -NdotNormal;
    } else { // from another medium to air
        cosIncident = NdotNormal;
        normal = -normal;
        std::swap(etaIncident, etaRefract);
    }

    float eta = etaIncident / etaRefract;

    float k = 1 - eta * eta * (1 - cosIncident * cosIncident);

    if(k < 0) { // totally reflection
        return 0;
    } else {
        return (eta * incidentDir + (eta * cosIncident - sqrt(k)) * normal).normalize();
    }
}

/**
 *  compute the portion for refraction and reflection
 *  Fresnel law
 */
void fresnel(const Vec3f &incidentDir, const Vec3f &Nhit, const float &ior, float &portionForReflection) {
    float cosIncident = clamp(-1, 1, incidentDir.dotProduct(Nhit));
    float etaIncident = 1, etaRefraction = ior;
    if (cosIncident > 0) { std::swap(etaIncident, etaRefraction); }
    else cosIncident = fabsf(cosIncident);
    float sinRefraction = etaIncident / etaRefraction * sqrtf(std::max(0.f, 1 - cosIncident * cosIncident));
    // Total internal reflection
    if (sinRefraction >= 1) {
        kr = 1;
    }
    else {
        float cosRefraction = sqrtf(std::max(0.f, 1 - sinRefraction * sinRefraction));
        float Rs = ((etaRefraction * cosIncident) - (etaIncident * cosRefraction)) / ((etaRefraction * cosIncident) + (etaIncident * cosRefraction));
        float Rp = ((etaIncident * cosIncident) - (etaRefraction * cosRefraction)) / ((etaIncident * cosIncident) + (etaRefraction * cosRefraction));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
}

/**
 * Test each object for intersection
 * @param orig
 * @param dir
 * @param objects
 * @param tNear
 * @param hitObject
 * @return
 */
bool trace(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects, IntersecInfo &intersecInfo, const RayType &rayType = kPrimaryRay) {
    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    for(; iter != objects.end(); ++iter) {
        float t = kInfinity;
        uint32_t indexForFace;
        Vec2f uvForFace;

        if((*iter)->intersect(orig, dir, t, indexForFace, uvForFace) && t < intersecInfo.tNear) {
//            std::cout << "intersect with something" << std::endl;
            if(rayType == kShadowRay && (*iter)->materialType == kReflectionAndRefraction) continue;
            intersecInfo.hitObject = iter->get();
            intersecInfo.tNear = t;
            intersecInfo.index = indexForFace;
            intersecInfo.uv = uvForFace;
        }
    }

    return (intersecInfo.hitObject != nullptr);
}

/**
 * Cast a ray with specific origin and direction, then compute the color in hit point.
 * @param orig
 * @param dir
 * @param objects
 * @return
 */
Vec3f castRay(const Vec3f &orig, const Vec3f &dir,
              const std::vector<std::unique_ptr<Object>> &objects,
              const std::vector<std::unique_ptr<Light>> &lights,
              const Options &options,
              const uint32_t depth = 0) {
    if(depth > options.maxDepth) return options.backgroundColor;

    Vec3f hitColor = 0;
    IntersecInfo intersecInfo;

    if(trace(orig, dir, objects, intersecInfo)) {
        Vec3f Phit = orig + dir * intersecInfo.tNear;
        Vec3f Nhit;
        Vec2f Texhit;
        intersecInfo.hitObject->getSurfaceData(Phit, dir, intersecInfo.index, intersecInfo.uv, Nhit, Texhit);

        // no lights
        if(lights.size() == 0) {
            float scale = 4;
            float pattern = (fmodf(Texhit.x * scale, 1) > 0.5) ^ (fmodf(Texhit.y * scale, 1) > 0.5);
            hitColor = hitColor + std::max(0.f, Nhit.dotProduct(-dir)) * mix(intersecInfo.hitObject->color, intersecInfo.hitObject->color * 0.8, pattern);
            return hitColor;
        }

        switch (intersecInfo.hitObject->materialType) {
            case kDiffuse: {
                for(uint32_t i=0; i<lights.size(); ++i) {
                    Vec3f lightDir, lightIntensity;
                    IntersecInfo tmpIntersectInfo;
                    lights[i]->getDirectionAndIntensity(Phit, lightDir, lightIntensity, tmpIntersectInfo.tNear);
                    // compute shadow
                    bool vis = !trace(Phit + Nhit * options.bias, -lightDir, objects, tmpIntersectInfo, kShadowRay);

                    float scale = 4;
                    float pattern = (fmodf(Texhit.x * scale, 1) > 0.5) ^ (fmodf(Texhit.y * scale, 1) > 0.5);
                    hitColor = hitColor + vis * lightIntensity * (intersecInfo.hitObject->albedo / M_PI * std::max(0.f, Nhit.dotProduct(-dir)) *
                            mix(intersecInfo.hitObject->color, intersecInfo.hitObject->color * 0.8, pattern));
//                    hitColor = hitColor + vis * lightIntensity * (intersecInfo.hitObject->albedo / M_PI *
//                                              std::max(0.f, Nhit.dotProduct(-lightDir)));
                }
                break;
            }

            case kReflection: {
                Vec3f lightDir = -reflection(dir, Nhit);
                IntersecInfo tmpIntersectInfo;
                hitColor = hitColor + 0.8 * castRay(Phit + Nhit * options.bias, -lightDir, objects, lights, options, depth + 1);
                break;
            }

            case kReflectionAndRefraction: {
                Vec3f refractionColor = 0, reflectionColor = 0;
                // compute fresnel
                float portionForReflection;
                fresnel(dir, Nhit, intersecInfo.hitObject->ior, portionForReflection);
                bool outside = dir.dotProduct(Nhit) < 0;
                Vec3f bias = options.bias * Nhit;
                // not total internal reflection
                if (portionForReflection < 1) {
                    Vec3f refractionDirection = refraction(dir, Nhit, intersecInfo.hitObject->ior);
                    Vec3f refractionRayOrig = outside ? Phit - bias : Phit + bias;
                    refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, depth + 1);
                }

                Vec3f reflectionDirection = reflection(dir, Nhit).normalize();
                Vec3f reflectionRayOrig = outside ? Phit + bias : Phit - bias;
                reflectionColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1);

                // add the color for reflection and refraction
                hitColor = hitColor + reflectionColor * portionForReflection + refractionColor * (1 - portionForReflection);
                break;
            }
            default: {
                std::cout << "exceptional type of material!!!" << std::endl;
                break;
            }
        }

    } else hitColor = options.backgroundColor;

    return hitColor;
}

/**
 * Render the whole scene with options and objects info
 * @param options
 * @param objects
 */
void render(const Options &options, const std::vector<std::unique_ptr<Object>> &objects, const std::vector<std::unique_ptr<Light>> &lights) {
    Vec3f *framebuffer = new Vec3f[options.width * options.height];
    Vec3f *pix = framebuffer;

    float scale = tan(deg2rad(options.fov * 0.5));
    float imageAspectRatio = options.width / (float) options.height;

    Vec3f orig;
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
    auto timeStart = std::chrono::high_resolution_clock::now();
    for(uint32_t j=0; j<options.height; ++j) {
        for(uint32_t i=0; i<options.width; ++i) {
            float x = (2 * (i + 0.5) / (float) options.width - 1) * scale * imageAspectRatio;
            float y = (1 - 2 * (j + 0.5) / (float) options.height) * scale;

            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();
            *(pix++) = castRay(orig, dir, objects, lights, options, 0);
        }
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
    fprintf(stderr, "\rDone: %.2f (sec)\n", passedTime / 1000);

    std::ofstream ofs("./out/out.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << options.width << " " << options.height << "\n225\n";
    for(uint32_t i=0; i<options.height * options.width; ++i) {
        char r = (char)(255*clamp(0, 1, framebuffer[i].x));
        char g = (char)(255*clamp(0, 1, framebuffer[i].y));
        char b = (char)(255*clamp(0, 1, framebuffer[i].z));
        ofs << r << g << b;
    }

    ofs.close();

    delete [] framebuffer;
}


int main(int argc, char **argv) {
    std::vector<std::unique_ptr<Object>> objects;
    std::vector<std::unique_ptr<Light>> lights;

    gen.seed(0);

    // setting up options
    Options options;
#if 1 // control the resolution
    options.width = 640;
    options.height = 480;
#else
    options.width = 2560;
    options.height = 1600;
#endif
    options.backgroundColor = kDefaultBackgroundColor;
    options.fov = 51.52;
    options.bias = 1e-4;

#if 0 // version 0.1
    options.cameraToWorld = Matrix44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);
    options.bias = 10;

    uint32_t numSpheres = 3;
    // sphere
    for (uint32_t i = 0; i < numSpheres; ++i) {
        Vec3f randPos((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, (0.5 + dis(gen) * 10));
        float randRadius = (0.5 + dis(gen) * 1);
        objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius)));
    }

    // triangle
    Vec3f v0((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, (1 + dis(gen) * 10));
    Vec3f v1(0.5 * 10, (0.5 - dis(gen)) * 10, (2 + dis(gen) * 10));
    Vec3f v2((0.5 - dis(gen)) * 10, 0.5 * 10, (3 + dis(gen) * 10));
    objects.push_back(std::unique_ptr<Object>(new Triangle(v0, v1, v2)));

    // disk
    Vec3f normal1(0.2, 0.5 * 10, (1 + dis(gen) * 10));
    Vec3f center1(0.2, 0.5 * 8, 0.1 );
    float radius1 = 2;
    objects.push_back(std::unique_ptr<Object>(new Disk(normal1, center1, radius1)));

    Vec3f normal2(0.2, 0.5 * 10, 0);
    Vec3f center2(0, 0.5 * 8, 0.1 );
    float radius2 = 2;
    objects.push_back(std::unique_ptr<Object>(new Disk(normal2, center2, radius2)));

#elif 0 // version 0.2

    Matrix44f tmp = Matrix44f(0.707107, -0.331295, 0.624695, 0, 0, 0.883452, 0.468521, 0, -0.707107, -0.331295, 0.624695, 0, -1.63871, -5.747777, -40.400412, 1);
    options.cameraToWorld = tmp.inverse();

    TriangleMesh *mesh = loadPolyMeshFromFile("./geometry/cow.geo");
    if (mesh != nullptr) objects.push_back(std::unique_ptr<Object>(mesh));

# elif 0// test Smooth Shading
    options.cameraToWorld[3][2] = 7;
    TriangleMesh *meshBall = generatePolyShphere(1.5, 16);
    if(meshBall != nullptr) objects.push_back(std::unique_ptr<Object>(meshBall));

    // version 0.3
# elif 0 // test distant light shadow
    Matrix44f l2w(1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w)));

    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f();
    options.bias = 0.1;

    Vec3f randPos(0, 2.5, -15);
    float randRadius = 2;
    objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius)));

    Vec3f randPos1(0, -4, -16.5);
    float randRadius1 = 4;
    objects.push_back(std::unique_ptr<Object>(new Sphere(randPos1, randRadius1)));
#elif 0 // test point light shadow
    Matrix44f l2w(2, 0, 0, 0, 0, 1, 0, 0, 0, 10, 0, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new PointLight(Matrix44f(), Vec3f(5, 2, 3), 500)));

    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f();
    options.bias = 0.15;

    Vec3f randPos(0, 1, -4);
    float randRadius = 1;
    objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius, 0.2)));

//    Vec3f randPos1(0, -4, -16.5);
//    float randRadius1 = 4;
//    objects.push_back(std::unique_ptr<Object>(new Sphere(randPos1, randRadius1)));
#elif 0 // test several lights
    Matrix44f l2w(1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w)));

    Matrix44f l2w1(3, 2, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(Matrix44f())));

    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f();
    options.bias = 0.1;


    Vec3f randPos(0, 3.5, -15);
    float randRadius = 2;
    objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius, 0.18, Matrix44f(), kDiffuse, "sphere1")));

    Vec3f randPos1(0, -12, -40);
    float randRadius1 = 25;
    objects.push_back(std::unique_ptr<Object>(new Sphere(randPos1, randRadius1,0.18, Matrix44f(), kDiffuse, "sphere2")));
#elif 1 // simple plane example (patterns)
    options.fov = 36.87;
    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f(0.707107, 0, -0.707107, 0, -0.331295, 0.883452, -0.331295, 0, 0.624695, 0.468521, 0.624695, 0, 28, 21, 28, 1);

    TriangleMesh *mesh = loadPolyMeshFromFile("./geometry/plane.geo");
    if (mesh != nullptr) {
        // mesh->smoothShading = false;
        objects.push_back(std::unique_ptr<Object>(mesh));
    }
    Matrix44f l2w(11.146836, -5.781569, -0.0605886, 0, -1.902827, -3.543982, -11.895445, 0, 5.459804, 10.568624, -4.02205, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 1000)));
#elif 0 // multiple glasses example
    options.fov = 36.87;
    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f(0.999945, 0, 0.0104718, 0, 0.00104703, 0.994989, -0.0999803, 0, -0.0104193, 0.0999858, 0.994934, 0, -0.978596, 17.911879, 75.483369, 1);

    TriangleMesh *mesh = loadPolyMeshFromFile("./geometry/glasses.geo");
    if (mesh != nullptr) {
        mesh->materialType = kReflectionAndRefraction;
        mesh->ior = 1.3;
        objects.push_back(std::unique_ptr<Object>(mesh));
    }

    TriangleMesh *mesh1 = loadPolyMeshFromFile("./geometry/backdrop1.geo");
    if (mesh1 != nullptr) {
        mesh1->materialType = kDiffuse;
        mesh1->albedo = 0.18;
        objects.push_back(std::unique_ptr<Object>(mesh1));
    }

    Matrix44f l2w(0.95292, 0.289503, 0.0901785, 0, -0.0960954, 0.5704, -0.815727, 0, -0.287593, 0.768656, 0.571365, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 15)));

#elif 0
    // glass and pen example
    // setting up options
    options.fov = 36.87;
    options.maxDepth = 10;
    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f(-0.972776, 0, -0.231748, 0, -0.114956, 0.8683, 0.482536, 0, 0.201227, 0.49604, -0.844661, 0, 6.696465, 22.721296, -30.097976, 1);

    TriangleMesh *mesh1 = loadPolyMeshFromFile("./geometry/backdrop.geo");
    if (mesh1 != nullptr) {
        mesh1->materialType = kDiffuse;
        objects.push_back(std::unique_ptr<Object>(mesh1));
    }

    TriangleMesh *mesh3 = loadPolyMeshFromFile("./geometry/cylinder.geo");
    if (mesh3 != nullptr) {
        mesh3->materialType = kReflectionAndRefraction;
        mesh3->ior = 1.5;
        objects.push_back(std::unique_ptr<Object>(mesh3));
    }

    TriangleMesh *mesh4 = loadPolyMeshFromFile("./geometry/pen.geo");
    if (mesh4 != nullptr) {
        mesh4->materialType = kDiffuse;
        mesh4->albedo = 0.18;
        objects.push_back(std::unique_ptr<Object>(mesh4));
    }

    Matrix44f xform1;
    xform1[3][0] = -1.2;
    xform1[3][1] = 6;
    xform1[3][2] = -3;
    Sphere *sph1 = new Sphere(0, 5, xform1);
    if(sph1 != nullptr) {
        sph1->materialType = kReflectionAndRefraction;
        sph1->ior = 1.7;
        objects.push_back(std::unique_ptr<Object>(sph1));
    }

    Matrix44f l2w(11.146836, -5.781569, -0.0605886, 0, -1.902827, -3.543982, -11.895445, 0, 5.459804, 10.568624, -4.02205, 0, 0, 0, 0, 1);
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 20)));

# endif

    render(options, objects, lights);

    return 0;
}

#endif