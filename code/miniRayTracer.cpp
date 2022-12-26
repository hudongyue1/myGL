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

#if !defined METHOD_GEOMETRY || !defined CULLING || !FACE_NORMAL
//#define METHOD_GEOMETRY
#define CULLING
//#define FACE_NORMAL

/**
 * Test each object for intersection
 * @param orig
 * @param dir
 * @param objects
 * @param tNear
 * @param hitObject
 * @return
 */
bool trace(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects, IntersecInfo &intersecInfo) {
    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    for(; iter != objects.end(); ++iter) {
        float t = kInfinity;
        uint32_t indexForFace;
        Vec2f uvForFace;

        if((*iter)->intersect(orig, dir, t, indexForFace, uvForFace) && t < intersecInfo.tNear) {
//            std::cout << "intersect with something" << std::endl;
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
              const Options &options) {
    Vec3f hitColor = 0;
    IntersecInfo intersecInfo;

    if(trace(orig, dir, objects, intersecInfo)) {
        Vec3f Phit = orig + dir * intersecInfo.tNear;
        Vec3f Nhit;
        Vec2f Texhit;
        intersecInfo.hitObject->getSurfaceData(Phit, dir, intersecInfo.index, intersecInfo.uv, Nhit, Texhit);


        for(uint32_t i=0; i<lights.size(); ++i) {
            Vec3f lightDir, lightIntensity;
            IntersecInfo tmpIntersecInfo;
            lights[i]->getDirectionAndIntensity(Phit, lightDir, lightIntensity, tmpIntersecInfo.tNear);

            //        std::cout << "test shadow ray" << std::endl;
            bool vis = !trace(Phit + Nhit * options.bias, -lightDir, objects, tmpIntersecInfo);

            float scale = 4;
            float pattern = (fmodf(Texhit.x * scale, 1) > 0.5) ^ (fmodf(Texhit.y * scale, 1) > 0.5);

#ifdef NO_LIGHT
            hitColor = hitColor + std::max(0.f, Nhit.dotProduct(-dir)) * mix(hitObject->color, hitObject->color * 0.8, pattern);
#else
            hitColor = vis * hitColor + lightIntensity * (vis * intersecInfo.hitObject->albedo / M_PI *
                    std::max(0.f, Nhit.dotProduct(-lightDir)) *
                    mix(intersecInfo.hitObject->color, intersecInfo.hitObject->color * 0.8, pattern));
//        hitColor = vis * hitObject->albedo / M_PI * light->intensity * std::max(0.f, Nhit.dotProduct(-LightDir));
#endif
//        hitColor = std::max(0.f, Nhit.dotProduct(-LightDir));
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
            *(pix++) = castRay(orig, dir, objects, lights, options);
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

    TriangleMesh *mesh = loadPolyMeshFromFile("./resource/cow.geo");
    if (mesh != nullptr) objects.push_back(std::unique_ptr<Object>(mesh));

# elif 0// test Smooth Shading
    options.cameraToWorld[3][2] = 7;
    TriangleMesh *meshBall = generatePolyShphere(1.5, 16);
    if(meshBall != nullptr) objects.push_back(std::unique_ptr<Object>(meshBall));

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
#else // test point light shadow
    lights.push_back(std::unique_ptr<Light>(new PointLight(Matrix44f(), Vec3f(10, 2, 3), 1000)));

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

# endif

    render(options, objects, lights);

    return 0;
}

#endif