#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <chrono>

#include "head/geometry.h"

static const float kInfinity = std::numeric_limits<float>::max();
static const float kEpsilon = 1e-8;
static const Vec3f kDefaultBackgroundColor = Vec3f(0.235294, 0.67451, 0.843137);

inline
float clamp(const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }

struct Options
{
    uint32_t width = 640;
    uint32_t height = 480;
    float fov = 90;
    Vec3f backgroundColor = kDefaultBackgroundColor;
    Matrix44f cameraToWorld;
    float bias = 0.0001;
    uint32_t maxDepth = 5;
};

enum MaterialType { kDiffuse, kReflection, kReflectionAndRefraction };

class Object
{
public:
    Object(const Matrix44f &o2w) : objectToWorld(o2w), worldToObject(o2w.inverse()) {}
    virtual ~Object() {}
    virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;
    virtual void getSurfaceProperties(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec3f &, Vec2f &) const = 0;
    Matrix44f objectToWorld, worldToObject;
    const char *name;
    MaterialType type = kDiffuse;
    float ior = 1;
    Vec3f albedo = 0.18;
};
bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0) return false;
    else if (discr == 0) {
        x0 = x1 = - 0.5 * b / a;
    }
    else {
        float q = (b > 0) ?
                  -0.5 * (b + sqrt(discr)) :
                  -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }

    return true;
}
class Sphere : public Object
{
public:
    Sphere(const Matrix44f &o2w, const float &r) : Object(o2w), radius(r), radius2(r *r)
    {
        o2w.multVecMatrix(Vec3f(0), center);
        this->name = "sphere";
    }
    bool intersect(
            const Vec3f &orig,
            const Vec3f &dir,
            float &tNear,
            uint32_t &triIndex,  //not used for sphere
            Vec2f &uv) const     //not used for sphere
    {
        float t0, t1;  //solutions for t if the ray intersects
        // analytic solution
        Vec3f L = orig - center;
        float a = dir.dotProduct(dir);
        float b = 2 * dir.dotProduct(L);
        float c = L.dotProduct(L) - radius2;
        if (!solveQuadratic(a, b, c, t0, t1)) return false;

        if (t0 > t1) std::swap(t0, t1);

        if (t0 < 0) {
            t0 = t1;  //if t0 is negative, let's use t1 instead
            if (t0 < 0) return false;  //both t0 and t1 are negative
        }

        tNear = t0;

        return true;
    }
    void getSurfaceProperties(
            const Vec3f &hitPoint,
            const Vec3f &viewDirection,
            const uint32_t &triIndex,
            const Vec2f &uv,
            Vec3f &hitNormal,
            Vec2f &hitTextureCoordinates) const
    {
        hitNormal = hitPoint - center;
        hitNormal.normalize();
        // In this particular case, the normal is simular to a point on a unit sphere
        // centred around the origin. We can thus use the normal coordinates to compute
        // the spherical coordinates of Phit.
        // atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
        // acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
        hitTextureCoordinates.x = (1 + atan2(hitNormal.z, hitNormal.x) / M_PI) * 0.5;
        hitTextureCoordinates.y = acosf(hitNormal.y) / M_PI;
    }
    float radius, radius2;
    Vec3f center;
};

bool rayTriangleIntersect(
        const Vec3f &orig, const Vec3f &dir,
        const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
        float &t, float &u, float &v)
{
    Vec3f v0v1 = v1 - v0;
    Vec3f v0v2 = v2 - v0;
    Vec3f pvec = dir.crossProduct(v0v2);
    float det = v0v1.dotProduct(pvec);

    // ray and triangle are parallel if det is close to 0
    if (fabs(det) < kEpsilon) return false;

    float invDet = 1 / det;

    Vec3f tvec = orig - v0;
    u = tvec.dotProduct(pvec) * invDet;
    if (u < 0 || u > 1) return false;

    Vec3f qvec = tvec.crossProduct(v0v1);
    v = dir.dotProduct(qvec) * invDet;
    if (v < 0 || u + v > 1) return false;

    t = v0v2.dotProduct(qvec) * invDet;

    return (t > 0) ? true : false;
}

class TriangleMesh : public Object
{
public:
    // Build a triangle mesh from a face index array and a vertex index array
    TriangleMesh(
            const Matrix44f &o2w,
            const uint32_t nfaces,
            const std::unique_ptr<uint32_t []> &faceIndex,
            const std::unique_ptr<uint32_t []> &vertsIndex,
            const std::unique_ptr<Vec3f []> &verts,
            std::unique_ptr<Vec3f []> &normals,
            std::unique_ptr<Vec2f []> &st) :
            Object(o2w),
            numTris(0)
    {
        this->name = "trianglemesh";
        uint32_t k = 0, maxVertIndex = 0;
        // find out how many triangles we need to create for this mesh
        for (uint32_t i = 0; i < nfaces; ++i) {
            numTris += faceIndex[i] - 2;
            for (uint32_t j = 0; j < faceIndex[i]; ++j)
                if (vertsIndex[k + j] > maxVertIndex)
                    maxVertIndex = vertsIndex[k + j];
            k += faceIndex[i];
        }
        maxVertIndex += 1;

        // allocate memory to store the position of the mesh vertices
        P = std::unique_ptr<Vec3f []>(new Vec3f[maxVertIndex]);
        for (uint32_t i = 0; i < maxVertIndex; ++i) {
            objectToWorld.multVecMatrix(verts[i], P[i]);
        }

        // allocate memory to store triangle indices
        trisIndex = std::unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);
        uint32_t l = 0;
        N = std::unique_ptr<Vec3f []>(new Vec3f[numTris * 3]);
        sts = std::unique_ptr<Vec2f []>(new Vec2f[numTris * 3]);
        Matrix44f transformNormals = worldToObject.transpose();
        // generate the triangle index array and set normals and st coordinates
        for (uint32_t i = 0, k = 0; i < nfaces; ++i) {  //for each  face
            for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) {  //for each triangle in the face
                trisIndex[l] = vertsIndex[k];
                trisIndex[l + 1] = vertsIndex[k + j + 1];
                trisIndex[l + 2] = vertsIndex[k + j + 2];
                transformNormals.multDirMatrix(normals[k], N[l]);
                transformNormals.multDirMatrix(normals[k + j + 1], N[l + 1]);
                transformNormals.multDirMatrix(normals[k + j + 2], N[l + 2]);
                N[l].normalize();
                N[l + 1].normalize();
                N[l + 2].normalize();
                sts[l] = st[k];
                sts[l + 1] = st[k + j + 1];
                sts[l + 2] = st[k + j + 2];
                l += 3;
            }
            k += faceIndex[i];
        }
    }
    // Test if the ray interesests this triangle mesh
    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tNear, uint32_t &triIndex, Vec2f &uv) const
    {
        uint32_t j = 0;
        bool isect = false;
        for (uint32_t i = 0; i < numTris; ++i) {
            const Vec3f &v0 = P[trisIndex[j]];
            const Vec3f &v1 = P[trisIndex[j + 1]];
            const Vec3f &v2 = P[trisIndex[j + 2]];
            float t = kInfinity, u, v;
            if (rayTriangleIntersect(orig, dir, v0, v1, v2, t, u, v) && t < tNear) {
                tNear = t;
                uv.x = u;
                uv.y = v;
                triIndex = i;
                isect = true;
            }
            j += 3;
        }

        return isect;
    }
    void getSurfaceProperties(
            const Vec3f &hitPoint,
            const Vec3f &viewDirection,
            const uint32_t &triIndex,
            const Vec2f &uv,
            Vec3f &hitNormal,
            Vec2f &hitTextureCoordinates) const
    {
        if (smoothShading) {
            // vertex normal
            const Vec3f &n0 = N[triIndex * 3];
            const Vec3f &n1 = N[triIndex * 3 + 1];
            const Vec3f &n2 = N[triIndex * 3 + 2];
            hitNormal = (1 - uv.x - uv.y) * n0 + uv.x * n1 + uv.y * n2;
        }
        else {
            // face normal
            const Vec3f &v0 = P[trisIndex[triIndex * 3]];
            const Vec3f &v1 = P[trisIndex[triIndex * 3 + 1]];
            const Vec3f &v2 = P[trisIndex[triIndex * 3 + 2]];
            hitNormal = (v1 - v0).crossProduct(v2 - v0);
        }

        // doesn't need to be normalized as the N's are normalized but just for safety
        hitNormal.normalize();

        // texture coordinates
        const Vec2f &st0 = sts[triIndex * 3];
        const Vec2f &st1 = sts[triIndex * 3 + 1];
        const Vec2f &st2 = sts[triIndex * 3 + 2];
        hitTextureCoordinates = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;
    }
    // member variables
    uint32_t numTris;                        //number of triangles
    std::unique_ptr<Vec3f []> P;             //triangles vertex position
    std::unique_ptr<uint32_t []> trisIndex;  //vertex index array
    std::unique_ptr<Vec3f []> N;             //triangles vertex normals
    std::unique_ptr<Vec2f []> sts;           //triangles texture coordinates
    bool smoothShading = true;               //smooth shading by default
};

TriangleMesh* loadPolyMeshFromFile(const char *file, const Matrix44f &o2w)
{
    std::ifstream ifs;
    try {
        ifs.open(file);
        if (ifs.fail()) throw;
        std::stringstream ss;
        ss << ifs.rdbuf();
        uint32_t numFaces;
        ss >> numFaces;
        std::unique_ptr<uint32_t []> faceIndex(new uint32_t[numFaces]);
        uint32_t vertsIndexArraySize = 0;
        // reading face index array
        for (uint32_t i = 0; i < numFaces; ++i) {
            ss >> faceIndex[i];
            vertsIndexArraySize += faceIndex[i];
        }
        std::unique_ptr<uint32_t []> vertsIndex(new uint32_t[vertsIndexArraySize]);
        uint32_t vertsArraySize = 0;
        // reading vertex index array
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> vertsIndex[i];
            if (vertsIndex[i] > vertsArraySize) vertsArraySize = vertsIndex[i];
        }
        vertsArraySize += 1;
        // reading vertices
        std::unique_ptr<Vec3f []> verts(new Vec3f[vertsArraySize]);
        for (uint32_t i = 0; i < vertsArraySize; ++i) {
            ss >> verts[i].x >> verts[i].y >> verts[i].z;
        }
        // reading normals
        std::unique_ptr<Vec3f []> normals(new Vec3f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> normals[i].x >> normals[i].y >> normals[i].z;
        }
        // reading st coordinates
        std::unique_ptr<Vec2f []> st(new Vec2f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> st[i].x >> st[i].y;
        }

        return new TriangleMesh(o2w, numFaces, faceIndex, vertsIndex, verts, normals, st);
    }
    catch (...) {
        ifs.close();
    }
    ifs.close();

    return nullptr;
}
class Light
{
public:
    Light(const Matrix44f &l2w, const Vec3f &c = 1, const float &i = 1) : lightToWorld(l2w), color(c), intensity(i) {}
    virtual ~Light() {}
    virtual void illuminate(const Vec3f &P, Vec3f &, Vec3f &, float &) const = 0;
    Vec3f color;
    float intensity;
    Matrix44f lightToWorld;
};
class DistantLight : public Light
{
    Vec3f dir;
public:
    DistantLight(const Matrix44f &l2w, const Vec3f &c = 1, const float &i = 1) : Light(l2w, c, i)
    {
        l2w.multDirMatrix(Vec3f(0, 0, -1), dir);
        dir.normalize();  //in case the matrix scales the light
    }
    void illuminate(const Vec3f &P, Vec3f &lightDir, Vec3f &lightIntensity, float &distance) const
    {
        lightDir = dir;
        lightIntensity = color * intensity;
        distance = kInfinity;
    }
};
class PointLight : public Light
{
    Vec3f pos;
public:
    PointLight(const Matrix44f &l2w, const Vec3f &c = 1, const float &i = 1) : Light(l2w, c, i)
    { l2w.multVecMatrix(Vec3f(0), pos); }
    // P: is the shaded point
    void illuminate(const Vec3f &P, Vec3f &lightDir, Vec3f &lightIntensity, float &distance) const
    {
        lightDir = (P - pos);
        float r2 = lightDir.norm();
        distance = sqrt(r2);
        lightDir.x /= distance, lightDir.y /= distance, lightDir.z /= distance;
        // avoid division by 0
        lightIntensity = color * (intensity / (4 * M_PI * r2));
    }
};

enum RayType { kPrimaryRay, kShadowRay };

struct IsectInfo
{
    const Object *hitObject = nullptr;
    float tNear = kInfinity;
    Vec2f uv;
    uint32_t index = 0;
};

bool trace(
        const Vec3f &orig, const Vec3f &dir,
        const std::vector<std::unique_ptr<Object>> &objects,
        IsectInfo &isect,
        RayType rayType = kPrimaryRay)
{
    isect.hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNear = kInfinity;
        uint32_t index = 0;
        Vec2f uv;
        if (objects[k]->intersect(orig, dir, tNear, index, uv) && tNear < isect.tNear) {
            if (rayType == kShadowRay && objects[k]->type == kReflectionAndRefraction) continue;
            isect.hitObject = objects[k].get();
            isect.tNear = tNear;
            isect.index = index;
            isect.uv = uv;
        }
    }

    return (isect.hitObject != nullptr);
}
Vec3f reflect(const Vec3f &I, const Vec3f &N)
{
    return I - 2 * I.dotProduct(N) * N;
}

Vec3f refract(const Vec3f &I, const Vec3f &N, const float &ior)
{
    float cosi = clamp(-1, 1, I.dotProduct(N));
    float etai = 1, etat = ior;
    Vec3f n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
}
void fresnel(const Vec3f &I, const Vec3f &N, const float &ior, float &kr)
{
    float cosi = clamp(-1, 1, I.dotProduct(N));
    float etai = 1, etat = ior;
    if (cosi > 0) { std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        kr = 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        kr = (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

inline float modulo(const float &f)
{
    return f - std::floor(f);
}

Vec3f castRay(
        const Vec3f &orig, const Vec3f &dir,
        const std::vector<std::unique_ptr<Object>> &objects,
        const std::vector<std::unique_ptr<Light>> &lights,
        const Options &options,
        const uint32_t & depth = 0)
{
    if (depth > options.maxDepth) return options.backgroundColor;
    Vec3f hitColor = 0;
    IsectInfo isect;
    if (trace(orig, dir, objects, isect)) {
        Vec3f hitPoint = orig + dir * isect.tNear;
        Vec3f hitNormal;
        Vec2f hitTexCoordinates;
        isect.hitObject->getSurfaceProperties(hitPoint, dir, isect.index, isect.uv, hitNormal, hitTexCoordinates);
        switch (isect.hitObject->type) {
            case kDiffuse:
            {
                for (uint32_t i = 0; i < lights.size(); ++i) {
                    Vec3f lightDir, lightIntensity;
                    IsectInfo isectShad;
                    lights[i]->illuminate(hitPoint, lightDir, lightIntensity, isectShad.tNear);
                    bool vis = !trace(hitPoint + hitNormal * options.bias, -lightDir, objects, isectShad, kShadowRay);
                    // compute the pattern
                    float angle = deg2rad(45);
                    float s = hitTexCoordinates.x * cos(angle) - hitTexCoordinates.y * sin(angle);
                    float t = hitTexCoordinates.y * cos(angle) + hitTexCoordinates.x * sin(angle);
                    float scaleS = 20, scaleT = 20;
                    float pattern = (cos(hitTexCoordinates.y * 2 * M_PI * scaleT) * sin(hitTexCoordinates.x * 2 * M_PI * scaleS) + 1) * 0.5; // isect.hitObject->albedo
                    // float pattern = (modulo(s * scaleS) < 0.5) ^ (modulo(t * scaleT) < 0.5);
                    // float pattern = (modulo(s * scaleS) < 0.5);
                    hitColor = hitColor + vis * lightIntensity * std::max(0.f, hitNormal.dotProduct(-lightDir));
                }
                break;
            }
            case kReflection:
            {
                Vec3f R = reflect(dir, hitNormal);
                R.normalize();
                break;
            }
            case kReflectionAndRefraction:
            {
                Vec3f refractionColor = 0, reflectionColor = 0;
                // compute fresnel
                float kr;
                fresnel(dir, hitNormal, isect.hitObject->ior, kr);
                bool outside = dir.dotProduct(hitNormal) < 0;
                Vec3f bias = options.bias * hitNormal;
                // compute refraction if it is not a case of total internal reflection
                if (kr < 1) {
                    Vec3f refractionDirection = refract(dir, hitNormal, isect.hitObject->ior).normalize();
                    Vec3f refractionRayOrig = outside ? hitPoint - bias : hitPoint + bias;
                    refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, depth + 1);
                }

                Vec3f reflectionDirection = reflect(dir, hitNormal).normalize();
                Vec3f reflectionRayOrig = outside ? hitPoint + bias : hitPoint - bias;
                reflectionColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1);

                // mix the two
                hitColor = hitColor + reflectionColor * kr + refractionColor * (1 - kr);
                break;
            }
            default:
                break;
        }
    }
    else {
        hitColor = options.backgroundColor;
    }

    return hitColor;
}
void render(
        const Options &options,
        const std::vector<std::unique_ptr<Object>> &objects,
        const std::vector<std::unique_ptr<Light>> &lights)
{
    std::unique_ptr<Vec3f []> framebuffer(new Vec3f[options.width * options.height]);
    Vec3f *pix = framebuffer.get();
    float scale = tan(deg2rad(options.fov * 0.5));
    float imageAspectRatio = options.width / (float)options.height;
    Vec3f orig;
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
    auto timeStart = std::chrono::high_resolution_clock::now();
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {
            // generate primary ray direction
            float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();
            *(pix++) = castRay(orig, dir, objects, lights, options);
        }
        fprintf(stderr, "\r%3d%c", uint32_t(j / (float)options.height * 100), '%');
    }
    auto timeEnd = std::chrono::high_resolution_clock::now();
    auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
    fprintf(stderr, "\rDone: %.2f (sec)\n", passedTime / 1000);

    // save framebuffer to file
    std::ofstream ofs;
    ofs.open("out.ppm");
    ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
    for (uint32_t i = 0; i < options.height * options.width; ++i) {
        char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
        char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
        char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
        ofs << r << g << b;
    }
    ofs.close();
}
int main(int argc, char **argv)
{
    // loading gemetry
    std::vector<std::unique_ptr<Object>> objects;
    // lights
    std::vector<std::unique_ptr<Light>> lights;
    Options options;

#if 1
    // glass and pen example
    // setting up options
    options.fov = 36.87;
    options.maxDepth = 10;
    options.bias = 0.001;
    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f(-0.972776, 0, -0.231748, 0, -0.114956, 0.8683, 0.482536, 0, 0.201227, 0.49604, -0.844661, 0, 6.696465, 22.721296, -30.097976, 1);

    TriangleMesh *mesh1 = loadPolyMeshFromFile("./geometry/backdrop.geo", Matrix44f());
    if (mesh1 != nullptr) {
        mesh1->type = kDiffuse;
        objects.push_back(std::unique_ptr<Object>(mesh1));
    }

    TriangleMesh *mesh3 = loadPolyMeshFromFile("./geometry/cylinder.geo", Matrix44f());
    if (mesh3 != nullptr) {
        mesh3->type = kReflectionAndRefraction;
        mesh3->ior = 1.5;
        objects.push_back(std::unique_ptr<Object>(mesh3));
    }

    TriangleMesh *mesh4 = loadPolyMeshFromFile("./geometry/pen.geo", Matrix44f());
    if (mesh4 != nullptr) {
        mesh4->type = kDiffuse;
        mesh4->albedo = 0.18;
        mesh4->smoothShading = false;
        objects.push_back(std::unique_ptr<Object>(mesh4));
    }

    Matrix44f xform1;
    xform1[3][0] = -1.2;
    xform1[3][1] = 6;
    xform1[3][2] = -3;
    Sphere *sph1 = new Sphere(xform1, 5);
    sph1->type = kReflectionAndRefraction;

    Matrix44f l2w(11.146836, -5.781569, -0.0605886, 0, -1.902827, -3.543982, -11.895445, 0, 5.459804, 10.568624, -4.02205, 0, 0, 0, 0, 1);

#elif 0
    // simple plane example (patterns)
    options.fov = 36.87;
    options.width = 1024;
    options.height = 747;
    options.bias = 1e-4;
    options.cameraToWorld = Matrix44f(0.707107, 0, -0.707107, 0, -0.331295, 0.883452, -0.331295, 0, 0.624695, 0.468521, 0.624695, 0, 28, 21, 28, 1);

    TriangleMesh *mesh = loadPolyMeshFromFile("./geometry/plane.geo", Matrix44f());
    if (mesh != nullptr) {
        mesh->type = kDiffuse;
        mesh->albedo = 0.18;
        mesh->smoothShading = false;
        objects.push_back(std::unique_ptr<Object>(mesh));
    }
    Matrix44f l2w(11.146836, -5.781569, -0.0605886, 0, -1.902827, -3.543982, -11.895445, 0, 5.459804, 10.568624, -4.02205, 0, 0, 0, 0, 1);
#elif 0
    // multiple glasses example
    options.fov = 36.87;
    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f(0.999945, 0, 0.0104718, 0, 0.00104703, 0.994989, -0.0999803, 0, -0.0104193, 0.0999858, 0.994934, 0, -0.978596, 17.911879, 75.483369, 1);

    TriangleMesh *mesh = loadPolyMeshFromFile("./geometry/glasses.geo", Matrix44f());
    if (mesh != nullptr) {
        mesh->type = kReflectionAndRefraction;
        mesh->ior = 1.3;
        objects.push_back(std::unique_ptr<Object>(mesh));
    }


    TriangleMesh *mesh1 = loadPolyMeshFromFile("./geometry/backdrop1.geo", Matrix44f());
    if (mesh1 != nullptr) {
        mesh1->type = kDiffuse;
        mesh1->albedo = 0.18;
        objects.push_back(std::unique_ptr<Object>(mesh1));
    }

    Matrix44f l2w(0.95292, 0.289503, 0.0901785, 0, -0.0960954, 0.5704, -0.815727, 0, -0.287593, 0.768656, 0.571365, 0, 0, 0, 0, 1);
#else
    // aliasing example
    options.fov = 36.87;
    options.width = 1024;
    options.height = 747;
    options.cameraToWorld = Matrix44f(0.999945, 0, 0.0104718, 0, 0.00104703, 0.994989, -0.0999803, 0, -0.0104193, 0.0999858, 0.994934, 0, -0.978596, 17.911879, 75.483369, 1);

    Matrix44f xform;
    xform[0][0] = 10;
    xform[1][1] = 10;
    xform[2][2] = 10;
    xform[3][2] = -40;
    TriangleMesh *mesh = loadPolyMeshFromFile("./geometry/plane.geo", xform);
    if (mesh != nullptr) {
        mesh->type = kDiffuse;
        mesh->albedo = 0.18;
        mesh->smoothShading = false;
        objects.push_back(std::unique_ptr<Object>(mesh));
    }
    Matrix44f l2w(11.146836, -5.781569, -0.0605886, 0, -1.902827, -3.543982, -11.895445, 0, 5.459804, 10.568624, -4.02205, 0, 0, 0, 0, 1);
#endif
    lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 1)));

    // finally, render
    render(options, objects, lights);

    return 0;
}