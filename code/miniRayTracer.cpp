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

# include "geometry.h"

#if !defined METHOD_GEOMETRY || !defined CULLING
//#define METHOD_GEOMETRY
//#define CULLING

/*
 *
 */
// get the max float as infinity
static const float kInfinity = std::numeric_limits<float>::max();

static constexpr float kEpsilon = 1e-8;

static const Vec3f kDefaultBackgroundColor = Vec3f(0.2, 0.7, 0.8);

// get random
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

// degree to radian
inline
float deg2rad(const float &deg) {
    return deg * M_PI / 180;
}

// limit the range
inline
float clamp(const float &lo, const float &hi, const float &v) {
    return std::max(lo, std::min(hi, v));
}

//
inline
Vec3f mix(const Vec3f &a, const Vec3f& b, const float &mixValue){
    return a * (1 - mixValue) + b * mixValue;
}

// options in command
struct Options {
    uint32_t width;
    uint32_t height;
    Vec3f backgroundColor;
    float fov;
    Matrix44f cameraToWorld;
};

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1) {
    int delta = b * b - 4 * a * c;
    if(delta < 0) return false;
    else if(delta == 0) { x0 = x1 = -0.5 * b / a; }
    else {
        float q = (b > 0) ?
                  -0.5 * (b + sqrt(delta)) :
                  -0.5 * (b - sqrt(delta));
        x0 = q / a;
        x1 = c / q;
    }
    return true;
}

//intersect(const Vec3f &orig, const Vec3f &dir, float &tNear, uint32_t &triIndex, Vec2f &uv)

// the virtual class for supported object
class Object {
public:
    Object() : color(dis(gen), dis(gen), dis(gen)) {}
    virtual ~Object() {}

    // weather intersect and the distance of hit point
    virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;

    // get the surface date for texture in the hit point
    virtual void getSurfaceData(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec3f &, Vec2f &) const = 0;
    Vec3f color;
};

/**
 *  support classes : Sphere, Triangle, Plane, Disk, Box
 *
 *  todo list: World-to-Object Matrix in Ray-Tracing

 */
class Sphere : public Object {
private:
    float radius, radius2;
    Vec3f center;
public:
    Sphere(const Vec3f &c, const float &r) : center(c), radius(r), radius2(r*r) {}

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &index, Vec2f &uv) const {
        uv = 0;
        float t0, t1;
#ifdef METHOD_GEOMETRY // geometric solution
        Vec3f  L = center - orig;
        float tca = L.dotProduct(dir);
        if(tca < 0) return false;
        float d2 = L.dotProduct(L) - tca * tca;
        if(d2 > radius2) return false;
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;
#else // analytic solution
        Vec3f L = orig - center;
        float a = dir.dotProduct(dir);
        float b = 2 * dir.dotProduct(L);
        float c = L.dotProduct(L) - radius2;
        if(!solveQuadratic(a, b, c, t0, t1)) return false;
#endif
        if(t0 > t1) std::swap(t0, t1);

        if(t0 < 0) {
            t0 = t1;
            if(t0 < 0) return false;
        }

        t = t0;

        return true;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
        Nhit = Phit - center;
        Nhit.normalize();
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};

bool intersectTriangle(const Vec3f &orig, const Vec3f &dir, const Vec3f &vertex0, const Vec3f &vertex1, const Vec3f &vertex2, float &t, float &u, float &v) {
#ifdef METHOD_GEOMETRY // geometry solution
    Vec3f v0v1 = vertex1 - vertex0;
    Vec3f v0v2 = vertex2 - vertex0;

    float denom = normal.dotProduct(normal);

    // exclude parallel
    float normalDotRayDirection = normal.dotProduct(dir);
    if(abs(normalDotRayDirection) < kEpsilon)
        return false;

    // compute t for intersection
    float D = -(normal.dotProduct(vertex0));
    t = -(normal.dotProduct(orig) + D) / normalDotRayDirection;

    // exclue t < 0, intersection behind the origin of the ray
    if(t < 0) return false;

    // get the hit point
    Vec3f Phit = orig + t*dir;

    Vec3f v0P = Phit - vertex0;
    Vec3f v1P = Phit - vertex1;
    Vec3f v2P = Phit - vertex2;

    // judge hit point inside triangle or not
    Vec3f J;
    J = v0v1.crossProduct(v0P);
    if(normal.dotProduct(J) < 0) return false;

    Vec3f v1v2 = vertex2 - vertex1;
    J = v1v2.crossProduct(v1P);
    if((u = normal.dotProduct(J))< 0) return false;

    Vec3f v2v0 = vertex0 - vertex2;
    J = (-v0v2).crossProduct(v2P);
    if((v = normal.dotProduct(J)) < 0) return false;

    u /= denom;
    v /= denom;

    uv.x = u;
    uv.y = v;

    return true;

#else // MOLLER_TRUMBORE
    Vec3f v0v1 = vertex1 - vertex0;
    Vec3f v0v2 = vertex2 - vertex0;
    // v0v1 ✖️ v0v2 -> normal
    // v0v2 ✖️ v0v1 -> -normal
    // dir . (v0v1 ✖️ v0v2) = -dir . (v0v2 ✖️ v0v1) = -(dir ✖️ v0v2) . v0v1
    // P = (dir ✖️ v0v2)
    Vec3f P = dir.crossProduct(v0v2);

    float det = v0v1.dotProduct(P);
#ifdef CULLING //  det < 0 should be culled
    if(det < kEpsilon) return false;
#else
    if(abs(det) < kEpsilon) return false;
#endif
    Vec3f T = orig - vertex0;

    float inverseDet = 1 / det;

    u = inverseDet * P.dotProduct(T);
    if(u < 0) return false;

    Vec3f Q = T.crossProduct(v0v1);

    v = inverseDet * Q.dotProduct(dir);
    if(v < 0 || u+v > 1) return false;

    t = inverseDet * Q.dotProduct(v0v2);

    return true;
#endif
}

class Triangle : public Object {
private:
    Vec3f vertex0, vertex1, vertex2;
    Vec3f normal;
public:
    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2) : vertex0(v0), vertex1(v1), vertex2(v2) {
        Vec3f v0v1 = vertex1 - vertex0;
        Vec3f v0v2 = vertex2 - vertex0;
        normal = v0v1.crossProduct(v0v2);
        normal.normalize();
    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &index, Vec2f &uv) const {
        float u, v;
        if(intersectTriangle(orig, dir, vertex0, vertex1, vertex2, t, u, v)) {
            uv.x = u;
            uv.y = v;
            return true;
        }
        return false;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
        Nhit = normal;
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};

class TriangleMesh : public Object {
private:

    std::unique_ptr<Vec3f []> P; // vertices position in theis Triangle Mesh
    std::unique_ptr<uint32_t []> triangleVerticesIndex; // vertex index array
    std::unique_ptr<Vec3f []> N; // normal for triangles
    std::unique_ptr<Vec2f []> texCoordinates; // texture for triangles
public:
    uint32_t numTris; // number of triangles in this Triangle Mesh
    TriangleMesh(const uint32_t numFaces, // the input mesh may not be triangulated
                 const std::unique_ptr<uint32_t []> &faceIndex,
                 const std::unique_ptr<uint32_t []> &verticesIndex,
                 const std::unique_ptr<Vec3f []> &vertices,
                 std::unique_ptr<Vec3f []> &normals,
                 std::unique_ptr<Vec2f []> &uv) : numTris(0) {
        uint32_t accumulateIndex = 0, maxVerticesIndex = 0;

        // find out the total num of triangles in this mesh
        for(uint32_t i=0; i<numFaces; ++i) {
            // record how many triangles are needed. the polygon with 5 vertices, need 5 - 2 = 3 triangles to represent
            numTris += faceIndex[i] - 2;
            for(uint32_t j=0; j<faceIndex[i]; ++j) {
                // find the largest index for vertices, then allocate memory to store them
                if(verticesIndex[accumulateIndex+j] > maxVerticesIndex)
                    maxVerticesIndex = verticesIndex[accumulateIndex+j];
            }
            accumulateIndex += faceIndex[i];
        }
        maxVerticesIndex += 1;

        // allocate memory to store this vertices
        P = std::unique_ptr<Vec3f []>(new Vec3f[maxVerticesIndex]);
        for(uint32_t i=0; i<maxVerticesIndex; ++i) {
            P[i] = vertices[i];
        }

        // allocate memory to store triangle indices
        triangleVerticesIndex = std::unique_ptr<uint32_t []>(new uint32_t[numTris*3]);
        uint32_t l = 0;
        N = std::unique_ptr<Vec3f []>(new Vec3f[numTris*3]);
        texCoordinates = std::unique_ptr<Vec2f []>(new Vec2f[numTris*3]);
        for(uint32_t i=0, accumulateIndex=0; i<numFaces; ++i) {
            for(uint32_t j=0; j<faceIndex[i] - 2; ++j) {
                triangleVerticesIndex[l] = verticesIndex[accumulateIndex];
                triangleVerticesIndex[l + 1] = verticesIndex[accumulateIndex + j + 1];
                triangleVerticesIndex[l + 2] = verticesIndex[accumulateIndex + j + 2];
                N[l] = normals[accumulateIndex];
                N[l + 1] = normals[accumulateIndex + j + 1];
                N[l + 2] = normals[accumulateIndex + j + 2];
                texCoordinates[l] = uv[accumulateIndex];
                texCoordinates[l + 1] = uv[accumulateIndex + j + 1];
                texCoordinates[l + 2] = uv[accumulateIndex + j + 2];
                l += 3;
            }
            accumulateIndex += faceIndex[i];
        }

    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tNear, uint32_t &index, Vec2f &uv) const {
        uint32_t  j = 0;
        bool isIntersect = false;
        for(uint32_t i=0; i<numTris; ++i) {
            Vec3f &v0 = P[triangleVerticesIndex[j]];
            Vec3f &v1 = P[triangleVerticesIndex[j + 1]];
            Vec3f &v2 = P[triangleVerticesIndex[j + 2]];
            float u, v, t = kInfinity;
            if(intersectTriangle(orig, dir, v0, v1, v2, t, u, v) && t < tNear) {
                tNear = t;
                uv.x = u;
                uv.y = v;
                isIntersect = true;
                index = i;
            }
            j += 3;
        }
        return isIntersect;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
        // face normal
        const Vec3f &v0 = P[triangleVerticesIndex[index * 3]];
        const Vec3f &v1 = P[triangleVerticesIndex[index * 3 + 1]];
        const Vec3f &v2 = P[triangleVerticesIndex[index * 3 + 2]];

        Nhit = (v1 - v0).crossProduct(v2 - v0);
        Nhit.normalize();

        // texture coordinates
        const Vec2f  &st0 = texCoordinates[index * 3];
        const Vec2f  &st1 = texCoordinates[index * 3 + 1];
        const Vec2f  &st2 = texCoordinates[index * 3 + 2];
        tex = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;

        // vertex normal
        const Vec3f &n0 = N[index * 3];
        const Vec3f &n1 = N[index * 3 + 1];
        const Vec3f &n2 = N[index * 3 + 2];
        Nhit = (1 - uv.x - uv.y) * n0 + uv.x * n1 + uv.y * n2;
    }
};

class Plane : public Object {
private:
    Vec3f normal;
public:
    Plane(const Vec3f &N) : normal(N) {}

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &index, Vec2f &uv) const {
        uv = 0;
        float det = dir.dotProduct(normal);
#ifdef CULLING
        if(det < kEpsilon) return false;
#else
        if(abs(det) < kEpsilon) return false;
#endif
        return true;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
        Nhit = normal;
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};


class Disk : public Object {
private:
    Vec3f normal;
    Vec3f center;
    float radius, radius2;
public:
    Disk(const Vec3f &N, const Vec3f &c, const float &r) : normal(N), center(c), radius(r), radius2(r*r) {}

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &index, Vec2f &uv) const {
        uv = 0;
        float det = dir.dotProduct(normal);
#ifdef CULLING
        if(det < kEpsilon) return false;
#else
        if(abs(det) < kEpsilon) return false;
#endif
        // compute the t
        float D = -(normal.dotProduct(center));
        t = -(normal.dotProduct(orig) + D) / det;

        Vec3f Phit = orig + t*dir;
        Vec3f pc = center - Phit;

        if(pc.dotProduct(pc) > radius2) return false;
        return true;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
        Nhit = normal;
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};

class Box : public Object {
private:

public:

};


/**
 * Test each object for intersection
 * @param orig
 * @param dir
 * @param objects
 * @param tNear
 * @param hitObject
 * @return
 */
bool trace(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects, float &tNear, uint32_t &index, Vec2f &uv, const Object *&hitObject) {
    tNear = kInfinity;
    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    for(; iter != objects.end(); ++iter) {
        float t = kInfinity;
        uint32_t indexForFace;
        Vec2f uvForFace;

        if((*iter)->intersect(orig, dir, t, indexForFace, uvForFace) && t < tNear) {
            hitObject = iter->get();
            tNear = t;
            index = indexForFace;
            uv = uvForFace;
        }
    }

    return (hitObject != nullptr);
}

/**
 * Cast a ray with specific origin and direction, then compute the color in hit point.
 * @param orig
 * @param dir
 * @param objects
 * @return
 */
Vec3f castRay(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects, const Options &options) {
    Vec3f hitColor = options.backgroundColor;
    const Object *hitObject = nullptr;
    float t;

    uint32_t index = 0;
    Vec2f uv;
    if(trace(orig, dir, objects, t, index, uv, hitObject)) {
        Vec3f Phit = orig + dir * t;
        Vec3f Nhit;
        Vec2f tex;
        hitObject->getSurfaceData(Phit, dir, index, uv,Nhit, tex);

        float scale = 4;
        float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y * scale, 1) > 0.5);
        hitColor = std::max(0.f, Nhit.dotProduct(-dir)) * mix(hitObject->color, hitObject->color * 0.8, pattern);
    }

    return hitColor;
}

/**
 * Render the whole scene with options and objects info
 * @param options
 * @param objects
 */
void render(const Options &options, const std::vector<std::unique_ptr<Object>> &objects) {
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
            *(pix++) = castRay(orig, dir, objects, options);
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

TriangleMesh* generatePolyShphere(float rad, uint32_t divs)
{
    // generate points
    uint32_t numVertices = (divs - 1) * divs + 2;
    std::unique_ptr<Vec3f []> P(new Vec3f[numVertices]);
    std::unique_ptr<Vec3f []> N(new Vec3f[numVertices]);
    std::unique_ptr<Vec2f []> st(new Vec2f[numVertices]);

    float u = -M_PI_2;
    float v = -M_PI;
    float du = M_PI / divs;
    float dv = 2 * M_PI / divs;

    P[0] = N[0] = Vec3f(0, -rad, 0);
    uint32_t k = 1;
    for (uint32_t i = 0; i < divs - 1; i++) {
        u += du;
        v = -M_PI;
        for (uint32_t j = 0; j < divs; j++) {
            float x = rad * cos(u) * cos(v);
            float y = rad * sin(u);
            float z = rad * cos(u) * sin(v) ;
            P[k] = N[k] = Vec3f(x, y, z);
            st[k].x = u / M_PI + 0.5;
            st[k].y = v * 0.5 / M_PI + 0.5;
            v += dv, k++;
        }
    }
    P[k] = N[k] = Vec3f(0, rad, 0);

    uint32_t npolys = divs * divs;
    std::unique_ptr<uint32_t []> faceIndex(new uint32_t[npolys]);
    std::unique_ptr<uint32_t []> vertsIndex(new uint32_t[(6 + (divs - 1) * 4) * divs]);

    // create the connectivity lists
    uint32_t vid = 1, numV = 0, l = 0;
    k = 0;
    for (uint32_t i = 0; i < divs; i++) {
        for (uint32_t j = 0; j < divs; j++) {
            if (i == 0) {
                faceIndex[k++] = 3;
                vertsIndex[l] = 0;
                vertsIndex[l + 1] = j + vid;
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid : j + vid + 1;
                l += 3;
            }
            else if (i == (divs - 1)) {
                faceIndex[k++] = 3;
                vertsIndex[l] = j + vid + 1 - divs;
                vertsIndex[l + 1] = vid + 1;
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
                l += 3;
            }
            else {
                faceIndex[k++] = 4;
                vertsIndex[l] = j + vid + 1 - divs;
                vertsIndex[l + 1] = j + vid + 1;
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 : j + vid + 2;
                vertsIndex[l + 3] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
                l += 4;
            }
            numV++;
        }
        vid = numV;
    }

    return new TriangleMesh(npolys, faceIndex, vertsIndex, P, N, st);
}

TriangleMesh* loadPolyMeshFromFile(const char *file)
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

        return new TriangleMesh(numFaces, faceIndex, vertsIndex, verts, normals, st);
    }
    catch (...) {
        ifs.close();
    }
    ifs.close();

    return nullptr;
}

int main(int argc, char **argv) {
    std::vector<std::unique_ptr<Object>> objects;

    uint32_t numSpheres = 3;
    gen.seed(0);


    // setting up options
    Options options;
    options.width = 640;
    options.height = 480;
//    options.width = 2560;
//    options.height = 1600;
    options.backgroundColor = kDefaultBackgroundColor;
    options.fov = 51.52;

#if 0

    options.cameraToWorld = Matrix44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

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

#else

    Matrix44f tmp = Matrix44f(0.707107, -0.331295, 0.624695, 0, 0, 0.883452, 0.468521, 0, -0.707107, -0.331295, 0.624695, 0, -1.63871, -5.747777, -40.400412, 1);
    options.cameraToWorld = tmp.inverse();

    TriangleMesh *mesh = loadPolyMeshFromFile("./resource/cow.geo");
    if (mesh != nullptr) objects.push_back(std::unique_ptr<Object>(mesh));

# endif

    render(options, objects);

    return 0;
}

#endif