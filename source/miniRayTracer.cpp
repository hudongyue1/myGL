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

// the virtual class for supported object
class Object {
public:
    Object() : color(dis(gen), dis(gen), dis(gen)) {}
    virtual ~Object() {}

    // weather intersect and the distance of hit point
    virtual bool intersect(const Vec3f &, const Vec3f &, float &) const = 0;

    // get the surface date for texture in the hit point
    virtual void getSurfaceData(const Vec3f &, Vec3f &, Vec2f &) const = 0;
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

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const {
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

    void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const {
        Nhit = Phit - center;
        Nhit.normalize();
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};

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

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const {
        float u, v;
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

    void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const {
        Nhit = normal;
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};


class Plane : public Object {
private:
    Vec3f normal;
public:
    Plane(const Vec3f &N) : normal(N) {}

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const {
        float det = dir.dotProduct(normal);
#ifdef CULLING
        if(det < kEpsilon) return false;
#else
        if(abs(det) < kEpsilon) return false;
#endif
        return true;
    }

    void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const {
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

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const {
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

    void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const {
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
bool trace(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects, float &tNear, const Object *&hitObject) {
    tNear = kInfinity;
    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    for(; iter != objects.end(); ++iter) {
        float t = kInfinity;
        if((*iter)->intersect(orig, dir, t) && t < tNear) {
            hitObject = iter->get();
            tNear = t;
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
Vec3f castRay(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects) {
    Vec3f hitColor = 0;
    const Object *hitObject = nullptr;
    float t;
    if(trace(orig, dir, objects, t, hitObject)) {
        Vec3f Phit = orig + dir * t;
        Vec3f Nhit;
        Vec2f tex;
        hitObject->getSurfaceData(Phit, Nhit, tex);

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

    for(uint32_t j=0; j<options.height; ++j) {
        for(uint32_t i=0; i<options.width; ++i) {
            float x = (2 * (i + 0.5) / (float) options.width - 1) * scale * imageAspectRatio;
            float y = (1 - 2 * (j + 0.5) / (float) options.height) * scale;

            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();
            *(pix++) = castRay(orig, dir, objects);
        }
    }

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

    uint32_t numSpheres = 3;
    gen.seed(0);

    // sphere
    for (uint32_t i = 0; i < numSpheres; ++i) {
        Vec3f randPos((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, (0.5 + dis(gen) * 10));
        float randRadius = (0.5 + dis(gen) * 0.5);
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

    // setting up options
    Options options;
    options.width = 640;
    options.height = 480;
    options.backgroundColor = kDefaultBackgroundColor;
//    options.width = 2560;
//    options.height = 1600;
    options.fov = 51.52;
    options.cameraToWorld = Matrix44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

    render(options, objects);

    return 0;
}

#endif