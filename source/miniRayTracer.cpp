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

/*
 *
 */
// get the max float as infinity
const float kInfinity = std::numeric_limits<float>::max();

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
    float fov;
    Matrix44f cameraToWorld;
};

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

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1) {

}

class Sphere : public Object {
private:
    float radius, radius2;
    Vec3f center;
public:
    Sphere(const Vec3f &c, const float &r) : center(c), radius(r), radius2(r*r) {}

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const {
        float t0, t1;

        // geometric solution
        Vec3f  L = center - orig;
        float tca = L.dotProduct(dir);
        if(tca < 0) return false;
        float d2 = L.dotProduct(L) - tca * tca;
        if(d2 > radius2) return false;
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;

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

    std::ofstream ofs("../out/out.ppm", std::ios::out | std::ios::binary);
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

    uint32_t numSpheres = 32;
    gen.seed(0);
    for (uint32_t i = 0; i < numSpheres; ++i) {
        Vec3f randPos((0.5 - dis(gen)) * 10, (0.5 - dis(gen)) * 10, (0.5 + dis(gen) * 10));
        float randRadius = (0.5 + dis(gen) * 0.5);
        objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius)));
    }

    // setting up options
    Options options;
    options.width = 640;
    options.height = 480;
    options.fov = 51.52;
    options.cameraToWorld = Matrix44f(0.945519, 0, -0.325569, 0, -0.179534, 0.834209, -0.521403, 0, 0.271593, 0.551447, 0.78876, 0, 4.208271, 8.374532, 17.932925, 1);

    render(options, objects);

    return 0;
}
