//
// Created by 胡栋月 on 7/12/22.
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

#include "geometry.h"

const float kInfinity = std::numeric_limits<float>::max();
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

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
    uint32_t width;
    uint32_t height;
    float fov;
    Matrix44f cameraToWorld;
};

class Object
{
public:
    Object() : color(dis(gen), dis(gen), dis(gen)) {}
    virtual ~Object() {}
    // Method to compute the intersection of the object with a ray
    // Returns true if an intersection was found, false otherwise
    // See method implementation in children class for details
    virtual bool intersect(const Vec3f &, const Vec3f &, float &) const = 0;
    // Method to compute the surface data such as normal and texture coordnates at the intersection point.
    // See method implementation in children class for details
    virtual void getSurfaceData(const Vec3f &, Vec3f &, Vec2f &) const = 0;
    Vec3f color;
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
    Sphere(const Vec3f &c, const float &r) : radius(r), radius2(r *r ), center(c) {}
    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t) const
    {
        float t0, t1;  //solutions for t if the ray intersects
#if 0
        // geometric solution
        Vec3f L = center - orig;
        float tca = L.dotProduct(dir);
        if (tca < 0) return false;
        float d2 = L.dotProduct(L) - tca * tca;
        if (d2 > radius2) return false;
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;
#else
        // analytic solution
        Vec3f L = orig - center;
        float a = dir.dotProduct(dir);
        float b = 2 * dir.dotProduct(L);
        float c = L.dotProduct(L) - radius2;
        if (!solveQuadratic(a, b, c, t0, t1)) return false;
#endif
        if (t0 > t1) std::swap(t0, t1);

        if (t0 < 0) {
            t0 = t1;  //if t0 is negative, let's use t1 instead
            if (t0 < 0) return false;  //both t0 and t1 are negative
        }

        t = t0;

        return true;
    }
    void getSurfaceData(const Vec3f &Phit, Vec3f &Nhit, Vec2f &tex) const
    {
        Nhit = Phit - center;
        Nhit.normalize();
        // In this particular case, the normal is simular to a point on a unit sphere
        // centred around the origin. We can thus use the normal coordinates to compute
        // the spherical coordinates of Phit.
        // atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
        // acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
    float radius, radius2;
    Vec3f center;
};

bool trace(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects, float &tNear, const Object *&hitObject)
{
    tNear = kInfinity;
    std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
    for (; iter != objects.end(); ++iter) {
        float t = kInfinity;
        if ((*iter)->intersect(orig, dir, t) && t < tNear) {
            hitObject = iter->get();
            tNear = t;
        }
    }

    return (hitObject != nullptr);
}

Vec3f castRay(
        const Vec3f &orig, const Vec3f &dir,
        const std::vector<std::unique_ptr<Object>> &objects)
{
    Vec3f hitColor = 0;
    const Object *hitObject = nullptr;  //this is a pointer to the hit object
    float t;  //this is the intersection distance from the ray origin to the hit point
    if (trace(orig, dir, objects, t, hitObject)) {
        Vec3f Phit = orig + dir * t;
        Vec3f Nhit;
        Vec2f tex;
        hitObject->getSurfaceData(Phit, Nhit, tex);
        // Use the normal and texture coordinates to shade the hit point.
        // The normal is used to compute a simple facing ratio and the texture coordinate
        // to compute a basic checker board pattern
        float scale = 4;
        float pattern = (fmodf(tex.x * scale, 1) > 0.5) ^ (fmodf(tex.y * scale, 1) > 0.5);
        hitColor = std::max(0.f, Nhit.dotProduct(-dir)) * mix(hitObject->color, hitObject->color * 0.8, pattern);
    }

    return hitColor;
}

void render(
        const Options &options,
        const std::vector<std::unique_ptr<Object>> &objects)
{
    Vec3f *framebuffer = new Vec3f[options.width * options.height];
    Vec3f *pix = framebuffer;
    float scale = tan(deg2rad(options.fov * 0.5));
    float imageAspectRatio = options.width / (float)options.height;
    Vec3f orig;
    options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
    for (uint32_t j = 0; j < options.height; ++j) {
        for (uint32_t i = 0; i < options.width; ++i) {
#ifdef MAYA_STYLE
            float x = (2 * (i + 0.5) / (float)options.width - 1) * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale * 1 / imageAspectRatio;
#elif

            float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
#endif
            Vec3f dir;
            options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
            dir.normalize();
            *(pix++) = castRay(orig, dir, objects);
        }
    }

    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./out.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
    for (uint32_t i = 0; i < options.height * options.width; ++i) {
        char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
        char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
        char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
        ofs << r << g << b;
    }

    ofs.close();

    delete [] framebuffer;
}

int main(int argc, char **argv)
{
    // creating the scene (adding objects and lights)
    std::vector<std::unique_ptr<Object>> objects;

    // generate a scene made of random spheres
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

    // finally, render
    render(options, objects);

    return 0;
}
