//
// Created by 胡栋月 on 26/12/22.
//

#ifndef MYGL_OBJECT_H
#define MYGL_OBJECT_H

#include "geometry.h"
#include "tools.h"

// options in command
struct Options {
    uint32_t width = 640;
    uint32_t height = 480;
    Vec3f backgroundColor = kDefaultBackgroundColor;
    float fov = 90;
    Matrix44f cameraToWorld;
    float bias = 1e-4;
    uint32_t maxDepth = 5;
    bool geometrySolution = false;
    bool culling = false;
    bool smoothShading = true;
};

enum MaterialType { kDiffuse, kReflection, kReflectionAndRefraction };

enum RayType { kPrimaryRay, kShadowRay };

// the virtual class for supported object
class Object {
public:
    Object(const float &al, const Matrix44f &o2w, const MaterialType &type, const char *na)
    : color(dis(gen), dis(gen), dis(gen)),
    albedo(al), objectToWorld(o2w), worldToObject(o2w.inverse()),
    materialType(type), name(na) {}
    virtual ~Object() {}

    // weather intersect and the distance of hit point
    virtual bool intersect(const Vec3f &orig, const Vec3f &dir, const Options &options, float &t, uint32_t &index, Vec2f &uv) const = 0;

    // get the surface date for texture in the hit point
    virtual void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, const Options &options, Vec3f &Nhit, Vec2f &tex) const = 0;
    Vec3f color;
    Matrix44f objectToWorld, worldToObject;
    float albedo;
    const char *name;
    MaterialType materialType;
    float ior;
};

// intersection Info : hitObject, tNear, uv, index
struct IntersecInfo {
    const Object *hitObject = nullptr;
    float tNear = kInfinity;
    Vec2f uv;
    uint32_t index = 0;
};

#endif //MYGL_OBJECT_H
