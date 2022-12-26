//
// Created by 胡栋月 on 26/12/22.
//

#ifndef MYGL_OBJECT_H
#define MYGL_OBJECT_H

#include "geometry.h"
#include "tools.h"

// the virtual class for supported object
class Object {
public:
    Object(const float &al, const Matrix44f &o2w) : color(dis(gen), dis(gen), dis(gen)), albedo(al), objectToWorld(o2w), worldToObject(o2w.inverse()) {}
    virtual ~Object() {}

    // weather intersect and the distance of hit point
    virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;

    // get the surface date for texture in the hit point
    virtual void getSurfaceData(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec3f &, Vec2f &) const = 0;
    Vec3f color;
    Matrix44f objectToWorld, worldToObject;
    float albedo;
};


#endif //MYGL_OBJECT_H
