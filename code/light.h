//
// Created by 胡栋月 on 23/12/22.
//

#ifndef MYGL_LIGHT_H
#define MYGL_LIGHT_H

#include "geometry.h"

class Light {
public:
    Light(const Matrix44f &l2w) : lightToWorld(l2w) {}
    virtual ~Light() {}
    Matrix44f lightToWorld;
    Vec3f color;
    float intensity;
};

class DistantLight : public Light {
public:
    DistantLight(const Matrix44f &l2w = Matrix44f(), const Vec3f &c = 1, const float &i = 15) : Light(l2w) {
        this->color = c;
        this->intensity = i;
        l2w.multDirMatrix(Vec3f(0, 0, -1), dir);
        dir.normalize();
    }
    Vec3f dir;
};

class PointLight : Light {
public:
    PointLight(const Matrix44f&l2w = Matrix44f(), const Vec3f &c = 1, const float  &i = 1) : Light(l2w) {
        this->color = c;
        this->intensity = i;
        l2w.multDirMatrix(Vec3f(0), pos);
    }
    Vec3f pos;
};

#endif //MYGL_LIGHT_H
