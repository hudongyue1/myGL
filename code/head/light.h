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
    PointLight(const Matrix44f&l2w = Matrix44f(), const Vec3f &c = 1, const float  &i = 20) : Light(l2w) {
        this->color = c;
        this->intensity = i;
        l2w.multDirMatrix(Vec3f(0), pos);
    }

    void getDirectionAndIntensity(const Vec3f &Phit, Vec3f &lightDirection, float &lightIntensity) const {
        lightDirection = Phit - pos;
        float r2 = lightDirection.norm();
        lightDirection.normalize();
        lightIntensity = intensity / (4 * M_PI * r2);
    }
    Vec3f pos;
};

#endif //MYGL_LIGHT_H
