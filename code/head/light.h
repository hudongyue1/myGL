//
// Created by 胡栋月 on 23/12/22.
//

#ifndef MYGL_LIGHT_H
#define MYGL_LIGHT_H

#include "geometry.h"
#include "tools.h"

class Light {
public:
    Light(const Matrix44f &l2w) : lightToWorld(l2w) {}
    virtual ~Light() {}
    virtual void getDirectionAndIntensity(const Vec3f &Phit, Vec3f &lightDirection, Vec3f &lightIntensity, float &distance) const = 0;
    Matrix44f lightToWorld;
    Vec3f color;
    float intensity;
};

class DistantLight : public Light {
public:
    DistantLight(const Matrix44f &l2w = Matrix44f(), const Vec3f &c = 1, const float &i = 10) : Light(l2w) {
        this->color = c;
        this->intensity = i;
        l2w.multDirMatrix(Vec3f(0, 0, -1), dir);
        dir.normalize();
    }
    void getDirectionAndIntensity(const Vec3f &Phit, Vec3f &lightDirection, Vec3f &lightIntensity, float &distance) const {
        lightDirection = dir;
        lightIntensity = color * intensity;
        distance = kInfinity;
    }
    Vec3f dir;
};

class PointLight : public Light {
public:
    PointLight(const Matrix44f&l2w = Matrix44f(), const Vec3f &c = 1, const float  &i = 1000) : Light(l2w) {
        this->color = c;
        this->intensity = i;
        l2w.multDirMatrix(Vec3f(0), pos);
    }

    void getDirectionAndIntensity(const Vec3f &Phit, Vec3f &lightDirection, Vec3f &lightIntensity, float &distance) const {
        lightDirection = Phit - pos;
        float r2 = lightDirection.norm();
        distance = sqrt(r2);
        lightDirection.normalize();
        lightIntensity = color * (intensity / (4 * M_PI * r2));
    }
    Vec3f pos;
};

#endif //MYGL_LIGHT_H
