//
// Created by 胡栋月 on 26/12/22.
//

#ifndef MYGL_TOOLS_H
#define MYGL_TOOLS_H

#include <random>

#include "geometry.h"

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

inline float modulo(const float &f)
{
    return f - std::floor(f);
}

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

#endif //MYGL_TOOLS_H
