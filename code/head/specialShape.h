//
// Created by 胡栋月 on 26/12/22.
//

#ifndef MYGL_SPECIALSHAPE_H
#define MYGL_SPECIALSHAPE_H

#include "geometry.h"
#include "object.h"
#include "tools.h"

/**
 *  support classes : Sphere
 *
 *  remove: Plane, Disk, Box
 */

class Sphere : public Object {
private:
    float radius, radius2;
    Vec3f center;
public:
    Sphere(const Vec3f &c, const float &r, const Matrix44f &o2w = Matrix44f(), const float &al = 0.18,
           const MaterialType &materialType = kReflection,
           const char *name = "Sphere") :
           Object (al, o2w, materialType, name) {
        radius = r;
        radius2 = r * r;
        objectToWorld.multVecMatrix(c, center);
    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, const Options &options, float &t, uint32_t &index, Vec2f &uv) const {
        uv = 0;
        float t0, t1;
        if(options.geometrySolution) { /// geometric solution
            Vec3f  L = center - orig;
            float tca = L.dotProduct(dir);
            if(tca < 0) return false;
            float d2 = L.dotProduct(L) - tca * tca;
            if(d2 > radius2) return false;
            float thc = sqrt(radius2 - d2);
            t0 = tca - thc;
            t1 = tca + thc;
        } else { /// analytic solution
            Vec3f L = orig - center;
            float a = dir.dotProduct(dir);
            float b = 2 * dir.dotProduct(L);
            float c = L.dotProduct(L) - radius2;
            if(!solveQuadratic(a, b, c, t0, t1)) return false;
        }
        if(t0 > t1) std::swap(t0, t1);

        if(t0 < 0) {
            t0 = t1;
            if(t0 < 0) return false;
        }

        t = t0;

        return true;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, const Options &options, Vec3f &Nhit, Vec2f &tex) const {
        Nhit = Phit - center;
        Nhit.normalize();
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};

/**
 * wrong
 * not just normal, but distance
 * Ax + By + Cz + d = 0;
 */
//class Plane : public Object {
//private:
//    Vec3f normal;
//public:
//    Plane(const Vec3f &N, const Matrix44f &o2w = Matrix44f(), const float &al = 0.18,
//          const MaterialType &materialType = kReflection,
//          const char *name = "Sphere") :
//            Object (al, o2w, materialType, name) {
//        Matrix44f transformNormals = worldToObject.transpose();
//        transformNormals.multDirMatrix(N, normal);
//    }
//
//    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &index, Vec2f &uv) const {
//        uv = 0;
//        float det = dir.dotProduct(normal);
//#ifdef CULLING
//        if(det < kEpsilon) return false;
//#else
//        if(abs(det) < kEpsilon) return false;
//#endif
//        return true;
//    }
//
//    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
//        Nhit = normal;
//        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
//        tex.y = acosf(Nhit.y) / M_PI;
//    }
//};


//class Disk : public Object {
//private:
//    Vec3f normal;
//    Vec3f center;
//    float radius, radius2;
//public:
//    Disk(const Vec3f &N, const Vec3f &c, const float &r, const Matrix44f &o2w = Matrix44f(), const float &al = 0.18,
//         const MaterialType &materialType = kDiffuse,
//         const char *name = "TriangleMesh") :
//         normal(N), center(c), radius(r), radius2(r*r),
//         Object(al, o2w, materialType, name) {
//        radius = r;
//        radius2 = r * r;
//        objectToWorld.multVecMatrix(c, center);
//        Matrix44f transformNormal = worldToObject.transpose();
//        transformNormal.multDirMatrix(N, normal);
//    }
//
//    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &index, Vec2f &uv) const {
//        uv = 0;
//        float det = dir.dotProduct(normal);
//#ifdef CULLING
//        if(det < kEpsilon) return false;
//#else
//        if(abs(det) < kEpsilon) return false;
//#endif
//        // compute the t
//        float D = -(normal.dotProduct(center));
//        t = -(normal.dotProduct(orig) + D) / det;
//
//        Vec3f Phit = orig + t*dir;
//        Vec3f pc = center - Phit;
//
//        if(pc.dotProduct(pc) > radius2) return false;
//        return true;
//    }
//
//    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
//        Nhit = normal;
//        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
//        tex.y = acosf(Nhit.y) / M_PI;
//    }
//};
//
//class Box : public Object {
//private:
//
//public:
//
//};

#endif //MYGL_SPECIALSHAPE_H
