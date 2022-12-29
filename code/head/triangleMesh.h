//
// Created by 胡栋月 on 26/12/22.
//

#ifndef MYGL_TRIANGLEMESH_H
#define MYGL_TRIANGLEMESH_H

#include "geometry.h"
#include "object.h"
#include "tools.h"
#include <cstdio>
#include <iostream>
#include <cstdlib>

bool intersectTriangle(const Vec3f &orig, const Vec3f &dir, const Vec3f &vertex0, const Vec3f &vertex1, const Vec3f &vertex2, float &t, float &u, float &v) {
#ifdef METHOD_GEOMETRY // geometry solution
    std::cout << "method geometry" << std::endl;
    Vec3f v0v1 = vertex1 - vertex0;
    Vec3f v0v2 = vertex2 - vertex0;

    Vec3f normal = v0v1.crossProduct(v0v2);
    normal.normalize();

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

    return (t > 0) ? true : false;

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
    if(fabs(det) < kEpsilon) return false;
#endif
    Vec3f T = orig - vertex0;

    float inverseDet = 1 / det;

    u = inverseDet * P.dotProduct(T);
    if(u < 0) return false;

    Vec3f Q = T.crossProduct(v0v1);

    v = inverseDet * Q.dotProduct(dir);
    if(v < 0 || u+v > 1) return false;

    t = inverseDet * Q.dotProduct(v0v2);

    return (t > 0) ? true : false;
#endif
}

class Triangle : public Object {
private:
    Vec3f vertex0, vertex1, vertex2;
public:
    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, const Matrix44f &o2w = Matrix44f(), const float &al = 0.18,
             const MaterialType &materialType = kDiffuse,
             const char *name = "TriangleMesh") : Object(al, o2w, materialType, name) {
        objectToWorld.multDirMatrix(v0, vertex0);
        objectToWorld.multDirMatrix(v1, vertex1);
        objectToWorld.multDirMatrix(v2, vertex2);
    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &t, uint32_t &index, Vec2f &uv) const {
        float u, v;
        if(intersectTriangle(orig, dir, vertex0, vertex1, vertex2, t, u, v)) {
            uv.x = u;
            uv.y = v;
            return true;
        }
        return false;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
        Vec3f v0v1 = vertex1 - vertex0;
        Vec3f v0v2 = vertex2 - vertex0;

        Vec3f normal = v0v1.crossProduct(v0v2);
        normal.normalize();
        Nhit = normal;
        tex.x = (1 + atan2(Nhit.z, Nhit.x) / M_PI) * 0.5;
        tex.y = acosf(Nhit.y) / M_PI;
    }
};

class TriangleMesh : public Object {
private:

    std::unique_ptr<Vec3f []> P; // vertices position in theis Triangle Mesh
    std::unique_ptr<uint32_t []> triangleVerticesIndex; // vertex index array
    std::unique_ptr<Vec3f []> N; // normal for triangles
    std::unique_ptr<Vec2f []> texCoordinates; // texture for triangles
public:
    uint32_t numTris; // number of triangles in this Triangle Mesh
    TriangleMesh(const uint32_t numFaces, // the input mesh may not be triangulated
                 const std::unique_ptr<uint32_t []> &faceIndex,
                 const std::unique_ptr<uint32_t []> &verticesIndex,
                 const std::unique_ptr<Vec3f []> &vertices,
                 std::unique_ptr<Vec3f []> &normals,
                 std::unique_ptr<Vec2f []> &uv, const Matrix44f &o2w = Matrix44f(), const float &al = 0.18,
                 const MaterialType &materialType = kDiffuse,
                 const char *name = "TriangleMesh") :
            numTris(0), Object(al, o2w, materialType, name) {
        uint32_t accumulateIndex = 0, maxVerticesIndex = 0;

        // find out the total num of triangles in this mesh
        for(uint32_t i=0; i<numFaces; ++i) {
            // record how many triangles are needed. the polygon with 5 vertices, need 5 - 2 = 3 triangles to represent
            numTris += faceIndex[i] - 2;
            for(uint32_t j=0; j<faceIndex[i]; ++j) {
                // find the largest index for vertices, then allocate memory to store them
                if(verticesIndex[accumulateIndex+j] > maxVerticesIndex)
                    maxVerticesIndex = verticesIndex[accumulateIndex+j];
            }
            accumulateIndex += faceIndex[i];
        }maxVerticesIndex += 1;

        // allocate memory to store this vertices
        P = std::unique_ptr<Vec3f []>(new Vec3f[maxVerticesIndex]);
        for(uint32_t i=0; i<maxVerticesIndex; ++i) {
            objectToWorld.multDirMatrix(vertices[i], P[i]);
        }

        // used to transform normals
        Matrix44f transformNormals = worldToObject.transpose();

        // allocate memory to store triangle indices
        triangleVerticesIndex = std::unique_ptr<uint32_t []>(new uint32_t[numTris*3]);
        uint32_t l = 0;
        N = std::unique_ptr<Vec3f []>(new Vec3f[numTris*3]);
        texCoordinates = std::unique_ptr<Vec2f []>(new Vec2f[numTris*3]);
        for(uint32_t i=0, accumulateIndex=0; i<numFaces; ++i) {
            for(uint32_t j=0; j<faceIndex[i] - 2; ++j) {
                triangleVerticesIndex[l] = verticesIndex[accumulateIndex];
                triangleVerticesIndex[l + 1] = verticesIndex[accumulateIndex + j + 1];
                triangleVerticesIndex[l + 2] = verticesIndex[accumulateIndex + j + 2];
                transformNormals.multDirMatrix(normals[accumulateIndex], N[l]);
                transformNormals.multDirMatrix(normals[accumulateIndex + j + 1], N[l + 1]);
                transformNormals.multDirMatrix(normals[accumulateIndex + j + 2], N[l + 2]);
                texCoordinates[l] = uv[accumulateIndex];
                texCoordinates[l + 1] = uv[accumulateIndex + j + 1];
                texCoordinates[l + 2] = uv[accumulateIndex + j + 2];
                l += 3;
            }
            accumulateIndex += faceIndex[i];
        }

    }

    bool intersect(const Vec3f &orig, const Vec3f &dir, float &tNear, uint32_t &index, Vec2f &uv) const {
        uint32_t  j = 0;
        bool isIntersect = false;
        for(uint32_t i=0; i<numTris; ++i) {
            Vec3f &v0 = P[triangleVerticesIndex[j]];
            Vec3f &v1 = P[triangleVerticesIndex[j + 1]];
            Vec3f &v2 = P[triangleVerticesIndex[j + 2]];
            float u, v, t = kInfinity;
            if(intersectTriangle(orig, dir, v0, v1, v2, t, u, v) && t < tNear) {
                tNear = t;
                uv.x = u;
                uv.y = v;
                isIntersect = true;
                index = i;
            }
            j += 3;
        }
        return isIntersect;
    }

    void getSurfaceData(const Vec3f &Phit, const Vec3f &dir, const uint32_t &index, const Vec2f &uv, Vec3f &Nhit, Vec2f &tex) const {
#ifdef SMOOTH_SHADING
        // vertex normal
        const Vec3f &n0 = N[index * 3];
        const Vec3f &n1 = N[index * 3 + 1];
        const Vec3f &n2 = N[index * 3 + 2];
        Nhit = (1 - uv.x - uv.y) * n0 + uv.x * n1 + uv.y * n2;
#else
        // face normal
        const Vec3f &v0 = P[triangleVerticesIndex[index * 3]];
        const Vec3f &v1 = P[triangleVerticesIndex[index * 3 + 1]];
        const Vec3f &v2 = P[triangleVerticesIndex[index * 3 + 2]];
        Nhit = (v1 - v0).crossProduct(v2 - v0);
#endif

        Nhit.normalize();

        // texture coordinates
        const Vec2f  &st0 = texCoordinates[index * 3];
        const Vec2f  &st1 = texCoordinates[index * 3 + 1];
        const Vec2f  &st2 = texCoordinates[index * 3 + 2];
        tex = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;
    }
};

TriangleMesh* generatePolyShphere(float rad, uint32_t divs)
{
    // generate points
    uint32_t numVertices = (divs - 1) * divs + 2;
    std::unique_ptr<Vec3f []> P(new Vec3f[numVertices]);
    std::unique_ptr<Vec3f []> N(new Vec3f[numVertices]);
    std::unique_ptr<Vec2f []> st(new Vec2f[numVertices]);

    float u = -M_PI_2;
    float v = -M_PI;
    float du = M_PI / divs;
    float dv = 2 * M_PI / divs;

    P[0] = N[0] = Vec3f(0, -rad, 0);
    uint32_t k = 1;
    for (uint32_t i = 0; i < divs - 1; i++) {
        u += du;
        v = -M_PI;
        for (uint32_t j = 0; j < divs; j++) {
            float x = rad * cos(u) * cos(v);
            float y = rad * sin(u);
            float z = rad * cos(u) * sin(v) ;
            P[k] = N[k] = Vec3f(x, y, z);
            st[k].x = u / M_PI + 0.5;
            st[k].y = v * 0.5 / M_PI + 0.5;
            v += dv, k++;
        }
    }
    P[k] = N[k] = Vec3f(0, rad, 0);

    uint32_t npolys = divs * divs;
    uint32_t normalSize = (6 + (divs - 1) * 4) * divs;
    std::unique_ptr<uint32_t []> faceIndex(new uint32_t[npolys]);
    std::unique_ptr<uint32_t []> vertsIndex(new uint32_t[normalSize]);

    // create the connectivity lists
    uint32_t vid = 1, numV = 0, l = 0;
    k = 0;
    for (uint32_t i = 0; i < divs; i++) {
        for (uint32_t j = 0; j < divs; j++) {
            if (i == 0) {
                faceIndex[k++] = 3;
                vertsIndex[l] = 0;
                vertsIndex[l + 1] = j + vid;
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid : j + vid + 1;
                l += 3;
            }
            else if (i == (divs - 1)) {
                faceIndex[k++] = 3;
                vertsIndex[l] = j + vid + 1 - divs;
                vertsIndex[l + 1] = vid + 1;
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
                l += 3;
            }
            else {
                faceIndex[k++] = 4;
                vertsIndex[l] = j + vid + 1 - divs;
                vertsIndex[l + 1] = j + vid + 1;
                vertsIndex[l + 2] = (j == (divs - 1)) ? vid + 1 : j + vid + 2;
                vertsIndex[l + 3] = (j == (divs - 1)) ? vid + 1 - divs : j + vid + 2 - divs;
                l += 4;
            }
            numV++;
        }
        vid = numV;
    }

    // unify with the situation in loadPolyMeshFromFile
    std::unique_ptr<Vec3f []> normal(new Vec3f[normalSize]);
    std::unique_ptr<Vec2f []> stExpand(new Vec2f[normalSize]);

    uint32_t accumulateIndex = 0;
    for(uint32_t i=0; i<npolys; ++i) {
        for(uint32_t j=0; j<faceIndex[i]; ++j) {
            normal[accumulateIndex + j] = N[vertsIndex[accumulateIndex + j]];
            stExpand[accumulateIndex + j] = st[vertsIndex[accumulateIndex + j]];
        }
        accumulateIndex += faceIndex[i];
    }

//    return new TriangleMesh(npolys, faceIndex, vertsIndex, P, N, st);
    return new TriangleMesh(npolys, faceIndex, vertsIndex, P, normal, stExpand);
}

TriangleMesh* loadPolyMeshFromFile(const char *file, const Matrix44f &objectToWorld = Matrix44f())
{
    std::ifstream ifs;
    try {
        ifs.open(file);
        if (ifs.fail()) throw;
        std::stringstream ss;
        ss << ifs.rdbuf();
        uint32_t numFaces;
        ss >> numFaces;
        std::unique_ptr<uint32_t []> faceIndex(new uint32_t[numFaces]);
        uint32_t vertsIndexArraySize = 0;
        // reading face index array
        for (uint32_t i = 0; i < numFaces; ++i) {
            ss >> faceIndex[i];
            vertsIndexArraySize += faceIndex[i];
        }
        std::unique_ptr<uint32_t []> vertsIndex(new uint32_t[vertsIndexArraySize]);
        uint32_t vertsArraySize = 0;
        // reading vertex index array
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> vertsIndex[i];
            if (vertsIndex[i] > vertsArraySize) vertsArraySize = vertsIndex[i];
        }
        vertsArraySize += 1;
        // reading vertices
        std::unique_ptr<Vec3f []> verts(new Vec3f[vertsArraySize]);
        for (uint32_t i = 0; i < vertsArraySize; ++i) {
            ss >> verts[i].x >> verts[i].y >> verts[i].z;
        }
        // reading normals
        std::unique_ptr<Vec3f []> normals(new Vec3f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> normals[i].x >> normals[i].y >> normals[i].z;
        }
        // reading st coordinates
        std::unique_ptr<Vec2f []> st(new Vec2f[vertsIndexArraySize]);
        for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
            ss >> st[i].x >> st[i].y;
        }

        return new TriangleMesh(numFaces, faceIndex, vertsIndex, verts, normals, st, objectToWorld);
    }
    catch (...) {
        ifs.close();
    }
    ifs.close();

    return nullptr;
}


#endif //MYGL_TRIANGLEMESH_H
