#ifndef __RAY_H__
#define __RAY_H__

#include "defines.h"
// RAY
typedef struct ray_s {

    point3 orig; //! start point of the ray
    vec3 dir; //! ray direction, normalized

    float tmax; //! if no intersection computed, should be a large value, else contains the t of the intersection, so do not need to considere intersection if the founded t is more than tmax
    float tmin; //! usefull to store entry point in a aabb
    int depth; //! number of reflection/refraction

    int sign[3]; //! sign of the x,y,z component of dir, 0 -> positive, 1->negative. To optimize aabb intersection
    vec3 invdir; //! =1/dir, optimize aabb

} Ray;

inline void rayInit(Ray *r, point3 o, vec3 d, float tmin=0, float tmax=100000, int depth=0) {
    r->orig = o;
    r->dir = d;
    r->tmin = tmin;
    r->tmax = tmax;
    r->depth = depth;
    r->sign[0] = r->dir.x>=0?0:1;
    r->sign[1] = r->dir.y>=0?0:1;
    r->sign[2] = r->dir.z>=0?0:1;
    r->invdir = 1.f/d;
}

inline point3 rayAt(const Ray r, float t) {
    return r.orig + t*r.dir;
}

#endif
