#pragma once

#include "defines.h"

class Ray {
public:
    Ray();
    Ray(const point3& o, const vec3& d, float tmin_=0, float tmax_=100000, int depth_=0);
    ~Ray();

    point3 orig; //! start point of the ray
    vec3 dir; //! ray direction, normalized

    float tmax; //! if no intersection computed, should be a large value, else contains the t of the intersection, so do not need to considere intersection if the founded t is more than tmax
    float tmin; //! usefull to store entry point in a aabb
    int depth; //! number of reflection/refraction

    int sign[3]; //! sign of the x,y,z component of dir, 0 -> positive, 1->negative. To optimize aabb intersection
    vec3 invdir; //! =1/dir, optimize aabb

    bool shadow = false; // is it a shadow ray? to optimize kdtree traversal
                         
    point3 dox;
    point3 doy;
    vec3 ddx;
    vec3 ddy;

    bool hasDifferentials = true;
private:
};
