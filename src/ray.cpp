#include "ray.h"

Ray::Ray(){}

Ray::Ray(point3 o, vec3 d, vec2 pixel_, float tmin_, float tmax_, int depth_) {
    orig = o;
    dir = d;
    tmin = tmin_;
    tmax = tmax_;
    depth = depth_;
    sign[0] = dir.x>=0?0:1;
    sign[1] = dir.y>=0?0:1;
    sign[2] = dir.z>=0?0:1;
    invdir = 1.f/d;
    pixel = pixel_;
}

Ray::~Ray(){}
