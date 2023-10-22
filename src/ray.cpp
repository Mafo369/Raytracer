#include "ray.h"

Ray::Ray(){}

Ray::Ray(const point3& o, const vec3& d, const vec2& pixel_, float tmin_, float tmax_, int depth_) :
    orig(o),
    dir(d),
    tmax(tmax_),
    tmin(tmin_),
    depth(depth_),
    invdir(1.f/d),
    pixel(pixel_)
{
    sign[0] = d.x>=0?0:1;
    sign[1] = d.y>=0?0:1;
    sign[2] = d.z>=0?0:1;
}

Ray::~Ray(){}
