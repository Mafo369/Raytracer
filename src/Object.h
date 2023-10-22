#pragma once

#include "defines.h"
#include "ray.h"
#include "intersection.h"
#include "geometry.h"
#include "textures.hpp"
#include "materials/CookTorrance.h"

class Object {
  public:
    Object(std::shared_ptr<Material> material, Transform transform);
    virtual ~Object();
    virtual bool intersect(Ray* ray, Intersection* intersection) const = 0;
    virtual Intersection sample(const Intersection& inter, const point2& u, float* pdf) const = 0;
    virtual Intersection sample(const point2& u) const = 0;
    virtual float pdf(const Intersection& inter, const vec3& wi) const;
    virtual float Area() const { return 1; }

    void transformRay(vec3& origin, vec3& direction) const;

    Transform transform;
    Geometry geom;
    std::shared_ptr<Material> mat;
};

