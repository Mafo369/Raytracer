#pragma once
#include "../Object.h"

class Sphere : public Object {
  public:
    Sphere(std::shared_ptr<Material> mat, Transform transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const override;
    Intersection sample(const Intersection& inter, const point2& u) const override;
    Intersection sample(const point2& u) const override;
    float pdf(const Intersection& inter, const vec3& wi) const override;
};
