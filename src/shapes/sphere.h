#pragma once
#include "../Object.h"

class Sphere : public Object {
  public:
    Sphere(std::shared_ptr<Material> mat, Transform transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const;
};
