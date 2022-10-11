#pragma once
#include "../Object.h"

class Plane : public Object {
  public:
    Plane(std::shared_ptr<Material> mat, Transform transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const;
};
