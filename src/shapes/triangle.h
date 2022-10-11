#pragma once
#include "../Object.h"

class Triangle : public Object {
  public:
    Triangle(std::shared_ptr<Material> mat, Transform transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const;
};
