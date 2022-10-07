#pragma once
#include "../Object.h"

class Triangle : public Object {
  public:
    Triangle(Material mat, glm::mat4 transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const;
};
