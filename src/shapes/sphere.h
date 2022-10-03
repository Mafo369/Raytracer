#pragma once
#include "../Object.h"

class Sphere : public Object {
  public:
    Sphere(Material mat, glm::mat4 transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const;
};
