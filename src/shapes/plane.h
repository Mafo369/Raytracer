#pragma once
#include "../Object.h"

class Plane : public Object {
  public:
    Plane(std::shared_ptr<Material> mat, glm::mat4 transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const;
};
