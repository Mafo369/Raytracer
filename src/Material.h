#pragma once

#include "defines.h"
#include <vector>
#include "raytracer.h"

class Material {
  public:
    Material() = default;
    virtual ~Material() = default;
    virtual color3 shade(Intersection* intersection, vec3 v, color3 lc, float intensity, std::vector<vec3> &samples) = 0;
    virtual color3 textureColor(float u, float v, int face) = 0;
    virtual color3 ambientColor(color3 lightColor) = 0;
    virtual color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) = 0;
};
