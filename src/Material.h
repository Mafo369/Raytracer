#pragma once

#include "defines.h"
#include <vector>
#include "raytracer.h"
#include "Light.h"

class Material {
  public:
    Material() = default;
    virtual ~Material() = default;
    virtual color3 shade(Intersection* intersection, vec3 v, Light* light, float intensity) = 0;
    virtual color3 textureColor(float u, float v, int face) = 0;
    virtual color3 ambientColor(Ray* ray, Intersection* intersection, color3 lightColor) = 0;
    virtual color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) = 0;
    texture* m_texture = nullptr;
    color3 m_emission = color3(0,0,0);
};
