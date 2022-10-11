#pragma once

#include "../Material.h"

class Blinn : public Material {
  public:
    Blinn();
    ~Blinn() override;

    color3 shade(Intersection *intersection, vec3 v, Light* light, float intensity) override;
    color3 textureColor(float u, float v, int face) override;
    color3 ambientColor(Intersection* intersection, color3 lightColor) override;
    color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) override;

    color3 m_diffuseColor;
    color3 m_specularColor;
    float m_shininess;
    vec3 m_reflection;
    vec3 m_refraction;
    vec3 m_absorption;
    float m_IOR;

    texture* m_texture = nullptr;

};
