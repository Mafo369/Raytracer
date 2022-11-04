#pragma once

#include "../Material.h"
#include <random>

class Blinn : public Material {
  public:
    Blinn();
    ~Blinn() override;

    color3 f(const vec3& wo, const vec3& wi);
    color3 sample_f(vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type) override;
    color3 shade(Intersection *intersection, vec3 v, Light* light, float intensity) override;
    color3 textureColor(float u, float v, int face) override;
    color3 ambientColor(Ray* ray, Intersection* intersection, color3 lightColor) override;
    color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) override;

    color3 refractionColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection, color3 reflectionShade, vec3 normal);
    color3 reflectionColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection, color3& color, vec3 normal);

    color3 m_diffuseColor;
    color3 m_specularColor;
    float m_shininess;
    vec3 m_reflection;
    vec3 m_refraction;
    vec3 m_absorption;
    float m_reflectionGloss;
    float m_refractionGloss;
};
