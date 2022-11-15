#pragma once

#include "defines.h"
#include <vector>
#include "raytracer.h"
//#include "Light.h"
#include "sampling/sampler.h"

class Material {
  public:
    Material() = default;
    virtual ~Material() = default;
    virtual color3 f(const vec3& wo, const vec3& wi, const vec3& n) = 0;
    virtual color3 sample_f(vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type) = 0;
    virtual float pdf(const vec3& wo, const vec3& wi, const vec3& n) = 0;
    virtual color3 shade(Intersection* intersection, vec3 v, Light* light, float intensity) = 0;
    virtual color3 textureColor(float u, float v, int face) = 0;
    virtual color3 ambientColor(Ray* ray, Intersection* intersection, color3 lightColor) = 0;
    virtual color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) = 0;
    texture* m_texture = nullptr;
    color3 m_emission = color3(0,0,0);
    color3 m_diffuseColor = color3(0,0,0);
    color3 m_specularColor = color3(0,0,0);

    float m_IOR;	//! Index of refraction (for dielectric)
};

color3 specularReflect(Ray* ray, Intersection* intersection, Scene* scene, KdTree* tree, Sampler* sampler);
color3 specularTransmission(Ray* ray, Intersection* intersection, Scene* scene, KdTree* tree, Sampler* sampler);
