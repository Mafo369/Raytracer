#pragma once

#include "../Material.h"

class CookTorrance : public Material {
  public:
    CookTorrance(bool transparent);
    ~CookTorrance() override;

    color3 sample_f(vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type) override;
    color3 shade(Intersection *intersection, vec3 v, Light* light, float intensity) override;
    color3 textureColor(float u, float v, int face) override;
    color3 ambientColor(Ray* ray, Intersection* intersection, color3 lightColor) override;
    color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) override;

    float m_roughness; //! 0.001 - 0.01 : very smooth finish with slight imperfections. 0.1 : relatively rough. 0.3-0.7 extremely rough 
    color3 m_specularColor;	//! Specular "albedo"
    color3 m_diffuseColor;	//! Base color
    bool m_transparent;
};
