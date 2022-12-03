#pragma once

#include "../Material.h"

class CookTorrance : public Material {
  public:
    CookTorrance(MatType type = DIFFUSE);
    ~CookTorrance() override;

    color3 f(const vec3& wo, const vec3& wi, const vec3& n) override;
    color3 sample_f(vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type) override;
    float pdf(const vec3& wo, const vec3& wi, const vec3& n) override;
    color3 shade(Intersection *intersection, vec3 v, Light* light, float intensity) override;
    color3 textureColor(float u, float v, int face) override;
    color3 ambientColor(Ray* ray, Intersection* intersection, color3 lightColor) override;
    color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) override;

    color3 eval(Ray* ray, Intersection* intersection, const vec3& wi, float* scatteringPdf) override;
    color3 sample(Ray* ray, Intersection* intersection, const vec2& uScattering, vec3* wi, float* scatteringPdf) override;

    float m_roughness; //! 0.001 - 0.01 : very smooth finish with slight imperfections. 0.1 : relatively rough. 0.3-0.7 extremely rough 
    float m_metalness;
    MatType m_type;
  private:
    color3 scratchAPixelScatter(Ray* ray, Scene* scene, KdTree* tree, Intersection* intersection);
    color3 myScatter(Ray* ray, Scene* scene, KdTree* tree, Intersection* intersection);
};
