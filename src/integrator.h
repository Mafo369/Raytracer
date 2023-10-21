#pragma once

#include "defines.h"
#include "Light.h"
#include "sampling/sampling.h"
#include "Sky.h"
#include "LightDistribution.h"

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection, bool show_lights=true, Sampler* sampler=nullptr);

// from PBR
color3 directIllumination(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection, Light* light, const vec2& uScattering, const vec2& uLight);

class Integrator {
  public:
    Integrator(Sampler* sampler, int depth) : maxDepth(depth), m_sampler(sampler) {}
    virtual void preprocess(Scene* scene, Sampler* sampler) {}
    virtual color3 trace_ray( Scene* scene,
                              Ray* ray,
                              KdTree* tree,
                              Intersection* intersection,
                              Sampler* sampler ) = 0;
  int maxDepth;

  Sampler* m_sampler;
};

class Pathtracer : public Integrator {
  public:
    Pathtracer(int maxDepth, Sampler* sampler);
    void preprocess(Scene* scene, Sampler* sampler) override;
    color3 trace_ray( Scene* scene,
                      Ray* ray,
                      KdTree* tree,
                      Intersection* intersection,
                      Sampler* sampler ) override;

    std::unique_ptr<LightDistribution> m_lightDistrib;
};
