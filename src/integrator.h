#include "defines.h"
#include "raytracer.h"
#include "Light.h"
#include "sampling/sampling.h"
#include "Sky.h"
#include "LightDistribution.h"

// from PBR
color3 directIllumination(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection, Light* light, const vec2& uScattering, const vec2& uLight);

class Integrator {
  public:
    Integrator(Sampler* sampler) : m_sampler(sampler) {}
    virtual void preprocess(Scene* scene, Sampler* sampler) {}
    virtual color3 trace_ray( Scene* scene,
                              Ray* ray,
                              KdTree* tree,
                              Intersection* intersection,
                              Sampler* sampler ) = 0;

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

    const int maxDepth;
    std::unique_ptr<LightDistribution> m_lightDistrib;
};
