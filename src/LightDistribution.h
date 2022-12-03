#pragma once

#include "scene.h"
#include "sampling/sampling.h"

// From PBR
class LightDistribution {
  public:
    virtual ~LightDistribution();

    virtual const Distribution1D* lookup(const point3& p) const = 0;
};

class UniformLightDistribution : public LightDistribution {
  public:
    UniformLightDistribution(Scene* scene);
    const Distribution1D* lookup(const point3& p) const;

  private:
    std::unique_ptr<Distribution1D> distrib;
};
