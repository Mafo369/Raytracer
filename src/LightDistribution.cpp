#include "LightDistribution.h"

LightDistribution::~LightDistribution() {}

UniformLightDistribution::UniformLightDistribution(Scene *scene) {
    std::vector<float> prob(scene->lights.size(), 1.f);
    distrib.reset(new Distribution1D(&prob[0], int(prob.size())));
}

const Distribution1D *UniformLightDistribution::lookup(const point3& p) const {
    return distrib.get();
}
