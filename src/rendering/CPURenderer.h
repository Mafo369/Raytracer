#pragma once

#include "ImageRenderer.h"
#include "../kdtree.h"
#include "../integrator.h"

class CPURenderer : ImageRenderer {
public:
    void init(std::string& name, glm::vec2 resolution, Scene* scene) override;

    void destroy() override;

    void run() override;

private:
    void scaleDifferentials( Ray* ray, float s );

    RenderImage* m_renderImage;

    KdTree* m_kdTree;
    Sampler* m_sampler;    
    Integrator* m_tracer;
};
