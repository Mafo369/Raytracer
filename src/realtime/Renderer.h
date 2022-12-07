#pragma once

#include "Walnut/Image.h"

#include "Camera.h"
#include "../ray.h"

#include <memory>
#include <glm/glm.hpp>

#include "../defines.h"
#include "../image.h"
#include "../raytracer.h"
#include "../scene.h"
#include "../textures.hpp"
#include "../kdtree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class Integrator;

class Renderer
{
public:
  struct Settings {
    bool accumulate = true;
  };

  Renderer();

	void OnResize(uint32_t width, uint32_t height);
	void Render(const CameraI& camera);

	std::shared_ptr<Walnut::Image> GetFinalImage() const { return m_FinalImage; }

  void setScene(Scene* renderScene){ scene = renderScene; 
                                      tree = initKdTree(scene);}
  void setIntegrator(Integrator* integrator) { m_integrator = integrator; }
  void setSampler(Sampler* s) { sampler = s; }
  Integrator* getIntegrator() { return m_integrator; }

  Scene* scene; 

  void ResetFrameIndex() { m_frameIndex = 1; }
  Settings& getSettings() { return m_settings; }

private:
	glm::vec4 TraceRay(const Ray& ray);
private:
	std::shared_ptr<Walnut::Image> m_FinalImage;
	uint32_t* m_ImageData = nullptr;

  Settings m_settings;

  KdTree* tree;
  RenderImage* img;
  glm::vec4* m_accumulationData = nullptr;
  uint32_t m_frameIndex = 1;

  Sampler* sampler;
  Integrator* m_integrator;
};
