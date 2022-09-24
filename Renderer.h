#pragma once

#include "Walnut/Image.h"

#include "Camera.h"
#include "ray.h"

#include <memory>
#include <glm/glm.hpp>

#include "defines.h"
#include "image.h"
#include "ray.h"
#include "raytracer.h"
#include "scene.h"
#include "textures.hpp"
#include "kdtree.h"
#include "scene_types.h"
#include "intersection.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


class Renderer
{
public:
  Renderer();

	void OnResize(uint32_t width, uint32_t height);
	void Render(const CameraI& camera);

	std::shared_ptr<Walnut::Image> GetFinalImage() const { return m_FinalImage; }
private:
	glm::vec4 TraceRay(const Ray& ray);
private:
	std::shared_ptr<Walnut::Image> m_FinalImage;
	uint32_t* m_ImageData = nullptr;

  KdTree* tree;
  Image* img;
  Scene* scene; 
};
