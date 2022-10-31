
#include <algorithm>

#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "scene.h"
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <iostream>
#include <chrono>

#include <glm/gtc/epsilon.hpp>
#include <random>
#include <omp.h>

#include "textures.hpp"
#include "bsdf.hpp"

#include "Light.h"
#include "Object.h"
#include "Camera.h"

#include "sampling/stratified.h"

// Mostly for debugging purposes
bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection)
{
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();

  //!\todo loop on each object of the scene to compute intersection

  float dist;

  for (size_t i = 0; i < objectCount; i++)
  {
    Intersection temp;
    if(scene->objects[i]->intersect(ray, &temp)){
      float temp_dist = ray->tmax;
      if (hasIntersection)
      {
        if (temp_dist < dist)
        {
          dist = temp_dist;
          *intersection = temp;
        }
      }
      else
      {
        hasIntersection = true;
        *intersection = temp;
        dist = temp_dist;
      }

    }
  }
  return hasIntersection;
}

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection)
{
  color3 ret = color3(0, 0, 0);

  if (ray->depth > 3)
    return color3(0.f);

  if (intersectKdTree(scene, tree, ray, intersection))
  {
    intersection->hit = true;
    // Compute necessary differential information for texture filtering
    intersection->computeDifferentials(ray);

    // if skybox return directly the corresponding color
    if(intersection->face != -1){
      return intersection->mat->textureColor(intersection->u, intersection->v, intersection->face);
    }

    // shading
    size_t lightsCount = scene->lights.size();
    for (size_t i = 0; i < lightsCount; i++)
    {
      vec3 v = ray->dir * -1.0f;
      if(scene->lights[i]->isAmbient()){
        ret += intersection->mat->ambientColor(intersection, scene->lights[i]->getColor());
        continue;
      }
      auto intensity = scene->lights[i]->intensityAt(intersection->position, scene, tree, v, intersection); 
      if(intensity > 0.0f){
        ret += intersection->mat->shade(intersection, v, scene->lights[i], intensity);
      }
    }

    // If max contribution, we stop
    if (ret.r >= 1.f && ret.g >= 1.f && ret.b >= 1.f && ray->depth > 0)
      return color3(1.f);

    // Scatter
    ret += intersection->mat->scatterColor(scene, tree, ray, intersection);
  }
  else
  {
    vec3 pixelUV = vec3((float)ray->pixel.x / scene->cam->imgWidth, (float)ray->pixel.y / scene->cam->imgHeight, 0.f);
    vec3 p = glm::mod(scene->m_skyTexture->m_transform.transformTo(pixelUV), 1.f);
    if(scene->m_skyTexture != nullptr && ray->depth == 0)
      ret = scene->skyColor * scene->m_skyTexture->value(p.x, p.y);
    else{
    ret = scene->skyColor;
    }
  }

  return glm::clamp(ret, color3(0,0,0), color3(1,1,1));
}

void scaleDifferentials(Ray* ray, float s){
  ray->dox = ray->orig + (ray->dox - ray->orig) * s;
  ray->doy = ray->orig + (ray->doy - ray->orig) * s;
  ray->ddx = ray->dir + (ray->ddx - ray->dir) * s;
  ray->ddy = ray->dir + (ray->ddy - ray->dir) * s;
}


static std::minstd_rand engineSamples(time(NULL));
static std::uniform_real_distribution<float> m_unifDistributionSamples{0.0f, 1.0f};

void renderImage(RenderImage *img, Scene *scene)
{

  KdTree *tree = initKdTree(scene);

  auto startTime = std::chrono::system_clock::now();

  auto sampler = 
    new StratifiedSampler(1, 1, true, 2);

  for (size_t j = 0; j < img->height; j++)
  {
    if (j != 0)
      printf("\033[A\r");
    float progress = (float)j / img->height * 100.f;
    printf("progress\t[");
    int cpt = 0;
    for (cpt = 0; cpt < progress; cpt += 5)
      printf(".");
    for (; cpt < 100; cpt += 5)
      printf(" ");
    printf("]\n");
#pragma omp parallel
    {
    int seed = omp_get_thread_num();
    std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < img->width; i++)
    {
      color3 pixel_color(0,0,0);
      color3 *ptr = getPixelPtr(img, i, j);
      auto pixel = vec2(i, j);

      tileSampler->StartPixel(pixel);
      do {
        CameraSample cameraSample = tileSampler->GetCameraSample(pixel);
        Ray* rx = new Ray;
        scene->cam->get_ray(cameraSample.xy.x, cameraSample.xy.y, cameraSample.uv.x, cameraSample.uv.y, rx, vec2(int(i), int(j)));
        scaleDifferentials(rx, 1.f / sqrt(tileSampler->samplesPerPixel));
        Intersection intersection;
        pixel_color += trace_ray(scene, rx, tree, &intersection);
        delete rx;
      }while(tileSampler->StartNextSample());

      *ptr = pixel_color / (float)tileSampler->samplesPerPixel;
    }
  }
  }
  auto stopTime = std::chrono::system_clock::now();
  std::cout << "Rendering took "<< std::chrono::duration_cast<std::chrono::duration<double>>(
                    stopTime - startTime).count() << "s" << std::endl;
}
