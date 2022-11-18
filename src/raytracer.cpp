
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

color3 EstimateDirect(Ray* ray, const Intersection& it, const point2& uScattering, Light* light,
                      const point2& uLight, Scene* scene, KdTree* tree, Sampler* sampler, bool specular) {
  color3 Ld(0.f);
  vec3 wi;
  float lightPdf = 0, scatteringPdf = 0;
  color3 Li = light->sample_Li(it, uLight, &wi, &lightPdf);
  if(lightPdf > 0 && !isBlack(Li)){
    color3 f;
    if(it.normal != vec3(0)){
      const Intersection& isect = (const Intersection&)it;
      vec3 wo = -ray->dir;
      vec3 normal = isect.isOutside ? isect.normal : -isect.normal;
      f = it.mat->f(wo, wi, normal) * abs(dot(wi, normal));
      scatteringPdf = it.mat->pdf(wo, wi, normal);
    }
    else
    {
      f = color3(0);
      scatteringPdf = 1.f;
    }
    if(!isBlack(f)){
      if(!isBlack(Li)){
        float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
        Ld += f * Li * weight / lightPdf;
      }
    }
  }

  color3 f;
  bool sampledSpecular = false;
  if(it.normal != vec3(0)){
    const Intersection& isect = (const Intersection&)it;
    vec3 wo = -ray->dir;
    vec3 normal = it.isOutside ? it.normal : -it.normal;
    f = isect.mat->sample_f(wo, &wi, normal, uScattering, &scatteringPdf, 2);
    f *= abs(dot(wi, it.normal));
  }
  else
  {
    f = color3(0);
    scatteringPdf = 0.f;
  }
  if(!isBlack(f) && scatteringPdf > 0.f){
    float weight = 1;
    if(!sampledSpecular){
      lightPdf = light->pdf_Li(it, wi);
      if(lightPdf == 0)
        return Ld;
      weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
    }
    Intersection lightIsect;
    Ray new_ray;
    new_ray.hasDifferentials = false;
    rayInit(&new_ray, it.position + (acne_eps * wi), wi, ray->pixel, 0, 10000);
    bool foundSurfaceInteraction = intersectKdTree(scene, tree, &new_ray, &lightIsect);
    
    color3 Li(0.f);
    if(!foundSurfaceInteraction) {
      Li = color3(light->getColor());
    }

    if(!isBlack(Li))
      Ld += f * Li * weight / scatteringPdf;
  }
  return Ld;
}

float int_exponential(float y0, float ysol, float beta, float s, float uy){
  float result = 0.1 * exp((-y0 + ysol) * beta) * (1.f - exp(-s * uy*beta)) / (uy*beta);
  return result;
}

// pathtracer
color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection, bool show_lights, Sampler* sampler)
{
  color3 ret = color3(0, 0, 0);

  if (ray->depth > scene->depth)
    return color3(0.f);

  if (intersectKdTree(scene, tree, ray, intersection))
  {
    intersection->hit = true;
    // Compute necessary differential information for texture filtering
    if(ray->hasDifferentials)
      intersection->computeDifferentials(ray);
    
    // Scatter
    if(!isBlack(intersection->mat->m_emission)){
      auto lightRadius = scene->objects[scene->objects.size()-1]->geom.sphere.radius;
      return show_lights ? (intersection->mat->m_emission / (lightRadius * lightRadius )) : vec3(0);
    }
    else
    {
      ret += intersection->mat->scatterColor(scene, tree, ray, intersection);
    }
  }
  else
  {
    if(scene->m_skyTexture != nullptr && ray->depth == 0){
      vec3 pixelUV = vec3((float)ray->pixel.x / scene->cam->imgWidth, (float)ray->pixel.y / scene->cam->imgHeight, 0.f);
      vec3 p = glm::mod(scene->m_skyTexture->m_transform.transformTo(pixelUV), 1.f);
      ret = scene->skyColor * scene->m_skyTexture->value(p.x, p.y);
    }
    else{
      ret = scene->skyColor;
    }
  }

  return ret;

  //if(ray->tmax < 0 || !intersection->hit || ray->dir.z == 0) return ret;

  //vec3 Lv(0.f);
  //
  //float ysol = scene->ysol;

  //bool is_uniform_fog = false;
  //float beta = is_uniform_fog ? 0.04 : 0.04;
  //float p_uniform = 0.25;

  //bool uniform_sampling_ray = true;
  //int phase = 0; // 0: uniform, 1: schlick, 2: rayleigh

  //float int_ext;
  //if(is_uniform_fog){
  //  int_ext = beta * ray->tmax;
  //}
  //else
  //{
  //  int_ext = int_exponential(ray->orig.z, ysol, beta, ray->tmax, ray->dir.z);
  //}
  //float T = exp(-int_ext);

  //float randt, probat;
  //float clamped_t = min(10000.f, ray->tmax);
  //if(uniform_sampling_ray){
  //  randt = uniform01(engine) * clamped_t;
  //  probat = 1. / clamped_t;
  //}
  //else
  //{
  //  float alpha = 5. / clamped_t;
  //  do{
  //    randt = -log(uniform01(engine)) / alpha;
  //  } while(randt > clamped_t);
  //  float normalization = 1.f / alpha * (1.f - exp(-alpha * clamped_t));
  //  probat = exp(-alpha * randt) / normalization;
  //}

  //float int_ext_partiel;
  //if(is_uniform_fog){
  //  int_ext_partiel = beta * randt;
  //}
  //else
  //{
  //  int_ext_partiel = int_exponential(ray->orig.z, ysol, beta, randt, ray->dir.z);
  //}
  //vec3 randP = ray->orig + randt * ray->dir;

  //vec3 randDir;
  //float probaDir;
  //auto sphereL = scene->objects[scene->objects.size()-1];
  //vec3 axePO = normalize(randP - sphereL->geom.sphere.center);
  //vec3 ptA;

  //bool is_uniform;
  //if(uniform01(engine) < p_uniform){
  //  randDir = random_uniform();
  //  is_uniform = true;
  //}
  //else
  //{
  //  vec3 dirA = random_dir(axePO);
  //  ptA = dirA * sphereL->geom.sphere.radius + sphereL->geom.sphere.center;
  //  randDir = normalize(ptA - randP);
  //  is_uniform = false;
  //}

  //float phase_f;
  //float k = 0.4;
  //switch(phase){
  //  case 0: 
  //    phase_f = 0.3f / (4.f * M_PI);
  //    break;
  //  case 1:
  //    phase_f = (1. - (k*k)) / (4.f * Pi * (1. + k * dot(randDir, -ray->dir)));
  //    break;
  //  case 2:
  //    phase_f = 3./(16. * Pi) * (1.f + sqr(dot(randDir, ray->dir)));
  //    break;
  //  default:
  //    phase_f = 0.3f / (4.f * M_PI);
  //    break;
  //}


  //Ray L_Ray;
  //L_Ray.hasDifferentials = false;
  //rayInit(&L_Ray, randP, randDir, ray->pixel, 0, 100000, ray->depth+1);
  //Intersection interL;
  //color3 L = trace_ray(scene, &L_Ray, tree, &interL);

  //float V;
  //if(is_uniform){
  //  V = 1;
  //}
  //else
  //{
  //  float d_light2 = length_sq(ptA - randP);
  //  if(interL.hit && L_Ray.tmax * L_Ray.tmax < d_light2*0.9){
  //    V = 0;
  //  }else
  //  {
  //    V = 1;
  //  }
  //}

  //if(V == 0) {
  //  Lv = vec3(0);
  //}
  //else
  //{
  //  vec3 interN = interL.normal;
  //  vec3 interP = interL.position;

  //  float pdf_uniform = 1.f / (4.f * Pi);
  //  float J = dot(interN, -randDir) / glm::length_sq(interP - randt);
  //  float pdf_light = (interL.hit && !isBlack(interL.mat->m_emission)) ? (dot(normalize(interP - sphereL->geom.sphere.center), axePO) / (Pi * sqr(sphereL->geom.sphere.radius)) / J) : 0.f;

  //  probaDir = p_uniform * pdf_uniform + (1.f - p_uniform) * pdf_light;

  //  float ext;
  //  if(is_uniform_fog){
  //    ext = beta;
  //  }
  //  else
  //  {
  //    ext = 0.1 * exp(-beta * (randP.z - ysol));
  //  }

  //  Lv = L * phase_f * ext * exp(-int_ext_partiel) / (probat * probaDir); 
  //}

  //return ret * T + Lv;
}

// whitted
//color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection, bool show_lights, Sampler* sampler)
//{
//  color3 ret = color3(0, 0, 0);
//
//  if (ray->depth > 3)
//    return color3(0.f);
//
//  if (intersectKdTree(scene, tree, ray, intersection))
//  {
//    intersection->hit = true;
//    // Compute necessary differential information for texture filtering
//    if(ray->hasDifferentials)
//      intersection->computeDifferentials(ray);
//
//    // if skybox return directly the corresponding color
//    if(intersection->face != -1){
//      return intersection->mat->textureColor(intersection->u, intersection->v, intersection->face);
//    }
//
//    // shading
//    size_t lightsCount = scene->lights.size();
//    for (size_t i = 0; i < lightsCount; i++)
//    {
//      vec3 v = ray->dir * -1.0f;
//      if(scene->lights[i]->isAmbient()){
//        ret += intersection->mat->ambientColor(ray, intersection, scene->lights[i]->getColor());
//        continue;
//      }
//      auto intensity = scene->lights[i]->intensityAt(intersection->position, scene, tree, v, intersection); 
//      if(intensity > 0.0f){
//        ret += intersection->mat->shade(intersection, v, scene->lights[i], intensity);
//      }
//    }
//
//    // If max contribution, we stop
//    if (ret.r >= 1.f && ret.g >= 1.f && ret.b >= 1.f && ray->depth > 0)
//      return color3(1.f);
//
//    // Scatter
//    //ret += intersection->mat->scatterColor(scene, tree, ray, intersection);
//    ret += specularReflect(ray, intersection, scene, tree, sampler);
//    ret += specularTransmission(ray, intersection, scene, tree, sampler);
//  }
//  else
//  {
//    if(scene->m_skyTexture != nullptr && ray->depth == 0){
//      vec3 pixelUV = vec3((float)ray->pixel.x / scene->cam->imgWidth, (float)ray->pixel.y / scene->cam->imgHeight, 0.f);
//      vec3 p = glm::mod(scene->m_skyTexture->m_transform.transformTo(pixelUV), 1.f);
//      ret = scene->skyColor * scene->m_skyTexture->value(p.x, p.y);
//    }
//    else{
//      ret = scene->skyColor;
//    }
//  }
//
//  return glm::clamp(ret, color3(0,0,0), color3(1,1,1));
//}

void scaleDifferentials(Ray* ray, float s){
  ray->dox = ray->orig + (ray->dox - ray->orig) * s;
  ray->doy = ray->orig + (ray->doy - ray->orig) * s;
  ray->ddx = ray->dir + (ray->ddx - ray->dir) * s;
  ray->ddy = ray->dir + (ray->ddy - ray->dir) * s;
}

void renderImage(RenderImage *img, Scene *scene)
{

  KdTree *tree = initKdTree(scene);

  auto startTime = std::chrono::system_clock::now();

  auto sampler = 
    new StratifiedSampler(8, 8, true, 1);
  std::cout << "Spp: " << sampler->samplesPerPixel << std::endl;

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
    std::unique_ptr<Sampler> tileSampler = sampler->Clone(time(NULL));
#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < img->width; i++)
    {
      color3 pixel_color(0,0,0);
      color3 *ptr = getPixelPtr(img, i, j);
      auto pixel = vec2(i, j);

      tileSampler->StartPixel(pixel);
      do {
        CameraSample cameraSample = tileSampler->GetCameraSample(pixel);
        Ray rx;
        rx.hasDifferentials = false;
        scene->cam->get_ray(cameraSample.xy.x, cameraSample.xy.y, cameraSample.uv.x, cameraSample.uv.y, &rx, vec2(int(i), int(j)));
        if(rx.hasDifferentials)
          scaleDifferentials(&rx, 1.f / sqrt(tileSampler->samplesPerPixel));
        Intersection intersection;
        pixel_color += trace_ray(scene, &rx, tree, &intersection, sampler);
      }while(tileSampler->StartNextSample());

      color3 avgColor = pixel_color / (float)tileSampler->samplesPerPixel;

      // gamma-correction
      avgColor.r = powf(avgColor.r, 1.0f  / 2.2);
      avgColor.g = powf(avgColor.g, 1.0f  / 2.2);
      avgColor.b = powf(avgColor.b, 1.0f  / 2.2);

      *ptr = avgColor;
    }
  }
  }
  auto stopTime = std::chrono::system_clock::now();
  std::cout << "Rendering took "<< std::chrono::duration_cast<std::chrono::duration<double>>(
                    stopTime - startTime).count() << "s" << std::endl;
}
