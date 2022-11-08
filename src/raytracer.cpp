
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

float Phong_BRDF(const vec3& wi, const vec3& wo, const vec3& N, float phongExponent){
  vec3 reflected = reflect(wo, N);
  float lobe = pow(dot(reflected, wi), phongExponent) * (phongExponent + 2) / (2.f * Pi);
  return lobe;
}

vec3 random_Phong(const vec3& R, float phong_exponent){
  float r1 = uniform01(engine);
  float r2 = uniform01(engine);

  float theta = 2.f * Pi * r1;
  float r = sqrt(1.f-pow(r2, 2./(phong_exponent + 1)));
  vec3 dir_rand_local = vec3(cos(theta)*r, sin(theta)*r, pow(r2, 1./(phong_exponent+1)));
  vec3 randV = vec3(uniform01(engine)-0.5, uniform01(engine)-0.5, uniform01(engine)-0.5);

  vec3 tangent1 = normalize(cross(R, randV));
  vec3 tangent2 = cross(tangent1, R);

  return dir_rand_local[2] * R + dir_rand_local[0] * tangent1 + dir_rand_local[1] * tangent2;
}

// pathtracer
color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection, Sampler* sampler)
{
  color3 ret = color3(0, 0, 0);

  if (ray->depth > 3)
    return color3(0.f);

  if (intersectKdTree(scene, tree, ray, intersection))
  {
    intersection->hit = true;
    // Compute necessary differential information for texture filtering
    if(ray->hasDifferentials)
      intersection->computeDifferentials(ray);

    // if skybox return directly the corresponding color
    if(intersection->face != -1){
      return intersection->mat->textureColor(intersection->u, intersection->v, intersection->face);
    }

    //if(!isBlack(intersection->mat->m_emission)){
    //  return intersection->mat->m_diffuseColor * intersection->mat->m_emission / (4.f * Pi  * (5.f * 5.f));
    //}

    vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;

    auto sphereL = scene->lights[0];
    vec3 axePO = intersection->position - sphereL->getPosition();
    vec3 dirA = normalize(random_dir(axePO));
    vec3 ptA = dirA * sphereL->getSize() + sphereL->getPosition();
    vec3 wi = normalize(ptA - intersection->position);

    float d_light2 = glm::length_sq(ptA - intersection->position);

    vec3 Np = dirA;

    Intersection temp_inter;
    Ray rayS;
    vec3 origin = intersection->position + (acne_eps * normal);
    rayInit(&rayS, origin, wi, vec2(0,0),0.f, distance(origin, ptA));
    rayS.shadow = true;
    rayS.dox = vec3(0.f);
    rayS.doy = vec3(0.f);
    rayS.ddx = vec3(0.f);
    rayS.ddy = vec3(0.f);
    if(intersectKdTree(scene, tree, &rayS, &temp_inter) && rayS.tmax * rayS.tmax < d_light2*0.99){
      ret += color3(0);
    }
    else
    {
      color3 brdf = intersection->mat->m_diffuseColor / Pi * 
                    (color3(1) - intersection->mat->m_specularColor) + 
                    intersection->mat->m_specularColor * 
                    Phong_BRDF(wi, ray->dir, normal, 100)*intersection->mat->m_diffuseColor;
      float J = 1. * dot(Np, -wi) / d_light2;
      float pdf = dot(axePO, dirA) / (Pi * sphereL->getSize() * sphereL->getSize());
      ret += sphereL->getColor() * max(0.f, dot(normal, wi)) * J * brdf / pdf;
    }

    //// Indirect illumination
    color3 iColor = color3(0,0,0);
    auto specularColor = intersection->mat->m_specularColor;
    float proba = 1. - max(specularColor.r, max(specularColor.g, specularColor.b));
    vec3 R = normalize(reflect(ray->dir, normal));
    bool sample_diffuse;
    if(uniform01(engine) < proba){
      sample_diffuse = true;
      dirA = normalize(random_dir(normal));
    }
    else
    {
      sample_diffuse = false;
      dirA = random_Phong(R, 100);
      if(dot(dirA, normal) < 0) return color3(0);
      if(dot(dirA, R) < 0) return color3(0);
    }


    Ray ray_ref;
    rayInit(&ray_ref, intersection->position + (acne_eps * intersection->normal), normalize(dirA), ray->pixel,0, 100000, ray->depth + 1);
    ray_ref.dox =  vec3(0);
    ray_ref.doy =  vec3(0);
    ray_ref.ddx =  vec3(0);
    ray_ref.ddy =  vec3(0);

    Intersection temp_intersection;
    auto reflColor = trace_ray(scene, &ray_ref, tree, &temp_intersection);

    if(!temp_intersection.hit){
      if(scene->m_skyTexture != nullptr){
        vec3 dir = dirA;
        float z = asin(-dir.z) / float(M_PI) + 0.5;
        float x = dir.x / (abs(dir.x) + abs(dir.y));
        float y = dir.y / (abs(dir.x) + abs(dir.y));
        point3 p = point3(0.5, 0.5, 0.0) + z * (x * point3(0.5, 0.5, 0.0) + y * point3(-0.5, 0.5, 0.0));
        // TODO: Multiply with intensity var
        color3 env = 0.7f * scene->m_skyTexture->value(p.x, p.y);
        iColor += intersection->mat->m_diffuseColor * env;
      }
    }
    else {
      float pdfPhong = (100+1.f) / (2.f * Pi) * pow(dot(R, dirA), 100);
      float pdf = proba * dot(normal, dirA)/(2.f*Pi) + (1.f - proba) * pdfPhong;
      if(sample_diffuse)
        iColor += reflColor * intersection->mat->m_diffuseColor * dot(normal, dirA)/(2.f*Pi) / pdf;
      else
        iColor += reflColor * dot(normal, dirA) * Phong_BRDF(dirA, ray->dir, normal, 100)/(2.f*Pi) / pdf;
    }
    ret += iColor;

    //// If max contribution, we stop
    //if (ret.r >= 1.f && ret.g >= 1.f && ret.b >= 1.f && ray->depth > 0)
    //  return color3(1.f);

    // Scatter
    //ret += intersection->mat->scatterColor(scene, tree, ray, intersection);
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

  return glm::clamp(ret, color3(0,0,0), color3(1,1,1));
}

// whitted
//color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection, Sampler* sampler)
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
    new StratifiedSampler(12, 12, true, 1);

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
        rx.hasDifferentials = true;
        scene->cam->get_ray(cameraSample.xy.x, cameraSample.xy.y, cameraSample.uv.x, cameraSample.uv.y, &rx, vec2(int(i), int(j)));
        if(rx.hasDifferentials)
          scaleDifferentials(&rx, 1.f / sqrt(tileSampler->samplesPerPixel));
        Intersection intersection;
        pixel_color += trace_ray(scene, &rx, tree, &intersection);
      }while(tileSampler->StartNextSample());

      color3 avgColor = pixel_color / (float)tileSampler->samplesPerPixel;

      // gamma-correction
      avgColor.r = pow(avgColor.r, 1.0  / 2.2);
      avgColor.g = pow(avgColor.g, 1.0  / 2.2);
      avgColor.b = pow(avgColor.b, 1.0  / 2.2);

      *ptr = avgColor;
    }
  }
  }
  auto stopTime = std::chrono::system_clock::now();
  std::cout << "Rendering took "<< std::chrono::duration_cast<std::chrono::duration<double>>(
                    stopTime - startTime).count() << "s" << std::endl;
}
