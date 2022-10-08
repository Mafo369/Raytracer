
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

color3 shade(vec3 n, vec3 v, vec3 intersectionPos, color3 lc, const Material *mat, float uTex, float vTex, bool outside, float intensity, std::vector<vec3> &samples, int face)
{
  color3 ret = color3(0.f);

  //! \todo compute bsdf, return the shaded color taking into account the
  //! lightcolor

  if (mat->mtype == TRANSPARENT){
    for(auto& sample : samples){
      vec3 lp = sample - intersectionPos;
      vec3 l = lp / length(lp);
      float LdotN = abs(dot(l, n));
      if(LdotN == 0.0f) continue;
      float VdotN = abs(dot(v, n));
      float extIOR, intIOR;
      if(outside){
        extIOR = 1.f;
        intIOR = mat->IOR;
      }
      else
      {
        extIOR = mat->IOR;
        intIOR = 1.f;
      }

      // REFLECTION
      vec3 hr = (v + l);
      hr = hr / length(hr);
      float LdotH = abs(dot(l, hr));
      float NdotH = abs(dot(n, hr));
      float VdotH = abs(dot(v, hr));
      auto brdfColor = RDM_brdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, extIOR, intIOR);

      // REFRACTION
      hr = -extIOR * l - intIOR * v;
      hr = hr / length(hr);
      LdotH = abs(dot(l, hr));
      NdotH = abs(dot(n, hr));
      VdotH = abs(dot(v, hr));
      auto btdfColor = RDM_btdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, extIOR, intIOR);

      // BSDF
      ret += ( brdfColor + btdfColor ) * LdotN;
    }
    ret = lc * (ret / float(samples.size())) * intensity;
  } else if (mat->mtype == TORSPAR){
    for(auto& sample : samples){
      vec3 lp = sample - intersectionPos;
      vec3 l = lp / length(lp);
      float LdotN = abs(dot(l, n));
      float VdotN = abs(dot(v, n));
      if(LdotN == 0.0f || VdotN == 0.0f) continue;
      float extIOR, intIOR;
      if(outside){
        extIOR = 1.f;
        intIOR = mat->IOR;
      }
      else
      {
        extIOR = mat->IOR;
        intIOR = 1.f;
      }

      // REFLECTION
      vec3 hr = (v + l);
      if(hr.x == 0.0f && hr.y == 0.0f && hr.z == 0) continue;
      hr = hr / length(hr);
      float LdotH = abs(dot(l, hr));
      float NdotH = abs(dot(n, hr));
      float VdotH = abs(dot(v, hr));
      auto brdfColor = RDM_brdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, extIOR, intIOR);

      // REFRACTION
      //hr = -extIOR * l - intIOR * v;
      //hr = hr / length(hr);
      //LdotH = abs(dot(l, hr));
      //NdotH = abs(dot(n, hr));
      //VdotH = abs(dot(v, hr));
      //auto btdfColor = RDM_btdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, extIOR, intIOR);

      //// BSDF
      //ret += ( brdfColor + btdfColor ) * LdotN;

      ret += brdfColor * LdotN;
    }
    ret = lc * (ret / float(samples.size())) * intensity;
  }
  else
  {
    for(auto& sample : samples){
      vec3 lp = sample - intersectionPos;
      vec3 l = lp / length(lp);
      float LdotN = abs(dot(l, n));
      vec3 vl = v + l;
      vec3 h = vl;
      if(h.x == 0 && h.y == 0 && h.z == 0) 
        continue;
      h = normalize(h);
      float LdotH = dot(l, h);
      float NdotH = dot(n, h);
      float VdotH = dot(v, h);
      float VdotN = abs(dot(v, n));
      ret += lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, uTex, vTex, face) * LdotN;
    }
    ret = (ret / float(samples.size())) * intensity;
  }

  return ret;
}

vec3 sphereRand(){
  while(true){
    vec3 p(m_unifDistributionRand(engine), m_unifDistributionRand(engine), 
          m_unifDistributionRand(engine));
    if (glm::length_sq(p) >= 1) continue;
    return p;
  }
}

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree)
{
  color3 ret = color3(0, 0, 0);
  Intersection intersection;

  if (ray->depth > 10)
    return color3(0.f);

  if (intersectKdTree(scene, tree, ray, &intersection))
  {
    if(intersection.face != -1)
      return intersection.mat->m_texture->value(intersection.u, intersection.v, intersection.face); 
    size_t lightsCount = scene->lights.size();
    for (size_t i = 0; i < lightsCount; i++)
    {
      vec3 v = ray->dir * -1.0f;
      if(scene->lights[i]->isAmbient()){
        ret += intersection.mat->diffuseColor * scene->lights[i]->getColor();
        continue;
      }
      auto intensity = scene->lights[i]->intensityAt(intersection.position, scene, tree, v, &intersection); 
      auto samples = scene->lights[i]->getSamples();
      if(intensity > 0.0f){
        ret += shade(intersection.normal, v, intersection.position, scene->lights[i]->getColor(), intersection.mat, intersection.u, intersection.v, intersection.isOutside, intensity, samples, intersection.face);
      }
    }

    if (ret.r >= 1.f && ret.g >= 1.f && ret.b >= 1.f && ray->depth > 0) // Si contribution maximale -> on arrete
      return color3(1.f);

    color3 refractionColor = color3(0.0f);

    if(intersection.mat->mtype == TRANSPARENT) {
      // REFRACTION + REFLECTION

      vec3 normal = intersection.isOutside ? intersection.normal : -intersection.normal;
      vec3 r = reflect(ray->dir, normal);
      Ray ray_ref;

      rayInit(&ray_ref, intersection.position + (acne_eps * r), r, 0, 100000, ray->depth + 1);
      color3 reflectionColor = trace_ray(scene, &ray_ref, tree);

      vec3 v = ray->dir * -1.0f;
      vec3 h = v + ray_ref.dir;
      h = h / length(h);
      float LdotH = dot(ray_ref.dir, h);

      float refractionRatio = 
        intersection.isOutside ? (1.f/intersection.mat->IOR) : intersection.mat->IOR;
      float f = intersection.isOutside ? 
          RDM_Fresnel(LdotH, 1.f, intersection.mat->IOR) : 
          RDM_Fresnel(LdotH, intersection.mat->IOR, 1.f);

      vec3 unit_direction = normalize( ray->dir );
      if( f < 1.0f ){
        vec3 refr = refract(unit_direction, normal, refractionRatio);
        Ray ray_refr;
        rayInit(&ray_refr, intersection.position + (acne_eps * refr), refr, 0, 100000, ray->depth + 1);

        refractionColor = trace_ray(scene, &ray_refr, tree);
      }
      
      ret += ( reflectionColor * f * intersection.mat->specularColor ) + refractionColor * (1.0f - f);
    } 
    else
    {
      //// REFLECTION
      vec3 r = reflect(ray->dir, intersection.normal);
      Ray ray_ref;
      rayInit(&ray_ref, intersection.position + (acne_eps * r), r, 0, 100000, ray->depth + 1);

      color3 cr = trace_ray(scene, &ray_ref, tree);
      vec3 v = ray->dir * -1.0f;
      vec3 h = v + ray_ref.dir;
      h = h / length(h);
      float LdotH = dot(ray_ref.dir, h);
      float f = RDM_Fresnel(LdotH, 1.f, intersection.mat->IOR);

      ret += (f * cr * intersection.mat->specularColor);

      //auto indirectColor = color3(0.f);
      //float fuzz = 0.0f;
      //for(int i = 0; i < 5; i++){
      //  Ray scattered;
      //  rayInit(&scattered, intersection.position+(acne_eps*r), normalize(r + fuzz * sphereRand()), 0, 10000, ray->depth+1);

      //  if(dot(scattered.invdir, intersection.normal) <= 0)
      //    continue;

      //  vec3 cr = trace_ray(scene, &scattered, tree);
      //  vec3 h = v + scattered.dir;
      //  h = h / length(h);
      //  float LdotH = dot(scattered.dir, h);
      //  float f = RDM_Fresnel(LdotH, 1.f, intersection.mat->IOR);

      //  indirectColor += (f * cr * intersection.mat->specularColor);
      //}
      //indirectColor /= 5.f;

      //ret += indirectColor;

    }
  }
  else
  {
    ret = scene->skyColor;
  }

  return glm::clamp(ret, color3(0,0,0), color3(1,1,1));
}

static std::minstd_rand engineSamples(time(NULL));
static std::uniform_real_distribution<float> m_unifDistributionSamples{0.0f, 1.0f};

void renderImage(RenderImage *img, Scene *scene)
{

  auto samplesPerPixel = 10;

  KdTree *tree = initKdTree(scene);

  auto startTime = std::chrono::system_clock::now();

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
#pragma omp parallel for
    for (size_t i = 0; i < img->width; i++)
    {
      color3 pixel_color(0,0,0);
      color3 *ptr = getPixelPtr(img, i, j);
      for (int s = 0; s < samplesPerPixel; ++s) {
        auto u = (i + m_unifDistributionSamples(engineSamples)) / (img->width-1);
        auto v = (j + m_unifDistributionSamples(engineSamples)) / (img->height-1);
        Ray r;
        scene->cam->get_ray(u, v, &r);
        pixel_color += trace_ray(scene, &r, tree);
      }
      // Divide the color by the number of samples and gamma-correct for gamma=2.0.
      pixel_color /= samplesPerPixel;
      //pixel_color = glm::sqrt(pixel_color);
      pixel_color = glm::clamp(pixel_color, 0.0f, 1.0f);

      *ptr = pixel_color;
    }
  }
  auto stopTime = std::chrono::system_clock::now();
  std::cout << "Rendering took "<< std::chrono::duration_cast<std::chrono::duration<double>>(
                    stopTime - startTime).count() << "s" << std::endl;
}
