
#include <algorithm>

#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "scene.h"
#include <stdio.h>
#include <cmath>
#include <string.h>
#include <iostream>

#include <glm/gtc/epsilon.hpp>
#include <random>
#include <omp.h>

#include "textures.hpp"
#include "intersection.hpp"
#include "bsdf.hpp"

#include "Light.h"

color3 shade(vec3 n, vec3 v, vec3 intersectionPos, color3 lc, Material *mat, float uTex, float vTex, bool outside, float intensity, std::vector<vec3> &samples)
{
  color3 ret = color3(0.f);

  //! \todo compute bsdf, return the shaded color taking into account the
  //! lightcolor

  if (mat->mtype == DIELECTRIC){
    float cpt = 0.f;
    for(auto& sample : samples){
      cpt++;
      vec3 lp = sample - intersectionPos;
      vec3 l = lp / length(lp);
      float LdotN = abs(dot(l, n));
      if(LdotN == 0.0f) continue;
      float VdotN = abs(dot(v, n));
      float extIOR, intIOR;
      if(outside){
        extIOR = mat->IOR;
        intIOR = 1.f;
      }
      else
      {
        extIOR = 1.f;
        intIOR = mat->IOR;
      }

      // REFLECTION
      vec3 hr = (v + l);
      hr = hr / length(hr);
      float LdotH = abs(dot(l, hr));
      float NdotH = abs(dot(n, hr));
      float VdotH = abs(dot(v, hr));
      auto brdfColor = RDM_brdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, extIOR, intIOR);

      // REFRACTION
      hr = -intIOR * l - extIOR * v;
      hr = hr / length(hr);
      LdotH = abs(dot(l, hr));
      NdotH = abs(dot(n, hr));
      VdotH = abs(dot(v, hr));
      auto btdfColor = RDM_btdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, extIOR, intIOR);

      // BSDF
      ret += ( brdfColor + btdfColor ) * LdotN;
    }
    ret = lc * (ret / cpt) * intensity;
  }
  else
  {
    float cpt = 0.f;
    for(auto& sample : samples){
      cpt++;
      vec3 lp = sample - intersectionPos;
      vec3 l = lp / length(lp);
      float LdotN = dot(l, n);
      if (LdotN > 0)
      {
        vec3 vl = v + l;
        vec3 h = vl / length(vl);
        float LdotH = dot(l, h);
        float NdotH = dot(n, h);
        float VdotH = dot(v, h);
        float VdotN = dot(v, n);
        ret += (RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat, uTex, vTex) * LdotN);
      }
    }
    ret = lc * (ret / cpt) * intensity;
  }

  return ret;
}

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree)
{
  color3 ret = color3(0, 0, 0);
  Intersection intersection;

  if (ray->depth > 3)
    return color3(0.f);

  if (intersectKdTree(scene, tree, ray, &intersection))
  {
    size_t lightsCount = scene->lights.size();
    for (size_t i = 0; i < lightsCount; i++)
    {
      vec3 v = ray->dir * -1.0f;
      auto intensity = scene->lights[i]->intensityAt(intersection.position, scene, tree, v, &intersection); 
      auto samples = scene->lights[i]->getSamples();
      if(intensity > 0.0f){
        ret += shade(intersection.normal, v, intersection.position, scene->lights[i]->getColor(), intersection.mat, intersection.u, intersection.v, intersection.isOutside, intensity, samples);
      }
    }

    if (ret.r > 1.f && ret.g > 1.f && ret.b > 1.f && ray->depth > 0) // Si contribution maximale -> on arrete
      return color3(1.f);

    color3 reflectionColor = color3(0.0f);
    color3 refractionColor = color3(0.0f);

    if(intersection.mat->mtype == DIELECTRIC) {
      // REFRACTION + REFLECTION

      vec3 normal = intersection.isOutside ? intersection.normal : -intersection.normal;
      vec3 r = reflect(ray->dir, normal);
      Ray ray_ref;

      rayInit(&ray_ref, intersection.position + (acne_eps * r), r, 0, 100000, ray->depth + 1);
      color3 reflectionColor = trace_ray(scene, &ray_ref, tree);

      float LdotH = dot(ray_ref.dir, normal);

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
      // REFLECTION
      vec3 r = reflect(ray->dir, intersection.normal);
      Ray ray_ref;
      rayInit(&ray_ref, intersection.position + (acne_eps * r), r, 0, 100000, ray->depth + 1);

      color3 cr = trace_ray(scene, &ray_ref, tree);
      float LdotH = dot(ray_ref.dir, intersection.normal);
      float f = RDM_Fresnel(LdotH, 1.f, intersection.mat->IOR);

      reflectionColor = (f * cr * intersection.mat->specularColor);
      
      ret += reflectionColor;
    }
  }
  else
  {
    ret = scene->skyColor;
  }

  return ret;
}


color3 trace_ray_multisampling(Scene *scene, KdTree *tree, int indexI, int indexJ, vec3 dx,
                               vec3 dy, vec3 ray_delta_x, vec3 ray_delta_y)
{

  color3 pixelColor = color3(0.f);

  // We use the same process as for one ray:
  /* vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                     float(i) * dx + float(j) * dy*/
  // We simply need to add a coefficient to i and j to subdivide the pixel in 9 different points/rays
  // from which we will use the average.
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                     (indexI + float(i) / 3.f) * dx + (indexJ + float(j) / 3.f) * dy;
      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      pixelColor += trace_ray(scene, &rx, tree);
    }
  }
  return (pixelColor / 9.f);
}

color3 trace_ray_4multisampling(Scene *scene, KdTree *tree, int indexI, int indexJ, vec3 dx,
                               vec3 dy, vec3 ray_delta_x, vec3 ray_delta_y)
{

  color3 pixelColor = color3(0.f);

  // We use the same process as for one ray:
  /* vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                     float(i) * dx + float(j) * dy*/
  // We simply need to add a coefficient to i and j to subdivide the pixel in 9 different points/rays
  // from which we will use the average.
  for (int i = 1; i <= 3; i+=2)
  {
    for (int j = 1; j <= 3; j+=2)
    {
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                     (indexI + float(i) / 4.f) * dx + (indexJ + float(j) / 4.f) * dy;
      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      pixelColor += trace_ray(scene, &rx, tree);
    }
  }
  return (pixelColor / 4.f);
}

void renderImage(RenderImage *img, Scene *scene)
{

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  float aspect = 1.f / scene->cam.aspect;

  KdTree *tree = initKdTree(scene);

  printf("End building tree\n");

  //! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f);   //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                     aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x =
      (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;

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
      color3 *ptr = getPixelPtr(img, i, j);
      *ptr = trace_ray_4multisampling(scene, tree, i, j, dx, dy, ray_delta_x, ray_delta_y);
    }
  }
}
