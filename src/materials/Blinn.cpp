#include "Blinn.h"
#include "../bsdf.hpp"

Blinn::Blinn(){
  m_IOR = 1.0;
  m_specularColor = color3(0.7f);
  m_diffuseColor = color3(0.5f);
  m_shininess = 20.f;
  m_reflection = vec3(0,0,0);
  m_refraction = vec3(0,0,0);
  m_absorption = vec3(0,0,0);
}

color3 Blinn::shade(Intersection *intersection, vec3 v, Light* light, float intensity) {
  color3 ret = color3(0.f);
  vec3 n = intersection->normal;
  auto lc = light->getColor();
  auto samples = light->getSamples(); 
  
  if(intersection->isOutside){
    for(auto& sample : samples){
      vec3 l;
      if(light->isDirectional()){
        l = -sample;
      }
      else{
        l = normalize(sample - intersection->position);
      }
      float LdotN = dot(l, n);
      v = normalize(v);
      vec3 vl = v + l;
      vec3 h = vl;
      h = normalize(h);
      float NdotH = dot(n, h);
      float s = std::pow(NdotH, m_shininess);

      if(LdotN >= 0){
        if(m_texture != nullptr){
          color3 texColor = m_texture->value(intersection->u, intersection->v, intersection->duv);
          ret += lc *  LdotN * ( texColor + s * m_specularColor);
        }
        else{
          ret += lc *  (LdotN * m_diffuseColor + s * m_specularColor) ;
        }
      }
    }
    ret = (ret / float(samples.size())) * intensity;

  }

  return ret;
}

Blinn::~Blinn() {}

color3 Blinn::textureColor(float u, float v, int face) {
  if(m_texture != nullptr)
    return m_texture->value(u, v);
  return color3(0.f);
}

color3 Blinn::ambientColor(Intersection* intersection, color3 lightColor) {
  if(m_texture != nullptr){
    color3 texColor = m_texture->value(intersection->u, intersection->v, intersection->duv);
    return texColor * lightColor;
  }
  return m_diffuseColor * lightColor;
}

color3 Blinn::scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) {
  auto ret = color3(0.0);
  auto reflectionShade = color3(0.f);
  auto refractionShade = color3(0.f);
  if((m_reflection.x > 0.f || m_reflection.y > 0.f || m_reflection.z > 0) || (m_refraction.x > 0.f || m_refraction.y > 0.f || m_refraction.z > 0))
  {
    vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
    vec3 r = reflect(ray->dir, normal);
    Ray ray_ref;
    ray_ref.difx = new Ray();
    ray_ref.dify = new Ray();
    rayInit(&ray_ref, intersection->position + (acne_eps * r), r, ray->pixel,0, 100000, ray->depth + 1);
    vec3 difx_dir = reflect(ray->difx->dir, normal);
    vec3 dify_dir = reflect(ray->dify->dir, normal);
    rayInit(ray_ref.difx, intersection->position + (acne_eps * difx_dir), difx_dir, ray->pixel,0, 100000, ray->depth + 1);
    rayInit(ray_ref.dify, intersection->position + (acne_eps * dify_dir), dify_dir, ray->pixel,0, 100000, ray->depth + 1);

    Intersection temp_intersection;
    reflectionShade = trace_ray(scene, &ray_ref, tree, &temp_intersection);

    if(reflectionShade == scene->skyColor){
        vec3 dir = r;
        float z = asin(-dir.z) / float(M_PI) + 0.5;
        float x = dir.x / (abs(dir.x) + abs(dir.y) + 0.00001);
        float y = dir.y / (abs(dir.x) + abs(dir.y) + 0.00001);
        point3 p = point3(0.5, 0.5, 0.0) + z * (x * point3(0.5, 0.5, 0.0) + y * point3(-0.5, 0.5, 0.0));
        color3 env = 0.3f * scene->m_skyTexture->value(p.x, p.y);
        ret += m_reflection * env;
    }
    else if(intersection->isOutside)
      ret += (reflectionShade * m_reflection);
    delete ray_ref.difx;
    delete ray_ref.dify;
  }
  if(m_refraction.x > 0.f || m_refraction.y > 0.f || m_refraction.z > 0) {

    float n1, n2;
    vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
    if(intersection->isOutside){
      n1 = 1.0;
      n2 = m_IOR;
    }else {
      n2 = 1.0;
      n1 = m_IOR;
    }
    vec3 v = -ray->dir;

    // calculate refraction ray direction
    float c1 = dot(normal, v);
    float s1 = sqrt(1.0 - c1 * c1);
    float s2 = n1 / n2 * s1;
    float c2 = sqrt(1.0 - s2 * s2);
    vec3 p = normalize(v - c1 * normal);
    vec3 pt = s2 * -p;
    vec3 nt = c2 * -normal;

    vec3 refractDir = normalize(pt + nt);

    if(s2 * s2 <= 1.0){
      Ray ray_refr;
      ray_refr.difx = new Ray();
      ray_refr.dify = new Ray();
      vec3 difx_dir = refract(ray->difx->dir, normal, n1 / n2);
      vec3 dify_dir = refract(ray->dify->dir, normal, n1 / n2);
      rayInit(&ray_refr, intersection->position + (acne_eps * refractDir), refractDir, ray->pixel,0, 100000, ray->depth + 1);
      rayInit(ray_refr.difx, intersection->position + (acne_eps * difx_dir), difx_dir, ray->pixel,0, 100000, ray->depth + 1);
      rayInit(ray_refr.dify, intersection->position + (acne_eps * dify_dir), dify_dir, ray->pixel,0, 100000, ray->depth + 1);

      Intersection temp_inter;
      refractionShade= trace_ray(scene, &ray_refr, tree, &temp_inter);

      if(refractionShade != scene->skyColor){
        // Schlick's approximation for transmittance vs. reflectance
        float r0 = (n1 - n2) / (n1 + n2);
        r0 *= r0;
        float r;
        if(n1 <= n2)
          r = r0 + (1.0 - r0) * (1 - c1) * (1 - c1) * (1 - c1) * (1 - c1) * (1 - c1);
        else
          r = r0 + (1.0 - r0) * (1 - c2) * (1 - c2) * (1 - c2) * (1 - c2) * (1 - c2);
        float t = 1.0 - r;
        
        // compute total refraction color
        color3 refractionColor = m_refraction * (t * refractionShade + r * reflectionShade);

        if(!temp_inter.isOutside){
          refractionColor.r *= exp(-m_absorption.r * ray_refr.tmax);
          refractionColor.g *= exp(-m_absorption.g * ray_refr.tmax);
          refractionColor.b *= exp(-m_absorption.b * ray_refr.tmax);
        }
        ret += refractionColor;
        delete ray_refr.difx;
        delete ray_refr.dify;
      }
      else{
        vec3 dir = refractDir;
        float z = asin(-dir.z) / float(M_PI) + 0.5;
        float x = dir.x / (abs(dir.x) + abs(dir.y) + 0.00001);
        float y = dir.y / (abs(dir.x) + abs(dir.y) + 0.00001);
        point3 p = point3(0.5, 0.5, 0.0) + z * (x * point3(0.5, 0.5, 0.0) + y * point3(-0.5, 0.5, 0.0));
        color3 env = 0.3f * scene->m_skyTexture->value(p.x, p.y);
        ret += m_refraction * env;
      }
    }
    else
    {
      ret += m_refraction * reflectionShade;
    }
  }
  return ret;
}
