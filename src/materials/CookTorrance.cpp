#include "CookTorrance.h"
#include "../bsdf.hpp"

CookTorrance::CookTorrance(bool transparent){
  m_IOR = 1.0;
  m_roughness = 0.1;
  m_specularColor = color3(1.f);
  m_diffuseColor = color3(1.f);
  m_transparent = transparent;
}

color3 CookTorrance::shade(Intersection *intersection, vec3 v, color3 lc, float intensity, std::vector<vec3> &samples) {
  color3 ret = color3(0.f);
  vec3 n = intersection->normal;
  vec3 intersectionPos = intersection->position;
  float uTex = intersection->u;
  float vTex = intersection->v;
  bool outside = intersection->isOutside;
  int face = intersection->face;

  //! \todo compute bsdf, return the shaded color taking into account the
  //! lightcolor

  if (m_transparent){
    for(auto& sample : samples){
      vec3 lp = sample - intersectionPos;
      vec3 l = lp / length(lp);
      float LdotN = abs(dot(l, n));
      if(LdotN == 0.0f) continue;
      float VdotN = abs(dot(v, n));
      float extIOR, intIOR;
      if(outside){
        extIOR = 1.f;
        intIOR = m_IOR;
      }
      else
      {
        extIOR = m_IOR;
        intIOR = 1.f;
      }

      // REFLECTION
      vec3 hr = (v + l);
      hr = hr / length(hr);
      float LdotH = abs(dot(l, hr));
      float NdotH = abs(dot(n, hr));
      float VdotH = abs(dot(v, hr));
      auto brdfColor = RDM_brdf(LdotH, NdotH, VdotH, LdotN, VdotN, extIOR, intIOR, m_roughness, m_specularColor);

      // REFRACTION
      hr = -extIOR * l - intIOR * v;
      hr = hr / length(hr);
      LdotH = abs(dot(l, hr));
      NdotH = abs(dot(n, hr));
      VdotH = abs(dot(v, hr));
      auto btdfColor = RDM_btdf(LdotH, NdotH, VdotH, LdotN, VdotN, extIOR, intIOR, m_roughness, m_specularColor);

      // BSDF
      ret += ( brdfColor + btdfColor ) * LdotN;
    }
    ret = lc * (ret / float(samples.size())) * intensity;
  } else {
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
      ret += lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, m_texture, 
             m_diffuseColor, m_specularColor, m_roughness, m_IOR, uTex, vTex, face) * LdotN;
    }
    ret = (ret / float(samples.size())) * intensity;
  }

  return ret;
}

CookTorrance::~CookTorrance() {}

color3 CookTorrance::textureColor(float u, float v, int face) {
  return m_texture->value(u, v, face);
}

color3 CookTorrance::ambientColor(Intersection* intersection, color3 lightColor) {
  if(m_texture != nullptr)
    return m_texture->value(intersection->u, intersection->v) * lightColor;
  return m_diffuseColor * lightColor;
}

color3 CookTorrance::scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) {
  auto ret = color3(0.0);
  if(m_transparent) {
    // REFRACTION + REFLECTION

    vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
    vec3 r = reflect(ray->dir, normal);
    Ray ray_ref;

    rayInit(&ray_ref, intersection->position + (acne_eps * r), r, 0, 100000, ray->depth + 1);
    Intersection inter;
    color3 reflectionColor = trace_ray(scene, &ray_ref, tree, &inter);

    vec3 v = ray->dir * -1.0f;
    vec3 h = v + ray_ref.dir;
    h = h / length(h);
    float LdotH = dot(ray_ref.dir, h);

    float refractionRatio = 
      intersection->isOutside ? (1.f/m_IOR) : m_IOR;
    float f = intersection->isOutside ? 
        RDM_Fresnel(LdotH, 1.f, m_IOR) : 
        RDM_Fresnel(LdotH, m_IOR, 1.f);

    vec3 unit_direction = normalize( ray->dir );
    auto refractionColor = color3(0.f);
    if( f < 1.0f ){
      vec3 refr = refract(unit_direction, normal, refractionRatio);
      Ray ray_refr;
      rayInit(&ray_refr, intersection->position + (acne_eps * refr), refr, 0, 100000, ray->depth + 1);

      Intersection refInter;
      refractionColor = trace_ray(scene, &ray_refr, tree, &refInter);
    }
    
    ret += ( reflectionColor * f * m_specularColor ) + refractionColor * (1.0f - f);
  } 
  else
  {
    //// REFLECTION
    vec3 r = reflect(ray->dir, intersection->normal);
    Ray ray_ref;
    rayInit(&ray_ref, intersection->position + (acne_eps * r), r, 0, 100000, ray->depth + 1);

    Intersection inter;
    color3 cr = trace_ray(scene, &ray_ref, tree, &inter);
    vec3 v = ray->dir * -1.0f;
    vec3 h = v + ray_ref.dir;
    h = h / length(h);
    float LdotH = dot(ray_ref.dir, h);
    float f = RDM_Fresnel(LdotH, 1.f, m_IOR);

    ret += (f * cr * m_specularColor);

    //auto indirectColor = color3(0.f);
    //float fuzz = 0.0f;
    //for(int i = 0; i < 5; i++){
    //  Ray scattered;
    //  rayInit(&scattered, intersection->position+(acne_eps*r), normalize(r + fuzz * sphereRand()), 0, 10000, ray->depth+1);

    //  if(dot(scattered.invdir, intersection->normal) <= 0)
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
  return ret;
}