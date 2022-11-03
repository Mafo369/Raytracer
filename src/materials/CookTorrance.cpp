#include "CookTorrance.h"
#include "../bsdf.hpp"
#include "../sampling/sampling.h"

CookTorrance::CookTorrance(bool transparent){
  m_IOR = 1.0;
  m_roughness = 0.1;
  m_specularColor = color3(1.f);
  m_diffuseColor = color3(1.f);
  m_transparent = transparent;
}

color3 CookTorrance::shade(Intersection *intersection, vec3 v, Light* light, float intensity) {
  color3 ret = color3(0.f);
  vec3 n = intersection->normal;
  vec3 intersectionPos = intersection->position;
  float uTex = intersection->u;
  float vTex = intersection->v;
  bool outside = intersection->isOutside;
  int face = intersection->face;
  auto lc = light->getColor();
  auto samples = light->getSamples();

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

static std::uniform_real_distribution<float> m_unifDistributionRand{-1.f, 1.0f};

vec3 random_in_unit_sphere() {
    while (true) {
        auto p = vec3(m_unifDistributionRand(engine), m_unifDistributionRand(engine), m_unifDistributionRand(engine));
        if (glm::length_sq(p) >= 1) continue;
        return p;
    }
}

color3 CookTorrance::ambientColor(Ray* ray, Intersection* intersection, color3 lightColor) {
  if(m_texture != nullptr)
    return m_texture->value(intersection->u, intersection->v) * lightColor;
  return m_diffuseColor * lightColor;
}

color3 CookTorrance::scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) {
  auto ret = color3(0.0);
  if(m_transparent) {
    // REFRACTION + REFLECTION

    vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
    point3 v0 = point3(0, 1, 0);
    if(dot(v0, normal))
      v0 = point3(0,0,1);
    point3 v1 = normalize(cross(v0, normal));
    float rnd = sqrt(uniform01(engine));
    float rndReflection = rnd * m_roughness;
    float rndRefraction = rnd * m_roughness;
    float factor = uniform01(engine) * 2.0 * M_PI;
    
    vec3 n1 = normalize(normal + (v0 * rndReflection * cos(factor)) + (v1 * rndReflection * sin(factor)));
    vec3 n2 = normalize(normal + (v0 * rndRefraction * cos(factor)) + (v1 * rndRefraction * sin(factor)));

    normal = n1;

    vec3 r = reflect(ray->dir, normal);
    Ray ray_ref;

    rayInit(&ray_ref, intersection->position + (acne_eps * r), r, ray->pixel,0, 100000, ray->depth + 1);
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

    float fuzz = 0.1f;
    auto refractionColor = color3(0.f);
    if( f < 1.0f ){
      normal = n2;
      vec3 refr = refract(unit_direction, normal, refractionRatio);
      for(int i = 0; i < 20; i++){
        Ray ray_refr;
        rayInit(&ray_refr, intersection->position + (acne_eps * refr), refr + fuzz*random_in_unit_sphere(), ray->pixel,0, 100000, ray->depth + 1);

        Intersection refInter;
        refractionColor += trace_ray(scene, &ray_refr, tree, &refInter);
      }
      refractionColor = refractionColor / 20.f;
    }
    
    ret += ( reflectionColor * f * m_specularColor ) + refractionColor * (1.0f - f);
  } 
  else
  {
    point3 v0 = point3(0, 1, 0);
    if(dot(v0, intersection->normal))
      v0 = point3(0,0,1);
    point3 v1 = normalize(cross(v0, intersection->normal));
    float rnd = sqrt(uniform01(engine));
    float rndReflection = rnd * m_roughness;
    float factor = uniform01(engine) * 2.0 * M_PI;
    
    vec3 n1 = normalize(intersection->normal + (v0 * rndReflection * cos(factor)) + (v1 * rndReflection * sin(factor)));
    //// REFLECTION
    vec3 r = reflect(ray->dir, n1);
    Ray ray_ref;
    rayInit(&ray_ref, intersection->position + (acne_eps * r), r, ray->pixel,0, 100000, ray->depth + 1);

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
  vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
  int isamples = 1;
  color3 iColor = color3(0,0,0);
  for(int i = 0; i < isamples; i++){
    vec3 dirA = cosineWeightedSampling(normal);
    Ray ray_ref;
    rayInit(&ray_ref, intersection->position + (acne_eps * normal), normalize(dirA), ray->pixel,0, 100000, ray->depth + 1);
    
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
        iColor += m_diffuseColor * env;
      }
    }
    else if(intersection->isOutside){
      if(m_texture != nullptr){
        iColor += reflColor * m_texture->value(intersection->u, intersection->v);
      }
      else
      {
        iColor += reflColor * m_diffuseColor;
      }
    }
  }
  iColor /= (float)isamples;
  ret += iColor;
 
  return ret;
}
