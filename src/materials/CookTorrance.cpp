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
    float LdotH = dot(ray_ref.dir, normal);

    float refractionRatio = 
      intersection->isOutside ? (1.f/m_IOR) : m_IOR;
    float f = intersection->isOutside ? 
        RDM_Fresnel(LdotH, 1.f, m_IOR) : 
        RDM_Fresnel(LdotH, m_IOR, 1.f);

    vec3 unit_direction = normalize( ray->dir );

    auto refractionColor = color3(0.f);
    if( f < 1.0f ){
      normal = n2;
      vec3 refr = refract(unit_direction, normal, refractionRatio);
      Ray ray_refr;
      rayInit(&ray_refr, intersection->position + (acne_eps * refr), refr, ray->pixel,0, 100000, ray->depth + 1);

      Intersection refInter;
      refractionColor += trace_ray(scene, &ray_refr, tree, &refInter);
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
    float LdotH = dot(ray_ref.dir, intersection->normal);
    float f = RDM_Fresnel(LdotH, 1.f, m_IOR);

    ret += (f * cr * m_specularColor);

    vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
    color3 iColor = color3(0,0,0);
    vec3 dirA = normalize(random_dir(normal));
    Ray ray_gi;
    rayInit(&ray_gi, intersection->position + (acne_eps * normal), normalize(dirA), ray->pixel,0, 100000, ray->depth + 1);
    
    ray_gi.dox =  vec3(0);
    ray_gi.doy =  vec3(0);
    ray_gi.ddx =  vec3(0);
    ray_gi.ddy =  vec3(0);

    Intersection temp_intersection;
    auto reflColor = trace_ray(scene, &ray_gi, tree, &temp_intersection);

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
    ret += iColor;
  }


  return ret;
}

color3 CookTorrance::sample_f(vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type) {
  if(type == 0){
    bool isOutside = dot(wo, normal) < 0;
    normal = isOutside ? normal : -normal;
    *wi = reflect(wo, normal);
    *pdf = 1;

    vec3 v = wo * -1.0f;
    vec3 h = v + *wi;
    h = h / length(h);
    float LdotH = dot(*wi, h);

    float f = isOutside ? RDM_Fresnel(LdotH, 1.f, m_IOR) : RDM_Fresnel(LdotH, m_IOR, 1.f);
    
    return f * m_specularColor;
  }
  else if(type == 1){
    if(m_transparent) {
      // REFRACTION + REFLECTION

      //vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
      bool isOutside = dot(wo, normal) < 0;
      normal = isOutside ? normal : -normal;

      vec3 r = reflect(-wo, normal);

      vec3 h = wo + r;
      h = h / length(h);
      float LdotH = dot(r, h);

      float refractionRatio = isOutside ? (1.f/m_IOR) : m_IOR;
      float f = isOutside ? RDM_Fresnel(LdotH, 1.f, m_IOR) : RDM_Fresnel(LdotH, m_IOR, 1.f);

      vec3 unit_direction = normalize( -wo );

      if( f < 1.0f ){
        *wi = refract(unit_direction, normal, refractionRatio);
        *pdf = 1;
        return vec3(1.0f - f);
      }
      
    } 

  }
  return color3(0);
}
