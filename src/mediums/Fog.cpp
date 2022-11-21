#include "Fog.h"
#include "../sampling/sampling.h"
#include "../raytracer.h"

float Fog::int_exponential(float y0, float ysol, float beta, float s, float uy) const {
  float result = 0.1 * exp((-y0 + ysol) * beta) * (1.f - exp(-s * uy*beta)) / (uy*beta);
  return result;
}

float Fog::tr(const Ray& ray, float yFloor) const {
  float int_ext;
  if(m_isUniform){
    int_ext = m_beta * ray.tmax;
  }
  else
  {
    int_ext = int_exponential(ray.orig.z, yFloor, m_beta, ray.tmax, ray.dir.z);
  }
  return exp(-int_ext);
}

color3 Fog::sample(const Ray &ray, Scene* scene, KdTree* tree, float yFloor) const {
  vec3 Lv(0.f);
  
  float randt, probat;
  float clamped_t = min(10000.f, ray.tmax);
  if(m_uniformSampleRay){
    randt = uniform01(engine) * clamped_t;
    probat = 1. / clamped_t;
  }
  else
  {
    float alpha = 5. / clamped_t;
    do{
      randt = -log(uniform01(engine)) / alpha;
    } while(randt > clamped_t);
    float normalization = 1.f / alpha * (1.f - exp(-alpha * clamped_t));
    probat = exp(-alpha * randt) / normalization;
  }

  float int_ext_partiel;
  if(m_isUniform){
    int_ext_partiel = m_beta * randt;
  }
  else
  {
    int_ext_partiel = int_exponential(ray.orig.z, yFloor, m_beta, randt, ray.dir.z);
  }
  vec3 randP = ray.orig + randt * ray.dir;

  vec3 randDir;
  float probaDir;
  auto sphereL = scene->objects[scene->objects.size()-1];
  vec3 axePO = normalize(randP - sphereL->geom.sphere.center);
  vec3 ptA;

  bool is_uniform;
  if(uniform01(engine) < m_pUniform){
    randDir = random_uniform();
    is_uniform = true;
  }
  else
  {
    vec3 dirA = random_dir(axePO);
    ptA = dirA * sphereL->geom.sphere.radius + sphereL->geom.sphere.center;
    randDir = normalize(ptA - randP);
    is_uniform = false;
  }

  float phase_f;
  float k = 0.4;
  switch(m_phase){
    case 0: 
      phase_f = 0.3f / (4.f * M_PI);
      break;
    case 1:
      phase_f = (1. - (k*k)) / (4.f * Pi * (1. + k * dot(randDir, -ray.dir)));
      break;
    case 2:
      phase_f = 3./(16. * Pi) * (1.f + sqr(dot(randDir, ray.dir)));
      break;
    default:
      phase_f = 0.3f / (4.f * M_PI);
      break;
  }


  Ray L_Ray;
  L_Ray.hasDifferentials = false;
  rayInit(&L_Ray, randP, randDir, ray.pixel, 0, 100000, ray.depth+1);
  Intersection interL;
  color3 L = trace_ray(scene, &L_Ray, tree, &interL);

  float V;
  if(is_uniform){
    V = 1;
  }
  else
  {
    float d_light2 = length_sq(ptA - randP);
    if(interL.hit && L_Ray.tmax * L_Ray.tmax < d_light2*0.9){
      V = 0;
    }else
    {
      V = 1;
    }
  }

  if(V == 0) {
    Lv = vec3(0);
  }
  else
  {
    vec3 interN = interL.normal;
    vec3 interP = interL.position;

    float pdf_uniform = 1.f / (4.f * Pi);
    float J = dot(interN, -randDir) / glm::length_sq(interP - randt);
    float pdf_light = (interL.hit && !isBlack(interL.mat->m_emission)) ? (dot(normalize(interP - sphereL->geom.sphere.center), axePO) / (Pi * sqr(sphereL->geom.sphere.radius)) / J) : 0.f;

    probaDir = m_pUniform * pdf_uniform + (1.f - m_pUniform) * pdf_light;

    float ext;
    if(m_isUniform){
      ext = m_beta;
    }
    else
    {
      ext = 0.1 * exp(-m_beta * (randP.z - yFloor));
    }

    Lv = L * phase_f * ext * exp(-int_ext_partiel) / (probat * probaDir); 
  }
  return Lv;
}
