#include "Light.h"
#include "intersection.hpp"
#include "raytracer.h"
#include "scene.h"

bool is_shadowed(vec3 lightPosition, vec3 point, Scene* scene, KdTree* tree){
  Intersection temp_inter;
  Ray ray;
  vec3 dir = (lightPosition - point);
  dir = dir / length(dir);
  rayInit(&ray, point + (acne_eps * dir), dir, 0.f, distance(point + (acne_eps * dir), lightPosition));
  ray.shadow = true;
  if(intersectKdTree(scene, tree, &ray, &temp_inter)){
    return true;
  }
  return false;
}


PointLight::PointLight(vec3 position, color3 color){
  m_position = position;
  m_color = color;
  m_samples.push_back(position);
}

PointLight::~PointLight(){

}

float PointLight::intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection){
  return is_shadowed(m_position, intersection->position, scene, tree) ? 0.f : 1.f;
}

AreaLight::AreaLight(vec3 corner, vec3 full_uvec, int usteps, vec3 full_vvec, int vsteps, vec3 color){
    m_corner = corner;
    m_usteps = usteps;
    m_vsteps = vsteps;
    uvec = full_uvec / float(usteps);
    vvec = full_vvec / float(vsteps);
    m_color = color;  
    nbSamples = usteps * vsteps;
    m_position = corner + (full_uvec / 2.f) + (full_vvec / 2.f);

    for(int v = 0; v < m_vsteps; v++) {
      for(int u = 0; u < m_usteps; u++){
          m_samples.push_back(pointOnLight(u, v));
      }
    }
}

#include <random>
static std::minstd_rand engine(time(NULL));
static std::uniform_real_distribution<float> m_unifDistributionRand{0.f, 1.0f};

float AreaLight::intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection){
  float intensity = 0.0f;
  for(auto& sample : m_samples) {
    if(!is_shadowed(sample, intersection->position, scene, tree)){
      intensity += 1.0;
    }
  }
  return intensity / float(nbSamples);
}
point3 AreaLight::pointOnLight(float u, float v){
  return m_corner +
    uvec * (u + m_unifDistributionRand(engine)) +
    vvec * (v + m_unifDistributionRand(engine));
}

AreaLight::~AreaLight() {

}

void AreaLight::setup(Scene *scene){
  addObject(scene, m_t1);
  addObject(scene, m_t2);
}
