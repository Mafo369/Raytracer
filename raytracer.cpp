
#include <chrono>
#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>
#include <cmath>
#include <string.h>

#include <omp.h>
#include <glm/gtc/epsilon.hpp>

#include <cstdlib>
#include <limits>

#include <iostream>

#include "camera.h"
#include <glm/glm.hpp>
#include <random>

namespace glm {
// not really needed, but it makes it easier to follow the book...
template <int N, typename T, qualifier P> T length_sq(const vec<N, T, P> &x) { return dot(x, x); }
} // namespace glm

/// acne_eps is a small constant used to prevent acne when computing
/// intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;


static std::minstd_rand engine((omp_get_thread_num()+1));
static std::uniform_real_distribution<float> m_unifDistributionRand{-1.0f, 1.0f};

vec3 sphereRand(){
  while(true){
    vec3 p(m_unifDistributionRand(engine), m_unifDistributionRand(engine), 
          m_unifDistributionRand(engine));
    if (glm::length_sq(p) >= 1) continue;
    return p;
  }
}

vec3 mDiskRand(){
  while(true){
    auto p = vec3(m_unifDistributionRand(engine), m_unifDistributionRand(engine), 0);
    if (glm::length_sq(p) >= 1) continue;
    return p;
  }
}

static float reflectance(float cosine, float ref_idx) {
  // Use Schlick's approximation for reflectance.
  float r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
  r0 = r0 * r0;
  return r0 + (1.0f - r0) * std::pow((1.0f - cosine), 5.0f);
}
        
bool scatter(Ray *r_in, Intersection rec, color3 &attenuation, Ray *scattered) {
  
  if(rec.mat->type == LAMBERTIAN){
    vec3 scatter_direction = rec.normal + sphereRand();
            
    // Catch degenerate scatter direction
    //if (near_zero(scatter_direction))
    if (glm::any(glm::epsilonEqual(scatter_direction, vec3{0, 0, 0}, std::numeric_limits<float>::epsilon())))
      scatter_direction = rec.normal;

    rayInit(scattered, rec.position, scatter_direction);
    attenuation = rec.mat->diffuseColor;
  } else if (rec.mat->type == METAL){
      vec3 reflected = glm::reflect(glm::normalize(r_in->dir), rec.normal);
      rayInit(scattered, rec.position, reflected + rec.mat->fuzz * sphereRand());
      attenuation = rec.mat->diffuseColor;
      return (glm::dot(scattered->dir, rec.normal) > 0);
  } else if (rec.mat->type == DIELECTRIC){
      attenuation = color3(1.0f, 1.0f, 1.0f);
      float refraction_ratio = rec.front_face ? (1.0f/rec.mat->IOR) : rec.mat->IOR;

      vec3 unit_direction = glm::normalize(r_in->dir);
      float cos_theta = std::min(glm::dot(-unit_direction, rec.normal), 1.0f);
      float sin_theta = std::sqrt(1.0f - cos_theta*cos_theta);

      bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
      vec3 direction;
      static std::minstd_rand m_rnGenerator{};
      static std::uniform_real_distribution<float> m_unifDistribution{0.0f, 1.0f};
      if (cannot_refract || reflectance(cos_theta, refraction_ratio) > m_unifDistribution(m_rnGenerator))
        direction = glm::reflect(unit_direction, rec.normal);
      else
        direction = glm::refract(unit_direction, rec.normal, refraction_ratio);

      rayInit(scattered, rec.position, direction);
  }
  
  return true;
}


bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj)
{

  //! \todo : compute intersection of the ray and the plane object

  vec3 n = obj->geom.plane.normal;
  float dn = dot(ray->dir, n);

  if (dn == 0.0)
  {
    //Pas des solutions
    return false;
  }
  else
  {
    float t = -(dot(ray->orig, n) + obj->geom.plane.dist) / dn;
    if (t >= ray->tmin && t <= ray->tmax)
    {
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = &(obj->mat);
      intersection->normal = n;
      ray->tmax = t;
      return true;
    }
  }
  return false;
}
    

inline void set_face_normal(Ray *r, const vec3& outward_normal, Intersection *rec) {
    rec->front_face = glm::dot(r->dir, outward_normal) < 0;
    rec->normal = rec->front_face ? outward_normal :-outward_normal;
}

bool hit(Ray *r, float t_min, float t_max,Intersection *rec, Object *obj) {
    vec3 oc = r->orig - obj->geom.sphere.center;
    auto a = glm::length_sq(r->dir);
    auto half_b = glm::dot(oc, r->dir);
    float radius = obj->geom.sphere.radius;
    auto c = glm::length_sq(oc) - radius*radius;

    auto discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }
    
    r->tmax = root;
    rec->position = r->orig + (r->tmax * r->dir);
    //vec3 outward_normal = (rec->position - obj->geom.sphere.center) / obj->geom.sphere.radius;
    vec3 outward_normal = glm::normalize(rec->position - obj->geom.sphere.center);

    if (radius < 0)
      outward_normal *= -1.0f;

    set_face_normal(r, outward_normal, rec);
    rec->mat = &(obj->mat);
    
    return true;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj)
{

  //! \todo : compute intersection of the ray and the sphere object

  vec3 oc = ray->orig - obj->geom.sphere.center;
  float b = 2 * (dot(ray->dir, oc));
  float c = dot(oc, oc) - (obj->geom.sphere.radius * obj->geom.sphere.radius);

  float delta = b * b - 4 * c;

  if (delta == 0)
  {
    //Une solution
    float t = -b / 2;
    if (t >= ray->tmin && t <= ray->tmax)
    {
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = &(obj->mat);
      intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
      ray->tmax = t;
      return true;
    }
  }
  else if (delta > 0)
  {
    //Deux solutions
    float t1 = (-b + sqrtf(delta)) / 2;
    float t2 = (-b - sqrtf(delta)) / 2;
    float t;
    if (t1 >= ray->tmin && t1 <= ray->tmax && t2 >= ray->tmin && t2 <= ray->tmax)
    {
      t = std::min(t1, t2);
    }
    else if (t1 >= ray->tmin && t1 <= ray->tmax)
    {
      t = t1;
    }
    else if (t2 >= ray->tmin && t2 <= ray->tmax)
    {
      t = t2;
    }
    else
    {
      return false;
    }
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = &(obj->mat);
    intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
    ray->tmax = t;
    return true;
  }
  else
  {
    //Pas de solutions -> pas d'intersection
  }
  return false;
}

//Moller-Trumbore Ray Triangle Intersection
bool intersectTriangle(Ray *ray, Intersection *intersection, Object *obj)
{
  vec3 v1v2 = obj->geom.triangle.p2 - obj->geom.triangle.p1;
  vec3 v1v3 = obj->geom.triangle.p3 - obj->geom.triangle.p1;

  vec3 cross_rayDir_v1v3 = cross(ray->dir, v1v3);

  float det = dot(v1v2, cross_rayDir_v1v3);

  if (det > -acne_eps && det < acne_eps)
    return false;

  float inv_det = 1.f / det;

  vec3 o_minus_p1 = ray->orig - obj->geom.triangle.p1;

  float u = dot(o_minus_p1, cross_rayDir_v1v3) * inv_det;

  if (u < 0.f || u > 1.f)
    return false;

  vec3 cross_oMinusp1_v1v2 = cross(o_minus_p1, v1v2);

  float v = dot(ray->dir, cross_oMinusp1_v1v2) * inv_det;

  if (v < 0.f || u + v > 1.f)
    return false;

  float t = dot(v1v3, cross_oMinusp1_v1v2) * inv_det;

  if (t >= ray->tmin && t <= ray->tmax)
  {
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = &(obj->mat);
    intersection->normal = normalize(cross(v1v2, v1v3));
    ray->tmax = t;
    return true;
  }
  return false;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection)
{
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();

  //!\todo loop on each object of the scene to compute intersection

  float dist;

  for (size_t i = 0; i < objectCount; i++)
  {
    Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
    if (scene->objects[i]->geom.type == PLANE)
    {
      if (intersectPlane(ray, temp, scene->objects[i]))
      {
        float temp_dist = ray->tmax;
        if (hasIntersection)
        {
          if (temp_dist < dist)
          {
            dist = temp_dist;
            *intersection = *temp;
          }
        }
        else
        {
          hasIntersection = true;
          *intersection = *temp;
          dist = temp_dist;
        }
      }
    }
    else if (scene->objects[i]->geom.type == SPHERE)
    {
      if (hit(ray, 0.001, std::numeric_limits<float>::infinity(),temp, scene->objects[i]))
      {
        float temp_dist = ray->tmax;
        if (hasIntersection)
        {
          if (temp_dist < dist)
          {
            dist = temp_dist;
            *intersection = *temp;
          }
        }
        else
        {
          hasIntersection = true;
          *intersection = *temp;
          dist = temp_dist;
        }
      }
    }
    else if (scene->objects[i]->geom.type == TRIANGLE)
    {
      if (intersectTriangle(ray, temp, scene->objects[i]))
      {
        float temp_dist = ray->tmax;
        if (hasIntersection)
        {
          if (temp_dist < dist)
          {
            dist = temp_dist;
            *intersection = *temp;
          }
        }
        else
        {
          hasIntersection = true;
          *intersection = *temp;
          dist = temp_dist;
        }
      }
    }
    free(temp);
  }
  return hasIntersection;
}

color3 trace_ray(Scene *scene, Ray *r, KdTree *tree) {

  Intersection rec;

  if(r->depth > 5)
    return color3(0,0,0);

  if (intersectKdTree(scene, tree, r, &rec)){
  //if(intersectScene(scene, r, &rec)){
    Ray scattered;
    color3 attenuation;
    if (scatter(r, rec, attenuation, &scattered)){
      scattered.depth =  r->depth + 1;
      return attenuation * trace_ray(scene, &scattered, tree);
    }
    return color3(0,0,0);
  }
  vec3 unit_direction = glm::normalize(r->dir);
  auto t = 0.5*(unit_direction.y + 1.0);
  auto op_t = 1.0 - t;
  return color3(op_t*1.0, op_t*1.0, op_t*1.0) + color3(t*0.5, t*0.7, t*1.0);
}
  
// Utility Functions

inline float degrees_to_radians(float degrees) {
  return degrees * M_PI / 180.0f;
}

void renderImage(Image *img, Scene *scene)
{

  // rng stuff

  auto samples_per_pixel = 50;
  
  //float dist_to_focus = glm::length(scene->cam.position-scene->cam.lookat);
  float dist_to_focus = 10.0;
  auto aperture = 0.1;
      

  //std::mt19937 m_rnGenerator{};
  std::minstd_rand engine((omp_get_thread_num()+1));
  std::uniform_real_distribution<float> m_unifDistribution{0.0f, 1.0f};
  
  camera cam(scene->cam.position, scene->cam.lookat, scene->cam.up, scene->cam.fov, scene->cam.aspect, aperture, dist_to_focus);

  KdTree *tree = initKdTree(scene);

  auto startTime = std::chrono::system_clock::now();
  //printf("End building tree\n");

  //! \todo initialize KdTree
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
      for (int s = 0; s < samples_per_pixel; ++s) {
        auto u = (i + m_unifDistribution(engine)) / (img->width-1);
        auto v = (j + m_unifDistribution(engine)) / (img->height-1);
        Ray r;
        cam.get_ray(u, v, &r);
        pixel_color += trace_ray(scene, &r, tree);
      }
      // Divide the color by the number of samples and gamma-correct for gamma=2.0.
      pixel_color /= samples_per_pixel;
      pixel_color = glm::sqrt(pixel_color);
      pixel_color = glm::clamp(pixel_color, 0.0f, 1.0f);

      *ptr = pixel_color;
    }
  }
  auto stopTime = std::chrono::system_clock::now();
  std::cout << "Rendering took "<< std::chrono::duration_cast<std::chrono::duration<double>>(
                    stopTime - startTime).count() << "s" << std::endl;
}
