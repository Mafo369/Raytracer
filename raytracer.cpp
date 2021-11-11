
#include "glm/fwd.hpp"
#include "glm/geometric.hpp"
#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>
#include <cmath>
#include <string.h>

#include <glm/gtc/epsilon.hpp>

#include <cstdlib>
#include <limits>

#include <iostream>

#include "camera.h"

/// acne_eps is a small constant used to prevent acne when computing
/// intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}

inline vec3 unit_vector(vec3 v) {
    int l = v.length();
    return vec3(v.x / l, v.y / l, v.z / l);
}
        
inline vec3 random(double min, double max) {
    return vec3(random_double(min,max), random_double(min,max), random_double(min,max));
}
        
double length_squared(vec3 e) {
  return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
}

vec3 random_in_unit_sphere() {
    while (true) {
        auto p = random(-1,1);
        if (length_squared(p) >= 1) continue;
        return p;
    }
}

vec3 random_unit_vector() {
    return unit_vector(random_in_unit_sphere());
}
        
bool near_zero(vec3 e) {
  // Return true if the vector is close to zero in all dimensions.
  const auto s = 1e-8;
  return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
}
        
bool scatter(Ray *r_in, Intersection rec, color3 &attenuation, Ray *scattered) {
  auto scatter_direction = rec.normal + random_unit_vector();
            
  // Catch degenerate scatter direction
  if (near_zero(scatter_direction))
    scatter_direction = rec.normal;

  rayInit(scattered, rec.position, normalize(scatter_direction), 0, 10000, r_in->depth+1);
  attenuation = rec.mat->diffuseColor;
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
    bool front_face = dot(r->dir, outward_normal) < 0;
    rec->normal = front_face ? outward_normal :-outward_normal;
}

bool hit(Ray *r, double t_min, double t_max,Intersection *rec, Object *obj) {
    vec3 oc = r->orig - obj->geom.sphere.center;
    auto a = length_squared(r->dir);
    auto half_b = dot(oc, r->dir);
    auto c = length_squared(oc) - obj->geom.sphere.radius*obj->geom.sphere.radius;

    auto discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < r->tmin || r->tmax < root) {
        root = (-half_b + sqrtd) / a;
        if (root < r->tmin || r->tmax < root)
            return false;
    }

    r->tmax = root;
    rec->position = r->orig + (r->tmax * r->dir);
    vec3 outward_normal = (rec->position - obj->geom.sphere.center) / obj->geom.sphere.radius;
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
      if (hit(ray, 0.001, std::numeric_limits<double>::infinity(),temp, scene->objects[i]))
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

/* ---------------------------------------------------------------------------
 */
/*
 *	The following functions are coded from Cook-Torrance bsdf model
 *description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba
 *renderer)
 */

// Shadowing and masking function. Linked with the NDF. Here, Smith function,
// suitable for Beckmann NDF
float RDM_chiplus(float c) { return (c > 0.f) ? 1.f : 0.f; }

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha)
{

  //! \todo compute Beckmann normal distribution

  float phi = RDM_chiplus(NdotH);

  float cos2 = NdotH * NdotH;
  float tan2 = (1 - (cos2)) / cos2;
  float alpha2 = alpha * alpha;
  float e = exp(-tan2 / (alpha2));
  float pi_alpha_cos = M_PI * alpha2 * (cos2 * cos2);

  float d = phi * (e / pi_alpha_cos);

  return d;
}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR)
{

  //! \todo compute Fresnel term
  if (LdotH < 0)
  {
    LdotH = -LdotH;
  }

  float n1_n2 = (extIOR / intIOR) * (extIOR / intIOR);
  float sin2_t = n1_n2 * (1 - (LdotH * LdotH));
  if (sin2_t >= 1.0f)
  {
    return 1.f;
  }
  float cos_t = sqrtf(1 - sin2_t);

  float ncos_minus_it = (extIOR * LdotH - intIOR * cos_t);
  float ncos_minus2_it = ncos_minus_it * ncos_minus_it;
  float ncos_plus_it = (extIOR * LdotH + intIOR * cos_t);
  float ncos_plus2_it = ncos_plus_it * ncos_plus_it;

  float ncos_minus_ti = (extIOR * cos_t - intIOR * LdotH);
  float ncos_minus2_ti = ncos_minus_ti * ncos_minus_ti;
  float ncos_plus_ti = (extIOR * cos_t + intIOR * LdotH);
  float ncos_plus2_ti = ncos_plus_ti * ncos_plus_ti;

  float rs = ncos_minus2_it / ncos_plus2_it;
  float rp = ncos_minus2_ti / ncos_plus2_ti;

  float f = (rs + rp) / 2;

  return f;
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha)
{

  //! \todo compute G1 term of the Smith fonction

  float tanx = sqrtf(1 - (DdotN * DdotN)) / DdotN;
  float b = 1 / (alpha * tanx);

  float k = DdotH / DdotN;

  float phi_k = RDM_chiplus(k);

  if (b < 1.6f)
  {
    float b2 = b * b;
    float fraction = (3.535 * b + 2.181 * b2) / (1 + 2.276 * b + 2.577 * b2);
    float g1 = phi_k * fraction;
    return g1;
  }

  return phi_k;
}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha)
{
  return RDM_G1(LdotH, LdotN, alpha) * RDM_G1(VdotH, VdotN, alpha);
}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, Material *m)
{

  //! \todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G
  //! = RDM_Smith

  float d = RDM_Beckmann(NdotH, m->roughness);
  float f = RDM_Fresnel(LdotH, 1.f, m->IOR);
  float g = RDM_Smith(LdotH, LdotN, VdotH, VdotN, m->roughness);

  return color3(m->specularColor * ((d * f * g) / (4 * LdotN * VdotN)));
}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m)
{

  float pi = M_PI;
  return color3(m->diffuseColor / pi);
}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m)
{

  //! \todo compute bsdf diffuse and specular term
  return color3(RDM_bsdf_d(m) + RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, m));
}

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat)
{
  color3 ret = color3(0.f);

  //! \todo compute bsdf, return the shaded color taking into account the
  //! lightcolor

  float ln = dot(l, n);

  if (ln > 0)
  {
    vec3 vl = v + l;
    vec3 h = vl / length(vl);

    float LdotH = dot(l, h);
    float NdotH = dot(n, h);
    float VdotH = dot(v, h);
    float LdotN = dot(l, n);
    float VdotN = dot(v, n);
    ret = lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat) * LdotN;
  }

  return ret;
}

//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scene

color3 ray_color(Scene *scene, Ray *ray, KdTree *tree)
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
      vec3 lp = scene->lights[i]->position - intersection.position;
      vec3 l = lp / length(lp);

      Ray *ray_ombre = new Ray();
      rayInit(ray_ombre, intersection.position + (acne_eps * l), l, 0.f, distance(intersection.position + (acne_eps * l), scene->lights[i]->position));

      Intersection temp_inter;
      if (!intersectKdTree(scene, tree, ray_ombre, &temp_inter))
      {
        ret += shade(intersection.normal, v, l, scene->lights[i]->color, intersection.mat);
      }
      else
      {
        ret += color3(0.f);
      }

      free(ray_ombre);
    }

    if (ret.r > 1.f && ret.g > 1.f && ret.b > 1.f && ray->depth > 0) // Si contribution maximale -> on arrete
      return color3(1.f);

    vec3 r = reflect(ray->dir, intersection.normal);
    Ray *ray_ref = (Ray *)malloc(sizeof(Ray));
    rayInit(ray_ref, intersection.position + (acne_eps * r), r, 0, 100000, ray->depth + 1);

    color3 cr = ray_color(scene, ray_ref, tree);
    float LdotH = dot(ray_ref->dir, intersection.normal);
    float f = RDM_Fresnel(LdotH, 1.f, intersection.mat->IOR);

    ret += (f * cr * intersection.mat->specularColor);

    free(ray_ref);
  }
  else
  {
    ret = scene->skyColor;
  }

  return ret;
}


color3 trace_ray(Scene *scene, Ray *r, KdTree *tree) {

  Intersection rec;

  if(r->depth > 2)
    return color3(0,0,0);

  if (intersectKdTree(scene, tree, r, &rec)){
  //if(intersectScene(scene, r, &rec)){
    Ray *scattered = (Ray *)malloc(sizeof(Ray));
    color3 attenuation;
    if (scatter(r, rec, attenuation, scattered))
      return attenuation * trace_ray(scene, scattered, tree);
    return color3(0,0,0);
  }
  vec3 unit_direction = r->dir;
  auto t = 0.5*(unit_direction.y + 1.0);
  return color3((1.0-t)*1.0, (1.0-t)*1.0, (1.0-t)*1.0) + color3(t*0.5, t*0.7, t*1.0);
}
  
// Utility Functions

inline double degrees_to_radians(double degrees) {
  return degrees * M_PI / 180.0;
}

vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (length_squared(p) >= 1) continue;
        return p;
    }
}
        
inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

void renderImage(Image *img, Scene *scene)
{

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  
  auto samples_per_pixel = 100;
  
  auto dist_to_focus = (scene->cam.position-scene->cam.lookat).length();
  auto aperture = 0.005;
  
  camera cam(scene->cam.position, scene->cam.lookat, scene->cam.up, scene->cam.fov, scene->cam.aspect, aperture, dist_to_focus);

  KdTree *tree = initKdTree(scene);

  printf("End building tree\n");

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
      for (int s = 0; s < samples_per_pixel; ++s) {
        auto u = (i + random_double()) / (img->width-1);
        auto v = (j + random_double()) / (img->height-1);
        Ray r;
        cam.get_ray(u, v, &r);
        //get_ray(u1, v1, &r, lens_radius, u, v, origin, lower_left_corner, horizontal, vertical);
        pixel_color += trace_ray(scene, &r, tree);
      }
      auto r = pixel_color.x;
      auto g = pixel_color.y;
      auto b = pixel_color.z;

      // Divide the color by the number of samples and gamma-correct for gamma=2.0.
      auto scale = 1.0 / samples_per_pixel;
      r = sqrt(scale * r);
      g = sqrt(scale * g);
      b = sqrt(scale * b);

      color3 def_color(clamp(r, 0.0, 0.999),clamp(g, 0.0, 0.999),clamp(b, 0.0, 0.999)); 
      color3 *ptr = getPixelPtr(img, i, j);
      *ptr = def_color;
      //color3 *ptr = getPixelPtr(img, i, j);
      //*ptr = trace_ray_multisampling(scene, tree, i, j, dx, dy, ray_delta_x, ray_delta_y);
    }
  }
}
