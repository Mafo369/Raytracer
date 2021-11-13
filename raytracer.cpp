
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

inline float random_float() {
    // Returns a random real in [0,1).
    return (float)rand() / (float)(RAND_MAX + 1.0);
}

inline float random_float(float min, float max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_float();
}

inline vec3 unit_vector(vec3 v) {
    float l = length(v);
    return vec3(v.x / l, v.y / l, v.z / l);
}
        
inline vec3 random(float min, float max) {
    return vec3(random_float(min,max), random_float(min,max), random_float(min,max));
}
        
float length_squared(vec3 e) {
  return e.x*e.x + e.y*e.y + e.z*e.z;
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
  return (fabs(e.x) < s) && (fabs(e.y) < s) && (fabs(e.z) < s);
}
        
bool scatter(Ray *r_in, Intersection rec, color3 &attenuation, Ray *scattered) {
  
  if(rec.mat->type == LAMBERTIAN){
    auto scatter_direction = rec.normal + random_unit_vector();
            
    // Catch degenerate scatter direction
    if (near_zero(scatter_direction))
      scatter_direction = rec.normal;

    rayInit(scattered, rec.position, scatter_direction);
    attenuation = rec.mat->diffuseColor;
  } else if (rec.mat->type == METAL){
      vec3 reflected = reflect(unit_vector(r_in->dir), rec.normal);
      rayInit(scattered, rec.position, reflected + rec.mat->fuzz*random_in_unit_sphere(), 0, 10000, r_in->depth+1);
      attenuation = rec.mat->diffuseColor;
      return (dot(scattered->dir, rec.normal) > 0);
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
    bool front_face = dot(r->dir, outward_normal) < 0;
    rec->normal = front_face ? outward_normal :-outward_normal;
}

bool hit(Ray *r, float t_min, float t_max,Intersection *rec, Object *obj) {
    vec3 oc = r->orig - obj->geom.sphere.center;
    auto a = length_squared(r->dir);
    auto half_b = dot(oc, r->dir);
    auto c = length_squared(oc) - obj->geom.sphere.radius*obj->geom.sphere.radius;

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

  if(r->depth > 4)
    return color3(0,0,0);

  //if (intersectKdTree(scene, tree, r, &rec)){
  if(intersectScene(scene, r, &rec)){
    Ray scattered;
    color3 attenuation;
    if (scatter(r, rec, attenuation, &scattered)){
      scattered.depth =  r->depth + 1;
      return attenuation * trace_ray(scene, &scattered, tree);
    }
    return color3(0,0,0);
  }
  vec3 unit_direction = unit_vector(r->dir);
  auto t = 0.5*(unit_direction.y + 1.0);
  return color3((1.0-t)*1.0, (1.0-t)*1.0, (1.0-t)*1.0) + color3(t*0.5, t*0.7, t*1.0);
}
  
// Utility Functions

inline float degrees_to_radians(float degrees) {
  return degrees * M_PI / 180.0;
}

vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_float(-1,1), random_float(-1,1), 0);
        if (length_squared(p) >= 1) continue;
        return p;
    }
}
        
inline float clamp(float x, float min, float max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}
float length(vec3 v) {
  return sqrt(length_squared(v));
}

void write_color(std::ostream &out, color3 pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x;
    auto g = pixel_color.y;
    auto b = pixel_color.z;

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / (float)samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

void renderImage(Image *img, Scene *scene)
{

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  
  auto samples_per_pixel = 100;
  
  float dist_to_focus = length(scene->cam.position-scene->cam.lookat);
  //std::cout << dist_to_focus << std::endl;
  auto aperture = 0.1;
  
  camera cam(scene->cam.position, scene->cam.lookat, scene->cam.up, scene->cam.fov, scene->cam.aspect, aperture, dist_to_focus);

  KdTree *tree = initKdTree(scene);

  //printf("End building tree\n");

  //! \todo initialize KdTree
  /*
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
        auto u = (i + random_float()) / (img->width-1);
        auto v = (j + random_float()) / (img->height-1);
        Ray r;
        cam.get_ray(u, v, &r);
        //get_ray(u1, v1, &r, lens_radius, u, v, origin, lower_left_corner, horizontal, vertical);
        pixel_color += trace_ray(scene, &r, tree);
      }
      auto r = pixel_color.x;
      auto g = pixel_color.y;
      auto b = pixel_color.z;

      // Divide the color by the number of samples and gamma-correct for gamma=2.0.
      auto scale = 1.f / samples_per_pixel;
      r = sqrt(scale * r);
      g = sqrt(scale * g);
      b = sqrt(scale * b);

      color3 def_color(clamp(r, 0.0, 0.999),clamp(g, 0.0, 0.999),clamp(b, 0.0, 0.999)); 
      color3 *ptr = getPixelPtr(img, i, j);
      *ptr = def_color;
      //color3 *ptr = getPixelPtr(img, i, j);
      // *ptr = trace_ray_multisampling(scene, tree, i, j, dx, dy, ray_delta_x, ray_delta_y);
    }
  }
  */

  const int image_width = 800;
  const int image_height = 600;

  std::cout << "P3\n" << 800 << ' ' << 600 << "\n255\n";

  for(int j = image_height-1; j >= 0; --j) {
    std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
    for(int i = 0; i < image_width; ++i) {
      color3 pixel_color(0, 0, 0);
#pragma omp parallel for
      for (int s = 0; s < samples_per_pixel; ++s) {
        auto u = (i + random_float()) / (image_width-1);
        auto v = (j + random_float()) / (image_height-1);
        Ray r;
        cam.get_ray(u, v, &r);
        pixel_color += trace_ray(scene, &r, tree);
      }
      write_color(std::cout, pixel_color, samples_per_pixel);
    }
  }
  std::cerr << "\nDone.\n";

}
