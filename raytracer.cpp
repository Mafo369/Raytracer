
#include "image.h"
#include "kdtree.h"
#include "ray.h"
#include "raytracer.h"
#include "scene_types.h"
#include <stdio.h>
#include <cmath>
#include <string.h>

#include <glm/gtc/epsilon.hpp>

/// acne_eps is a small constant used to prevent acne when computing
/// intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj) {
  
  //! \todo : compute intersection of the ray and the plane object
  
  vec3 n = obj->geom.plane.normal;
  float dn = dot(ray->dir, n);

  if(dn == 0.0){
    //Pas des solutions
    return false; 
  }else{
    float t = -(dot(ray->orig, n)+ obj->geom.plane.dist)/dn;
    if(t >= ray->tmin && t <= ray->tmax){
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = &(obj->mat);
      intersection->normal = n;
      return true;
    }
  }  

  return false;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj) {

  //! \todo : compute intersection of the ray and the sphere object

  vec3 oc = ray->orig - obj->geom.sphere.center;
  float b = 2*(dot(ray->dir, oc));
  float c = dot(oc, oc) - (obj->geom.sphere.radius * obj->geom.sphere.radius);

  float delta = b*b - 4*c;

  if(delta == 0){
    //Une solution
    float t = -b/2;
    if(t >= ray->tmin && t <= ray->tmax){
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = &(obj->mat);
      intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
      return true;
    }
  }else if(delta > 0){
    //Deux solutions
    float t1 = (-b + sqrtf(delta))/2;
    float t2 = (-b - sqrtf(delta))/2;
    float t;
    if(t1 >= ray->tmin && t1 <= ray->tmax && t2 >= ray->tmin && t2 <= ray->tmax){
      t = std::min(t1, t2);
    }else if(t1 >= ray->tmin && t1 <= ray->tmax){
      t = t1;
    }else if(t2 >= ray->tmin && t2 <= ray->tmax){
      t = t2;
    }else{
      return false;
    }
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = &(obj->mat);
    intersection->normal = normalize(intersection->position - obj->geom.sphere.center);
    return true;
  }else{
    //Pas de solutions -> pas d'intersection
  }
  return false;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();

//!\todo loop on each object of the scene to compute intersection

  float dist;

  for(size_t i = 0; i<objectCount; i++){
    Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
    if(scene->objects[i]->geom.type == PLANE){
      if(intersectPlane(ray, temp, scene->objects[i])){
        float temp_dist = distance(ray->orig, temp->position);
        if(hasIntersection){
          if(temp_dist < dist){
            dist = temp_dist;
            memcpy(intersection, temp, sizeof(Intersection));
          }
        }
        else{
          hasIntersection = true;
          memcpy(intersection, temp, sizeof(Intersection));
          dist = temp_dist; 
        }
      }
    }
    else if(scene->objects[i]->geom.type == SPHERE){
      if(intersectSphere(ray, temp, scene->objects[i])){
        float temp_dist = distance(ray->orig, temp->position);
        if(hasIntersection){
          if(temp_dist < dist){
            dist = temp_dist;
            memcpy(intersection, temp, sizeof(Intersection));
          }
        }
        else{
          hasIntersection = true;
          memcpy(intersection, temp, sizeof(Intersection));
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
float RDM_Beckmann(float NdotH, float alpha) {


  //! \todo compute Beckmann normal distribution
  return 0.5f;

}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {

  //! \todo compute Fresnel term
  return 0.5f;

}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {

  //! \todo compute G1 term of the Smith fonction
  return 0.5f;

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha) {

  //! \todo the Smith fonction
  return 0.5f;

}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, Material *m) {

  //! \todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G
  //! = RDM_Smith
  return color3(.5f);

}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {

  //! \todo compute diffuse component of the bsdf
  return color3(.5f);

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                Material *m) {

  //! \todo compute bsdf diffuse and specular term
  return color3(0.f);

}



color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat) {
  color3 ret = color3(0.f);

//! \todo compute bsdf, return the shaded color taking into account the
//! lightcolor

  float ln = dot(l, n);
  float pi = M_PI;

  if(ln > 0)
    ret = (mat->diffuseColor/pi) * ln * lc;

  return ret;
}

//! if tree is not null, use intersectKdTree to compute the intersection instead
//! of intersect scene

color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree) {

  color3 ret = color3(0, 0, 0);
  Intersection intersection;

  if(intersectScene(scene, ray, &intersection)){
    //ret = (0.5f * intersection.normal) + 0.5f;
    
    size_t lightsCount = scene->lights.size();
    
    for(size_t i=0; i<lightsCount; i++){
      vec3 n = ray->dir * -1.0f;
      vec3 lp = scene->lights[i]->position - intersection.position;
      vec3 l = lp/length(lp);
      
      Ray *ray_ombre = (Ray *)malloc(sizeof(Ray));
      ray_ombre->orig = intersection.position + (acne_eps * l);
      ray_ombre->dir = l;
      ray_ombre->tmax = ray->tmax;
      ray_ombre->tmin = ray->tmin;

      Intersection temp_inter;
      if(!intersectScene(scene, ray_ombre, &temp_inter)){
        ret += shade(intersection.normal, n, l, scene->lights[i]->color, intersection.mat);
      }
      else
        ret += 0.f;

      free(ray_ombre);
    }
  }else{
    ret = scene->skyColor;
  }


  return ret;
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing
  //! and kdtree initializaion
  float aspect = 1.f / scene->cam.aspect;

  KdTree *tree = NULL;


//! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f);   //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) *
                     aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x =
      (0.5f - img->width * 0.5f) / (img->width * 0.5f) * scene->cam.xdir;


  for (size_t j = 0; j < img->height; j++) {
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
    for (size_t i = 0; i < img->width; i++) {
      color3 *ptr = getPixelPtr(img, i, j);
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y +
                     float(i) * dx + float(j) * dy;

      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      *ptr = trace_ray(scene, &rx, tree);

    }
  }
}
