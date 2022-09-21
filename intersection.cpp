#include "intersection.hpp"

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
      intersection->u = fmod(intersection->position.x, 1.0f);
      intersection->v = fmod(intersection->position.z, 1.0f);
      ray->tmax = t;
      return true;
    }
  }

  return false;
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
      vec3 normal = normalize(intersection->position - obj->geom.sphere.center);
      intersection->isOutside = dot(ray->dir, normal) < 0;
      //intersection->normal = intersection->isOutside ? normal : -normal;
      intersection->normal = normal;
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
    vec3 normal = normalize(intersection->position - obj->geom.sphere.center);
    intersection->isOutside = dot(ray->dir, normal) < 0;
    //intersection->normal = intersection->isOutside ? normal : -normal;
    intersection->normal = normal;
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
    intersection->isOutside = dot(ray->dir, obj->geom.triangle.normal) < 0;
    //intersection->normal = intersection->isOutside ? obj->geom.triangle.normal : -obj->geom.triangle.normal;
    intersection->u = u;
    intersection->v = v;

    vec3 normal = obj->geom.triangle.n2 * u + obj->geom.triangle.n3 * v + obj->geom.triangle.n1 * (1.f - u - v);

    intersection->normal = normal;
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
      if (intersectSphere(ray, temp, scene->objects[i]))
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
