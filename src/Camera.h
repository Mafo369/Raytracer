#pragma once

#include "defines.h"
#include "ray.h"
#include "raytracer.h"

#include <random>
#include <iostream>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

static std::minstd_rand engine(time(NULL));
static std::uniform_real_distribution<float> m_unifDistributionRand{-1.f, 1.0f};

class Camera {
  public:
    Camera() = default;
    virtual ~Camera() = default;

    virtual void get_ray(float s, float t, Ray* r, vec2 pixel) const = 0;
    float imgWidth;
    float imgHeight;
};

class SimpleCamera : public Camera {
  public:
    SimpleCamera(point3 lookfrom, point3 lookat, vec3 vup, float vfov, float aspect_ratio, 
                 int width, int height) : Camera() {
        pos = lookfrom;
        fov = vfov;
        aspect = aspect_ratio;

        imgWidth = width;
        imgHeight = height;

        auto theta = degrees_to_radians(vfov);
        auto h = std::tan(theta/2.f);
        auto viewport_height = 2.0f * h;
        auto viewport_width = aspect_ratio * viewport_height;

        auto w = normalize(lookfrom - lookat);
        u = normalize(cross(vup, w));
        v = cross(w, u);

        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        lower_left_corner = pos - horizontal/2.f - vertical/2.f - w;

        dXPixel = viewport_width / (float)width;
        dYPixel = viewport_height / (float)height;
    }

    ~SimpleCamera() {}

    void get_ray(float s, float t, Ray *r, vec2 pixel) const override {
      //vec3 d = (nearPlaneTopLeft + (s + 0.5f) * dXPixel + (t + .5f) * dYPixel) ;
      //vec3 ray_dir = d - pos;
      vec3 view = lower_left_corner - pos;
      vec3 d = view + s * horizontal + t*vertical;
      rayInit(r, pos, normalize(d), pixel);
      r->dox = vec3(0.f);
      r->doy = vec3(0.f);
      r->ddx = (dot(d,d) * horizontal - dot(d, horizontal) * d) / glm::pow(dot(d,d),1.5f);
      r->ddy = (dot(d,d) * vertical - dot(d, vertical) * d) / glm::pow(dot(d,d),1.5f);
      //std::cout << glm::to_string(r->ddx) << std::endl;
      r->dXPixel = horizontal;
      r->dYPixel = vertical;
    }

    vec3 pos;
    float fov;
    float aspect;
    vec3 horizontal;
    vec3 vertical;
    vec3 u;
    vec3 v;

    vec3 lower_left_corner;
    float dXPixel;
    float dYPixel;
};

class CameraFOV : public Camera {
    public:
        CameraFOV(point3 lookfrom, point3 lookat, vec3 vup, float vfov, float width, float height,
               float aperture, float focus_dist ) : Camera() {
            float aspect_ratio = width / height;
            auto theta = degrees_to_radians(vfov);
            auto h = std::tan(theta/2.f);
            auto viewport_height = 2.0f * h;
            auto viewport_width = aspect_ratio * viewport_height;
            lens_radius = aperture / 2.0f;
            imgWidth = width;
            imgHeight = height;

            w = glm::normalize(lookfrom - lookat);
            u = glm::normalize(glm::cross(vup, w));
            v = glm::cross(w, u);

            origin = lookfrom;
            horizontal = focus_dist * viewport_width * u;
            vertical = focus_dist * viewport_height * v;
            lower_left_corner = origin - horizontal/2.0f - vertical/2.0f - focus_dist*w;
        }

        ~CameraFOV() {};

        vec3 mDiskRand() const {
          while(true){
            auto p = vec3(m_unifDistributionRand(engine), m_unifDistributionRand(engine), 0);
            if (glm::length_sq(p) >= 1) continue;
            return p;
          }
        }

        void get_ray(float s, float t, Ray *r, vec2 pixel) const override {
            glm::vec2 rd = lens_radius * mDiskRand();
            vec3 offset = u * rd.x + v * rd.y;

            vec3 dir = lower_left_corner + s*horizontal + t*vertical - origin - offset;
            rayInit(r, origin + offset, normalize(dir), pixel);
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        float lens_radius;

};
