#pragma once

#include "defines.h"
#include "ray.h"
#include "raytracer.h"

#include <random>
#include <iostream>

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

        vec3 forward = normalize(lookat - lookfrom); 
        vec3 x = cross(forward, vup);
        vec3 up = normalize(cross(x, forward));
        h = tan(vfov / 2.f * (float)M_PI / 180.f) * (2.0f);
        w = h * width / float(height);

        vec3 zAxis = -forward; 
        yAxis = up;
        xAxis = cross(yAxis, zAxis);

        topLeft = lookfrom - zAxis + yAxis * h /2.f - xAxis * w / 2.f;

        pos = lookfrom;
        fov = vfov;
        aspect = aspect_ratio;

        imgWidth = width;
        imgHeight = height;


        dXPixel = xAxis * w / (float)width;
        dXPixel = yAxis * h / (float)height;

    }

    ~SimpleCamera() {}

    void get_ray(float s, float t, Ray *r, vec2 pixel) const override {
      vec3 d = topLeft + (s + 0.5f) * w / imgWidth * xAxis - (t + 0.5f) * h / imgHeight * yAxis;
      vec3 ray_dir = d - pos;
      rayInit(r, pos, normalize(ray_dir), pixel);
      r->dox = vec3(0.f);
      r->doy = vec3(0.f);
      //r->ddx = (dot(d,d) * dXPixel - dot(d, dXPixel) * d) / glm::pow(dot(d,d),1.5f);
      //r->ddy = (dot(d,d) * dYPixel - dot(d, dYPixel) * d) / glm::pow(dot(d,d),1.5f);
      r->ddx = dXPixel;
      r->ddy = dYPixel;
    }

    vec3 pos;
    float fov;
    float aspect;

    float h, w;
    vec3 yAxis;
    vec3 xAxis;
    vec3 topLeft;
    vec3 dXPixel;
    vec3 dYPixel;
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
