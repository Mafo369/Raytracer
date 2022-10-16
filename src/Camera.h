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
        zdir = normalize(lookat - lookfrom); 
        xdir = cross(zdir, vup);
        ydir = normalize(cross(xdir, zdir));

        pos = lookfrom;
        fov = vfov;
        aspect = aspect_ratio;

        imgWidth = width;
        imgHeight = height;

        vec3 left = normalize(cross(ydir, zdir));
        float radfov = fov * float(M_PI) / 180.f;
        float h = tan(radfov * .5f) * (2.0f);
        float w = h * width / float(height);

        point3 B = pos + zdir + (h / 2.0f * ydir);
        nearPlaneTopLeft = B + w / 2.0f * left;
        dXPixel = -left * (w / float(width));
        dYPixel = -ydir * (w / float(width));

    }

    ~SimpleCamera() {}

    void get_ray(float s, float t, Ray *r, vec2 pixel) const override {
      vec3 ray_dir = (nearPlaneTopLeft + (s + 0.5f) * dXPixel + (t + .5f) * dYPixel) - pos;
      rayInit(r, pos, normalize(ray_dir), pixel);
    }

    vec3 pos;
    float fov;
    float aspect;
    vec3 zdir;
    vec3 xdir;
    vec3 ydir;

    vec3 nearPlaneTopLeft;
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
