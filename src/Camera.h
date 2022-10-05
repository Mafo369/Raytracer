#pragma once

#include "defines.h"
#include "ray.h"
#include "raytracer.h"

#include <random>
#include <iostream>

namespace glm {
// not really needed, but it makes it easier to follow the book...
template <int N, typename T, qualifier P> T length_sq(const vec<N, T, P> &x) { return dot(x, x); }
} // namespace glm

// Utility Functions
inline float degrees_to_radians(float degrees) {
  return degrees * M_PI / 180.0f;
}

static std::minstd_rand engine(time(NULL));
static std::uniform_real_distribution<float> m_unifDistributionRand{-1.f, 1.0f};


class Camera {
    public:
        Camera(
            point3 lookfrom,
            point3 lookat,
            vec3   vup,
            float vfov, // vertical field-of-view in degrees
            float aspect_ratio,
            float aperture,
            float focus_dist 
        ) {
            auto theta = degrees_to_radians(vfov);
            auto h = std::tan(theta/2.f);
            auto viewport_height = 2.0f * h;
            auto viewport_width = aspect_ratio * viewport_height;
            lens_radius = aperture / 2.0f;

            w = glm::normalize(lookfrom - lookat);
            u = glm::normalize(glm::cross(vup, w));
            v = glm::cross(w, u);

            origin = lookfrom;
            horizontal = focus_dist * viewport_width * u;
            vertical = focus_dist * viewport_height * v;
            lower_left_corner = origin - horizontal/2.0f - vertical/2.0f - focus_dist*w;
        }

        ~Camera() {};

        vec3 mDiskRand() const {
          while(true){
            auto p = vec3(m_unifDistributionRand(engine), m_unifDistributionRand(engine), 0);
            if (glm::length_sq(p) >= 1) continue;
            return p;
          }
        }

        void get_ray(float s, float t, Ray *r) const {
            glm::vec2 rd = lens_radius * mDiskRand();
            vec3 offset = u * rd.x + v * rd.y;

            rayInit(r, origin + offset, normalize(lower_left_corner + s*horizontal + t*vertical - origin - offset));
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        float lens_radius;
};
