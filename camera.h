#ifndef CAMERA_H
#define CAMERA_H

#include "defines.h"
#include "ray.h"
#include "raytracer.h"

#include <random>
#include <glm/gtc/random.hpp>
#include <iostream>

class camera {
    public:
        camera(
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

            std::cerr << "lookfrom:" << lookfrom.x << " " << lookfrom.y << " " << lookfrom.z << std::endl;
            std::cerr << "lookat:" << lookat.x << " " << lookat.y << " " << lookat.z << std::endl;

            w = glm::normalize(lookfrom - lookat);
            u = glm::normalize(glm::cross(vup, w));
            v = glm::cross(w, u);

            origin = lookfrom;
            horizontal = focus_dist * viewport_width * u;
            vertical = focus_dist * viewport_height * v;
            lower_left_corner = origin - horizontal/2.0f - vertical/2.0f - focus_dist*w;

            std::cerr << "w:" << w.x << " " << w.y << " " << w.z << std::endl;
            std::cerr << "u:" << u.x << " " << u.y << " " << u.z << std::endl;
            std::cerr << "v:" << v.x << " " << v.y << " " << v.z << std::endl;

            std::cerr << "o:" << origin.x << " " << origin.y << " " << origin.z << std::endl;
            std::cerr << "h:" << horizontal.x << " " << horizontal.y << " " << horizontal.z << std::endl;
            std::cerr << "v:" << vertical.x << " " << vertical.y << " " << vertical.z << std::endl;
            std::cerr << "llc:" << lower_left_corner.x << " " << lower_left_corner.y << " " << lower_left_corner.z << std::endl;
            std::cerr << "lr:" << lens_radius << std::endl;
        }

        camera() {};

        void get_ray(float s, float t, Ray *r) const {
            glm::vec2 rd = glm::diskRand(lens_radius);
            vec3 offset = u * rd.x + v * rd.y;

            rayInit(r, origin + offset, lower_left_corner + s*horizontal + t*vertical - origin - offset);
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        float lens_radius;
};
#endif
