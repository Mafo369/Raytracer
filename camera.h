#ifndef CAMERA_H
#define CAMERA_H

#include "defines.h"
#include "ray.h"
#include "raytracer.h"

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
            auto h = tan(theta/2.f);
            auto viewport_height = 2.f * h;
            auto viewport_width = aspect_ratio * viewport_height;

            std::cerr << "lookfrom:" << lookfrom.x << " " << lookfrom.y << " " << lookfrom.z << std::endl;
            std::cerr << "lookat:" << lookat.x << " " << lookat.y << " " << lookat.z << std::endl;

            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);

            origin = lookfrom;
            horizontal = focus_dist * viewport_width * u;
            vertical = focus_dist * viewport_height * v;
            lower_left_corner = origin - horizontal/2.f - vertical/2.f - focus_dist*w;
            lens_radius = aperture / 2.f;

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
            vec3 rd = lens_radius * random_in_unit_disk();
            vec3 offset = u * rd.x + v * rd.y;

            rayInit(r, origin + offset, lower_left_corner + s*horizontal + t*vertical - origin - offset);
            //rayInit(r, origin,normalize(lower_left_corner + (float)s*horizontal + (float)t*vertical - origin));
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
