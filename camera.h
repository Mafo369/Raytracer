#ifndef CAMERA_H
#define CAMERA_H

#include "defines.h"
#include "ray.h"
#include "raytracer.h"

class camera {
    public:
        camera(
            point3 lookfrom,
            point3 lookat,
            vec3   vup,
            double vfov, // vertical field-of-view in degrees
            double aspect_ratio,
            double aperture,
            double focus_dist 
        ) {
            auto theta = degrees_to_radians(vfov);
            auto h = tan(theta/2);
            auto viewport_height = 2.0 * h;
            auto viewport_width = aspect_ratio * viewport_height;

            w = unit_vector(lookfrom - lookat);
            u = unit_vector(cross(vup, w));
            v = cross(w, u);

            origin = lookfrom;
            horizontal = (float)viewport_width * u;
            vertical = (float)viewport_height * v;
            lower_left_corner = origin - horizontal/2.f - vertical/2.f - w;
            //horizontal = (float)(focus_dist * viewport_width) * u;
            //vertical = (float)(focus_dist * viewport_height) * v;
            //lower_left_corner = origin - horizontal/2.f - vertical/2.f - (float)focus_dist*w;

            //lens_radius = aperture / 2;
        }

        camera() {};

        void get_ray(double s, double t, Ray *r) const {
            //vec3 rd = (float)lens_radius * random_in_unit_disk();
            //vec3 offset = u * rd.x + v * rd.y;
            
            //rayInit(r, origin + offset, normalize(lower_left_corner + (float)s*horizontal + (float)t*vertical - origin - offset));
            rayInit(r, origin,normalize(lower_left_corner + (float)s*horizontal + (float)t*vertical - origin));
        }

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        double lens_radius;
};
#endif
