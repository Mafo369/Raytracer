#pragma once

#include "defines.h"
#include "ray.h"
#include "sampling/sampling.h"

#include <random>
#include <iostream>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

static std::uniform_real_distribution<float> m_unifDistributionRand{-1.f, 1.0f};

class Camera {
  public:
    Camera() = default;
    virtual ~Camera() = default;

    virtual Ray get_ray(float s, float t, float uLens, float vLens, vec2 pixel, bool hasDifferentials = true) const = 0;
    float imgWidth;
    float imgHeight;
};

class SimpleCamera : public Camera {
  public:
    SimpleCamera(point3 lookfrom, point3 lookat, vec3 vup, float vfov, float aspect_ratio, 
                 int width, int height, float lensRadius, float focalDistance);

    ~SimpleCamera() override;

    Ray get_ray(float s, float t, float uLens, float vLens, vec2 pixel, bool hasDifferentials = true) const override;

    vec3 pos = vec3(0,0,0);
    vec3 dir = vec3(0,0,-1);
    vec3 up = vec3(0,1,0);
    vec3 left;
    float fov;
    float aspect;
    float l;
    float h;
    float w;
    vec3 nearPlaneTopLeft;
    vec3 dXPixel;
    vec3 dYPixel;
    float m_lensRadius;
    float m_focalDistance;
    vec3 xDir;
};

class CameraFOV : public Camera {
    public:
        CameraFOV(point3 lookfrom, point3 lookat, vec3 vup, float vfov, float width, float height,
               float aperture, float focus_dist );

        ~CameraFOV() override;

        vec3 mDiskRand() const;

        Ray get_ray(float s, float t, float uLens, float vLens, vec2 pixel, bool hasDifferentials = true) const override;

    private:
        point3 origin;
        point3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v, w;
        float lens_radius;

};
