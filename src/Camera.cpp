#include "Camera.h"
#include "ray.h"
#include "sampling/sampling.h"

SimpleCamera::SimpleCamera(point3 lookfrom, point3 lookat, vec3 vup, float vfov, float aspect_ratio, 
             int width, int height, float lensRadius, float focalDistance) {
    pos = lookfrom;
    fov = vfov;
    aspect = aspect_ratio;
    m_lensRadius = lensRadius;
    m_focalDistance = focalDistance;

    imgWidth = width;
    imgHeight = height;

    dir = lookat;
    dir -= pos;

    dir = normalize(dir);
    xDir = cross(dir, vup);
    up = normalize(cross(xDir, dir));

    left =  normalize(cross(up, dir));
    float radfov = fov * (float)M_PI / 180.0;
    l = m_focalDistance;
    h = tan(radfov * .5f) * (2.0f * l);
    w = h * imgWidth / (float)imgHeight;

    point3 B = pos + (l * dir) + (h / 2.0f * up);
    nearPlaneTopLeft = B + w / 2.0f * (left);
    dXPixel = -left * (w / (float)imgWidth);
    dYPixel = -up * (w / (float)imgWidth);
}

SimpleCamera::~SimpleCamera() {}

Ray SimpleCamera::get_ray(float s, float t, float uLens, float vLens, vec2 pixel, bool hasDifferentials) const {
  Ray r;
  if(m_lensRadius > 0) {
    point2 pLens = m_lensRadius * ConcentricSampleDisk(vec2(uLens, vLens));
    vec3 offset = xDir * pLens.x + up * pLens.y;
    vec3 d = (nearPlaneTopLeft + s * dXPixel + t  * dYPixel) - pos - offset;
    vec3 dnorm = normalize(d);
    r = Ray(pos + offset, dnorm, pixel);

    if(hasDifferentials){
      r.dox = pos + offset;
      r.doy = pos + offset;
      r.ddx = (nearPlaneTopLeft + (s+1.f) * dXPixel + t * dYPixel) - pos - offset;
      r.ddy = (nearPlaneTopLeft + s * dXPixel + (t+1.f) * dYPixel) - pos - offset;
    }
  }
  else{
    vec3 d = (nearPlaneTopLeft + s * dXPixel + t  * dYPixel) - pos;
    vec3 dnorm = normalize(d);
    r = Ray(pos, dnorm, pixel);
    if(hasDifferentials){
      r.dox = pos;
      r.doy = pos;
      r.ddx = (nearPlaneTopLeft + (s+1.f) * dXPixel + t * dYPixel) - pos;
      r.ddy = (nearPlaneTopLeft + s * dXPixel + (t+1.f) * dYPixel) - pos;
    }
  }
  
  r.hasDifferentials = hasDifferentials;
  return r;
}

CameraFOV::CameraFOV(point3 lookfrom, point3 lookat, vec3 vup, float vfov, float width, float height,
       float aperture, float focus_dist ) {
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

CameraFOV::~CameraFOV() {};

vec3 CameraFOV::mDiskRand() const {
  while(true){
    auto p = vec3(m_unifDistributionRand(engine), m_unifDistributionRand(engine), 0);
    if (glm::length_sq(p) >= 1) continue;
    return p;
  }
}

Ray CameraFOV::get_ray(float s, float t, float uLens, float vLens, vec2 pixel, bool hasDifferentials) const {
    glm::vec2 rd = lens_radius * mDiskRand();
    vec3 offset = u * rd.x + v * rd.y;
    vec3 dir = lower_left_corner + s*horizontal + t*vertical - origin - offset;

    return Ray(origin + offset, normalize(dir), pixel);
}
