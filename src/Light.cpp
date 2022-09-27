//#include "kdtree.h"
//#include "intersection.hpp"
//#define GLM_ENABLE_EXPERIMENTAL
//#include <glm/gtx/norm.hpp>
//
//point2 UniformSampleTriangle(const point2 &u) {
//    float su0 = std::sqrt(u[0]);
//    return point2(1 - su0, u[1] * su0);
//}
//
//#include "kdtree.h"
//
//    float Pdf(const Intersection &ref, const vec3& wi, Scene* scene, KdTree* tree, float area){
//      // Intersect sample ray with area light geometry
//      Ray ray;
//      rayInit(&ray, ref.position + (wi*acne_eps), wi);
//      Intersection isectLight;
//      // Ignore any alpha textures used for trimming the shape when performing
//      // this intersection. Hack for the "San Miguel" scene, where this is used
//      // to make an invisible area light.
//      if (!intersectKdTree(scene, tree, &ray, &isectLight)) return 0;
//
//      // Convert light sample weight to solid angle measure
//      float pdf = glm::distance2(ref.position, isectLight.position) /
//                  (abs(dot(isectLight.normal, -wi) * area));
//      if (std::isinf(pdf)) pdf = 0.f;
//      return pdf;
//    }
//
//    Intersection Sample(const point2& u, float *pdf, float area, Object* obj){
//      point2 b = UniformSampleTriangle(u);
//      // Get triangle vertices in _p0_, _p1_, and _p2_
//      const point3 &p0 = obj->geom.triangle.p1;
//      const point3 &p1 = obj->geom.triangle.p2;
//      const point3 &p2 = obj->geom.triangle.p3;
//      Intersection it;
//      it.position = b[0] * p0 + b[1] * p1 + (1 - b[0] - b[1]) * p2;
//      // Compute surface normal for sampled point on triangle
//      it.normal = normalize(vec3(cross(p1 - p0, p2 - p0)));
//      // Ensure correct orientation of the geometric normal; follow the same
//      // approach as was used in Triangle::Intersect().
//
//      // Compute error bounds for sampled point on triangle
//      *pdf = 1 / area;
//      return it;
//    }
//
//Intersection Sample(const Intersection &ref, const point2 &u,
//                          float *pdf, float area, Object* obj) {
//    Intersection intr = Sample(u, pdf, area, obj);
//    vec3 wi = intr.position - ref.position;
//    if ((wi.length() * wi.length()) == 0)
//        *pdf = 0;
//    else {
//        wi = normalize(wi);
//        // Convert from area measure, as returned by the Sample() call
//        // above, to solid angle measure.
//        *pdf *= distance2(ref.position, intr.position) / abs(dot(intr.normal, -wi));
//        if (std::isinf(*pdf)) *pdf = 0.f;
//    }
//    return intr;
//}
//
//
//DiffuseAreaLight::DiffuseAreaLight(const color3 &Lemit, int nSamples,
//                                   Object* shape,
//                                   bool twoSided)
//    : AreaLight(),
//      Lemit(Lemit),
//      shape(shape),
//      twoSided(twoSided),
//      area(0.5 * cross(shape->geom.triangle.p2 - shape->geom.triangle.p1, shape->geom.triangle.p3 - shape->geom.triangle.p1).length()) {
//}
//
//color3 DiffuseAreaLight::Power() const {
//    return (twoSided ? 2.f : 1.f) * Lemit * area * float(M_PI);
//}
//
//color3 DiffuseAreaLight::Sample_Li(const Intersection &ref, const point2 &u,
//                                     vec3 *wi, float *pdf ) const {
//    Intersection pShape = Sample(ref, u, pdf, area, shape );
//    auto l = (pShape.position - ref.position).length();
//    if (*pdf == 0 || (l*l) == 0) {
//        *pdf = 0;
//        return color3(0.f);
//    }
//    *wi = normalize(pShape.position - ref.position);
//    return L(pShape, -*wi);
//}
//
//float DiffuseAreaLight::Pdf_Li(const Intersection &ref,
//                               const vec3 &wi, Scene* scene, KdTree* tree) const {
//    return Pdf(ref, wi, scene, tree, area);
//}
//
//point2 ConcentricSampleDisk(const point2 &u) {
//    // Map uniform random numbers to $[-1,1]^2$
//    point2 uOffset = 2.f * u - vec2(1, 1);
//
//    // Handle degeneracy at the origin
//    if (uOffset.x == 0 && uOffset.y == 0) return point2(0, 0);
//
//    // Apply concentric mapping to point
//    float theta, r;
//    float PiOver4 = M_PI/4;
//    if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
//        r = uOffset.x;
//        theta = PiOver4 * (uOffset.y / uOffset.x);
//    } else {
//        r = uOffset.y;
//        theta = (M_PI/2) - PiOver4 * (uOffset.x / uOffset.y);
//    }
//    return r * point2(std::cos(theta), std::sin(theta));
//}
//
//inline vec3 CosineSampleHemisphere(const point2 &u) {
//    point2 d = ConcentricSampleDisk(u);
//    float z = std::sqrt(std::max((float)0, 1 - d.x * d.x - d.y * d.y));
//    return vec3(d.x, d.y, z);
//}
//
//inline float CosineHemispherePdf(float cosTheta) { return cosTheta * (1/M_PI); }
//
//inline void CoordinateSystem(const vec3 &v1, vec3 *v2,
//                             vec3 *v3) {
//    if (std::abs(v1.x) > std::abs(v1.y))
//        *v2 = vec3(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
//    else
//        *v2 = vec3(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
//    *v3 = cross(v1, *v2);
//}
//
//color3 DiffuseAreaLight::Sample_Le(const point2 &u1, const point2 &u2,
//                                     float time, Ray *ray, vec3 *nLight,
//                                     float *pdfPos, float *pdfDir) const {
//    // Sample a point on the area light's _Shape_, _pShape_
//    Intersection pShape = Sample(u1, pdfPos, area, shape);
//    *nLight = pShape.normal;
//
//    // Sample a cosine-weighted outgoing direction _w_ for area light
//    vec3 w;
//    if (twoSided) {
//        point2 u = u2;
//        // Choose a side to sample and then remap u[0] to [0,1] before
//        // applying cosine-weighted hemisphere sampling for the chosen side.
//        if (u[0] < .5) {
//            u[0] = std::min(u[0] * 2, 1-acne_eps);
//            w = CosineSampleHemisphere(u);
//        } else {
//            u[0] = std::min((u[0] - .5f) * 2, 1-acne_eps);
//            w = CosineSampleHemisphere(u);
//            w.z *= -1;
//        }
//        *pdfDir = 0.5f * CosineHemispherePdf(std::abs(w.z));
//    } else {
//        w = CosineSampleHemisphere(u2);
//        *pdfDir = CosineHemispherePdf(w.z);
//    }
//
//    vec3 v1, v2, n(pShape.normal);
//    CoordinateSystem(n, &v1, &v2);
//    w = w.x * v1 + w.y * v2 + w.z * n;
//    rayInit(ray, pShape.position +(w*acne_eps), w);
//    return L(pShape, w);
//}
//
//void DiffuseAreaLight::Pdf_Le(const Ray &ray, const vec3 &n, float *pdfPos,
//                              float *pdfDir) const {
//    Intersection it;
//    it.position = ray.orig;
//    it.normal = n;
//    *pdfPos = 1 / area;
//    *pdfDir = twoSided ? (.5 * CosineHemispherePdf(abs(dot(n, ray.dir))))
//                       : CosineHemispherePdf(dot(n, ray.dir));
//}

#include "Light.h"
#include "intersection.hpp"
#include "raytracer.h"

bool is_shadowed(vec3 lightPosition, vec3 point, Scene* scene, KdTree* tree){
  Intersection temp_inter;
  Ray ray;
  vec3 dir = (lightPosition - point);
  dir = dir / length(dir);
  rayInit(&ray, point + (acne_eps * dir), dir, 0.f, distance(point + (acne_eps * dir), lightPosition));
  ray.shadow = true;
  if(intersectKdTree(scene, tree, &ray, &temp_inter)){
    return true;
  }
  return false;
}


PointLight::PointLight(vec3 position, color3 color){
  m_position = position;
  m_color = color;
}

PointLight::~PointLight(){

}

float PointLight::intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection){
  return is_shadowed(m_position, intersection->position, scene, tree) ? 0.f : 1.f;
}

AreaLight::AreaLight(vec3 corner, vec3 full_uvec, int usteps, vec3 full_vvec, int vsteps, vec3 color){
    m_corner = corner;
    m_usteps = usteps;
    m_vsteps = vsteps;
    uvec = full_uvec / float(usteps);
    vvec = full_vvec / float(vsteps);
    m_color = color;  
    samples = usteps * vsteps;
    m_position = corner + (full_uvec / 2.f) + (full_vvec / 2.f);
}

#include <random>
static std::minstd_rand engine(time(NULL));
static std::uniform_real_distribution<float> m_unifDistributionRand{0.f, 1.0f};

float AreaLight::intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection){
  float intensity = 0.0f;
  for(int v = 0; v < m_vsteps; v++) {
    for(int u = 0; u < m_usteps; u++){
        auto lightPosition = pointOnLight(u, v);
        if(!is_shadowed(lightPosition, intersection->position, scene, tree)){
          intensity += 1.0;
        }
    }
  }
  return intensity / float(samples);
}
point3 AreaLight::pointOnLight(float u, float v){
  return m_corner +
         uvec * (u + m_unifDistributionRand(engine)) +
         vvec * (v + m_unifDistributionRand(engine));
}

AreaLight::~AreaLight() {

}
