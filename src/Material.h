#pragma once

#include "defines.h"
#include <vector>
#include "raytracer.h"
//#include "Light.h"
#include "sampling/sampler.h"

enum MatType {
  DIFFUSE=0,
  SPECULAR=1,
  TRANSPARENT=2
};

typedef struct s_BrdfData {
  float LdotH;
  float NdotH;
  float VdotH;
  float LdotN;
  float VdotN;
  vec3 F;
  color3 specularF0;
  color3 diffuseReflectance;
  vec2 uv;
  float roughness;
  float alpha;
  float alphaSq;

} BrdfData;

// Specifies minimal reflectance for dielectrics (when metalness is zero)
// Nothing has lower reflectance than 2%, but we use 4% to have consistent results with UE4, Frostbite, et al.
#define MIN_DIELECTRICS_F0 0.04f

// Schlick's approximation to Fresnel term
// f90 should be 1.0, except for the trick used by Schuler (see 'shadowedF90' function)
inline vec3 evalFresnelSchlick(vec3 f0, float f90, float NdotS)
{
	return f0 + (f90 - f0) * pow(clamp(1.0f - NdotS, 0.0f, 1.0f), 5.0f);
}

inline vec3 evalFresnel(vec3 f0, float f90, float NdotS)
{
	// Default is Schlick's approximation
	return evalFresnelSchlick(f0, f90, NdotS);
}

inline vec3 baseColorToSpecularF0(const vec3& baseColor, const float& metalness) {
	return lerp(vec3(MIN_DIELECTRICS_F0, MIN_DIELECTRICS_F0, MIN_DIELECTRICS_F0), baseColor, metalness);
}

inline vec3 baseColorToDiffuseReflectance(vec3 baseColor, float metalness)
{
	return baseColor * (1.0f - metalness);
}

inline float luminance(vec3 rgb)
{
	return dot(rgb, vec3(0.2126f, 0.7152f, 0.0722f));
}

// Attenuates F90 for very low F0 values
// Source: "An efficient and Physically Plausible Real-Time Shading Model" in ShaderX7 by Schuler
// Also see section "Overbright highlights" in Hoffman's 2010 "Crafting Physically Motivated Shading Models for Game Development" for discussion
// IMPORTANT: Note that when F0 is calculated using metalness, it's value is never less than MIN_DIELECTRICS_F0, and therefore,
// this adjustment has no effect. To be effective, F0 must be authored separately, or calculated in different way. See main text for discussion.
inline float shadowedF90(vec3 F0) {
	// This scaler value is somewhat arbitrary, Schuler used 60 in his article. In here, we derive it from MIN_DIELECTRICS_F0 so
	// that it takes effect for any reflectance lower than least reflective dielectrics
	//const float t = 60.0f;
	const float t = (1.0f / MIN_DIELECTRICS_F0);
	return min(1.0f, t * luminance(F0));
}

bool computeBrdfData(BrdfData& data, const vec3& v, const vec3& wi, const vec3& n, const vec2& uv, 
                     const color3& baseColor, const float& metalness, const float& roughness);

class Material {
  public:
    Material() = default;
    virtual ~Material() = default;
    virtual color3 f(const vec3& wo, const vec3& wi, const vec3& n) = 0;
    virtual color3 sample_f(vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type) = 0;
    virtual float pdf(const vec3& wo, const vec3& wi, const vec3& n) = 0;
    virtual color3 shade(Intersection* intersection, vec3 v, Light* light, float intensity) = 0;
    virtual color3 textureColor(float u, float v, int face) = 0;
    virtual color3 ambientColor(Ray* ray, Intersection* intersection, color3 lightColor) = 0;
    virtual color3 scatterColor(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection) = 0;
    texture* m_texture = nullptr;
    color3 m_emission = color3(0,0,0);
    color3 m_albedo = color3(0,0,0);

    float m_IOR;	//! Index of refraction (for dielectric)
};

color3 specularReflect(Ray* ray, Intersection* intersection, Scene* scene, KdTree* tree, Sampler* sampler);
color3 specularTransmission(Ray* ray, Intersection* intersection, Scene* scene, KdTree* tree, Sampler* sampler);
