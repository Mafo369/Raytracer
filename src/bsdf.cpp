#include "bsdf.hpp"

#include <cmath>
#include <glm/glm.hpp>
#include "Object.h"

/* ---------------------------------------------------------------------------
 */
/*
 *	The following functions are coded from Cook-Torrance bsdf model
 *description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba
 *renderer)
 */

// Shadowing and masking function. Linked with the NDF. Here, Smith function,
// suitable for Beckmann NDF
float RDM_chiplus(float c) { return (c > 0.f) ? 1.f : 0.f; }

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha)
{

  //! \todo compute Beckmann normal distribution

  float phi = RDM_chiplus(NdotH);

  float cos2 = NdotH * NdotH;
  float tan2 = (1 - (cos2)) / cos2;
  float alpha2 = alpha * alpha;
  float e = exp(-tan2 / (alpha2));
  float pi_alpha_cos = M_PI * alpha2 * (cos2 * cos2);

  float d = phi * (e / pi_alpha_cos);

  return d;
}

// Fresnel term computation. Implantation of the exact computation. we can use
// the Schlick approximation
// LdotH : Light . Half
//float RDM_Fresnel(float LdotH, float extIOR, float intIOR)
//{
//
//  //! \todo compute Fresnel term
//  if (LdotH < 0)
//  {
//    LdotH = -LdotH;
//  }
//
//  float n1_n2 = (extIOR / intIOR) * (extIOR / intIOR);
//  float sin2_t = n1_n2 * (1 - (LdotH * LdotH));
//  if (sin2_t >= 1.0f)
//  {
//    return 1.0f;
//  }
//  float cos_t = sqrtf(1 - sin2_t);
//
//  float ncos_minus_it = (extIOR * LdotH - intIOR * cos_t);
//  float ncos_minus2_it = ncos_minus_it * ncos_minus_it;
//  float ncos_plus_it = (extIOR * LdotH + intIOR * cos_t);
//  float ncos_plus2_it = ncos_plus_it * ncos_plus_it;
//
//  float ncos_minus_ti = (extIOR * cos_t - intIOR * LdotH);
//  float ncos_minus2_ti = ncos_minus_ti * ncos_minus_ti;
//  float ncos_plus_ti = (extIOR * cos_t + intIOR * LdotH);
//  float ncos_plus2_ti = ncos_plus_ti * ncos_plus_ti;
//
//  float rs = ncos_minus2_it / ncos_plus2_it;
//  float rp = ncos_minus2_ti / ncos_plus2_ti;
//
//  float f = (rs + rp) / 2;
//
//  return f;
//}

// extIOR = etaI
// intIOR = etaT
float RDM_Fresnel(float LdotH, float extIOR, float intIOR){
  LdotH = std::clamp(LdotH, -1.f, 1.f);
  bool entering = LdotH > 0.f;
  if(!entering){
    std::swap(extIOR, intIOR);
    LdotH = std::abs(LdotH);
  }
  float sinThetaI = std::sqrt(std::max(0.f, 1.f - LdotH * LdotH));
  float sinThetaT = extIOR / intIOR * sinThetaI;
  if(sinThetaT >= 1)
    return 1.f;

  float cosThetaT = std::sqrt(std::max(0.f, 1.f - sinThetaT * sinThetaT));

  float Rparl = ((intIOR * LdotH) - (extIOR * cosThetaT)) /
                ((intIOR * LdotH) + (extIOR * cosThetaT));
  float Rperp = ((extIOR * LdotH) - (intIOR * cosThetaT)) /
                ((extIOR * LdotH) + (intIOR * cosThetaT));
  return (Rparl * Rparl + Rperp * Rperp) / 2.f;
}

// extIOR and intIOR are inversed here compared to fresnel
// extIOR = etaT
// intIOR = etaI
float schlick(float VdotN, float extIOR, float intIOR){
  float cos = VdotN;

  if( intIOR > extIOR ){
      float n = intIOR / extIOR;
      float sin2_t = (n*n) * (1.0f - (cos*cos)); 
      if( sin2_t > 1.0f )
        return 1.0f;
      
      float cos_t = sqrtf(1.0f - sin2_t);

      cos = cos_t;
  }

  float r0 = ((intIOR - extIOR) / (intIOR + extIOR));
  float r02 = r0 * r0;
  return r02 + (1.f - r02) * powf(1.f - cos, 5);
}

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha)
{

  //! \todo compute G1 term of the Smith fonction

  float tanx = sqrtf(1.f - (DdotN * DdotN)) / DdotN;
  float b = 1.f / (alpha * tanx);

  float k = DdotH / DdotN;

  float phi_k = RDM_chiplus(k);

  if (b < 1.6f)
  {
    float b2 = b * b;
    float fraction = (3.535f * b + 2.181f * b2) / (1.f + 2.276f * b + 2.577f * b2);
    float g1 = phi_k * fraction;
    return g1;
  }

  return phi_k;
}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN,
                float alpha)
{
  return RDM_G1(LdotH, LdotN, alpha) * RDM_G1(VdotH, VdotN, alpha);
}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, float roughness, float IOR, color3 specularColor)
{

  //! \todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G
  //! = RDM_Smith

  float d = RDM_Beckmann(NdotH, roughness);
  float f = RDM_Fresnel(LdotH, 1.f, IOR);
  //float f = schlick(LdotH, m->IOR, 1.f);
  float g = RDM_Smith(LdotH, LdotN, VdotH, VdotN, roughness);

  return color3(specularColor * ((d * f * g) / (4.f * LdotN * VdotN)));
}

color3 RDM_btdf(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, float extIOR, float intIOR, float roughness, color3 specularColor)
{

  float d = RDM_Beckmann(NdotH, roughness);
  float f = RDM_Fresnel(LdotH, extIOR, intIOR);
  //float f = schlick(LdotH, extIOR, intIOR);
  float g = RDM_Smith(LdotH, LdotN, VdotH, VdotN, roughness);

  float denom = (extIOR*LdotH + intIOR*VdotH);

  float no = intIOR;
  return color3(specularColor * (LdotH * VdotH)/(LdotN * VdotN)*((no*no)*(1.0f-f)*g*d)/(denom*denom));
}

// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(color3 diffuseColor)
{

  float pi = M_PI;
  return color3(diffuseColor / pi);
}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                texture* texture, color3 diffuseColor, color3 specularColor, float roughness, float IOR, float uTex, float vTex, int face)
{

  //! \todo compute bsdf diffuse and specular term
  if(texture != nullptr){
    color3 texColor;
    if(face == -1)
      texColor = (texture->value(uTex, vTex));
    else
      texColor = (texture->value(uTex, vTex, face));
    return color3((texColor / float(M_PI)) + RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, roughness, IOR, specularColor));
  }
  return color3(RDM_bsdf_d(diffuseColor) + RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, roughness, IOR, specularColor));
}

color3 RDM_brdf(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, float extIOR, float intIOR, float roughness, color3 specularColor)
{
  if(VdotN == 0.f)
    return color3(0,0,0);

  float d = RDM_Beckmann(NdotH, roughness);
  float f = RDM_Fresnel(LdotH, extIOR, intIOR);
  //float f = schlick(LdotH, extIOR, intIOR);
  float g = RDM_Smith(LdotH, LdotN, VdotH, VdotN, roughness);

  return specularColor * ((d * f * g) / (4.f * LdotN * VdotN));
}

#include "sampling/sampling.h"

inline float ErfInv(float x) {
    float w, p;
    x = clamp(x, -.99999f, .99999f);
    w = -std::log((1 - x) * (1 + x));
    if (w < 5) {
        w = w - 2.5f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p * w;
        p = -3.5233877e-06f + p * w;
        p = -4.39150654e-06f + p * w;
        p = 0.00021858087f + p * w;
        p = -0.00125372503f + p * w;
        p = -0.00417768164f + p * w;
        p = 0.246640727f + p * w;
        p = 1.50140941f + p * w;
    } else {
        w = std::sqrt(w) - 3;
        p = -0.000200214257f;
        p = 0.000100950558f + p * w;
        p = 0.00134934322f + p * w;
        p = -0.00367342844f + p * w;
        p = 0.00573950773f + p * w;
        p = -0.0076224613f + p * w;
        p = 0.00943887047f + p * w;
        p = 1.00167406f + p * w;
        p = 2.83297682f + p * w;
    }
    return p * x;
}

inline float Erf(float x) {
    // constants
    float a1 = 0.254829592f;
    float a2 = -0.284496736f;
    float a3 = 1.421413741f;
    float a4 = -1.453152027f;
    float a5 = 1.061405429f;
    float p = 0.3275911f;

    // Save the sign of x
    int sign = 1;
    if (x < 0) sign = -1;
    x = std::abs(x);

    // A&S formula 7.1.26
    float t = 1 / (1 + p * x);
    float y =
        1 -
        (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

    return sign * y;
}

// Microfacet Utility Functions
static void BeckmannSample11(float cosThetaI, float U1, float U2,
                             float *slope_x, float *slope_y) {
    /* Special case (normal incidence) */
    if (cosThetaI > .9999) {
        float r = std::sqrt(-std::log(1.0f - U1));
        float sinPhi = std::sin(2 * Pi * U2);
        float cosPhi = std::cos(2 * Pi * U2);
        *slope_x = r * cosPhi;
        *slope_y = r * sinPhi;
        return;
    }

    /* The original inversion routine from the paper contained
       discontinuities, which causes issues for QMC integration
       and techniques like Kelemen-style MLT. The following code
       performs a numerical inversion with better behavior */
    float sinThetaI =
        std::sqrt(std::max((float)0, (float)1 - cosThetaI * cosThetaI));
    float tanThetaI = sinThetaI / cosThetaI;
    float cotThetaI = 1 / tanThetaI;

    /* Search interval -- everything is parameterized
       in the Erf() domain */
    float a = -1, c = Erf(cotThetaI);
    float sample_x = std::max(U1, (float)1e-6f);

    /* Start with a good initial guess */
    // float b = (1-sample_x) * a + sample_x * c;

    /* We can do better (inverse of an approximation computed in
     * Mathematica) */
    float thetaI = std::acos(cosThetaI);
    float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
    float b = c - (1 + c) * std::pow(1 - sample_x, fit);

    /* Normalization factor for the CDF */
    static const float SQRT_PI_INV = 1.f / std::sqrt(Pi);
    float normalization =
        1 /
        (1 + c + SQRT_PI_INV * tanThetaI * std::exp(-cotThetaI * cotThetaI));

    int it = 0;
    while (++it < 10) {
        /* Bisection criterion -- the oddly-looking
           Boolean expression are intentional to check
           for NaNs at little additional cost */
        if (!(b >= a && b <= c)) b = 0.5f * (a + c);

        /* Evaluate the CDF and its derivative
           (i.e. the density function) */
        float invErf = ErfInv(b);
        float value =
            normalization *
                (1 + b + SQRT_PI_INV * tanThetaI * std::exp(-invErf * invErf)) -
            sample_x;
        float derivative = normalization * (1 - invErf * tanThetaI);

        if (std::abs(value) < 1e-5f) break;

        /* Update bisection intervals */
        if (value > 0)
            c = b;
        else
            a = b;

        b -= value / derivative;
    }

    /* Now convert back into a slope value */
    *slope_x = ErfInv(b);

    /* Simulate Y component */
    *slope_y = ErfInv(2.0f * std::max(U2, (float)1e-6f) - 1.0f);

}

inline float CosPhi(const vec3& w){
  auto phi = std::atan2(w.y, w.x);
  return cos(phi);
}

inline float SinPhi(const vec3& w){
  auto phi = std::atan2(w.y, w.x);
  return sin(phi);
}

static vec3 BeckmannSample(const vec3 &wi, float alpha,
                               float U1, float U2, const vec3& normal) {
    // 1. stretch wi
    vec3 wiStretched =
        normalize(vec3(alpha * wi.x, alpha * wi.y, wi.z));

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    float slope_x, slope_y;
    BeckmannSample11(dot(wiStretched, normal), U1, U2, &slope_x, &slope_y);

    // 3. rotate
    float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
    slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha * slope_x;
    slope_y = alpha * slope_y;

    // 5. compute normal
    return normalize(vec3(-slope_x, -slope_y, 1.f));
}

float BeckmannPdf(const vec3& wo, const vec3& wh, const vec3& n, float roughness, float LdotH, float LdotN) {
  float NdotH = dot(n, wh);
  float VdotN = dot(wo, n);
  float VdotH = dot(wo, n);
  return RDM_Beckmann(NdotH, roughness) * RDM_Smith(LdotH, LdotN, VdotH, VdotN, roughness) * abs(dot(wo, wh)) / abs(VdotN);
}

vec3 Beckmann_Sample_wh(const vec3& wo, const point2& u, float roughness, const vec3& normal) {
  vec3 wh = BeckmannSample(wo, roughness, u[0], u[1], normal);
  wh = -wh;
  return wh;
}
