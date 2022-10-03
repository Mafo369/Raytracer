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
                  float VdotN, const Material *m)
{

  //! \todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G
  //! = RDM_Smith

  float d = RDM_Beckmann(NdotH, m->roughness);
  float f = RDM_Fresnel(LdotH, 1.f, m->IOR);
  //float f = schlick(LdotH, m->IOR, 1.f);
  float g = RDM_Smith(LdotH, LdotN, VdotH, VdotN, m->roughness);

  return color3(m->specularColor * ((d * f * g) / (4 * LdotN * VdotN)));
}

color3 RDM_btdf(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, const Material *m, float extIOR, float intIOR)
{

  float d = RDM_Beckmann(NdotH, m->roughness);
  float f = RDM_Fresnel(LdotH, extIOR, intIOR);
  //float f = schlick(LdotH, extIOR, intIOR);
  float g = RDM_Smith(LdotH, LdotN, VdotH, VdotN, m->roughness);

  float denom = (extIOR*LdotH + intIOR*VdotH);

  float no = intIOR;
  return color3(m->specularColor * (LdotH * VdotH)/(LdotN * VdotN)*((no*no)*(1.0f-f)*g*d)/(denom*denom));
}

// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(const Material *m)
{

  float pi = M_PI;
  return color3(m->diffuseColor / pi);
}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                const Material *m, float uTex, float vTex, int face)
{

  //! \todo compute bsdf diffuse and specular term
  if(m->m_texture != nullptr){
    color3 texColor;
    if(face == -1)
      texColor = (m->m_texture->value(uTex, vTex));
    else
      texColor = (m->m_texture->value(uTex, vTex, face));
    return color3((texColor / float(M_PI)) + RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, m));
  }
  return color3(RDM_bsdf_d(m) + RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, m));
}

color3 RDM_brdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                const Material *m, float extIOR, float intIOR )
{
  if(VdotN == 0.f)
    return color3(0,0,0);

  float d = RDM_Beckmann(NdotH, m->roughness);
  float f = RDM_Fresnel(LdotH, extIOR, intIOR);
  //float f = schlick(LdotH, extIOR, intIOR);
  float g = RDM_Smith(LdotH, LdotN, VdotH, VdotN, m->roughness);

  return m->specularColor * ((d * f * g) / (4.f * LdotN * VdotN));
}
