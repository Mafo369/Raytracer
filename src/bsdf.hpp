#pragma once

#include "Material.h"
#include "scene.h"

color3 RDM_bsdf( BrdfData& data, texture* texture, int face );

color3 RDM_bsdf_s( BrdfData& data );

color3 RDM_brdf( float LdotH,
                 float NdotH,
                 float VdotH,
                 float LdotN,
                 float VdotN,
                 float extIOR,
                 float intIOR,
                 float roughness,
                 color3 specularColor );

color3 RDM_btdf( float LdotH,
                 float NdotH,
                 float VdotH,
                 float LdotN,
                 float VdotN,
                 float extIOR,
                 float intIOR,
                 float roughness,
                 color3 specularColor );

float RDM_Fresnel( float LdotH, float extIOR, float intIOR );

float schlick( float VdotN, float extIOR, float intIOR );

float BeckmannPdf( const vec3& wo,
                   const vec3& wh,
                   const vec3& n,
                   float roughness,
                   float LdotH,
                   float LdotN );

vec3 Beckmann_Sample_wh( const vec3& wo, const point2& u, float roughness, const vec3& normal );
