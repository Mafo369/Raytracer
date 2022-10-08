#pragma once

#include "scene.h"

color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                texture* texture, color3 diffuseColor, color3 specularColor, float roughness, float IOR, float uTex, float vTex, int face);

color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, float roughness, float IOR, color3 specularColor);

color3 RDM_brdf(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, float extIOR, float intIOR, float roughness, color3 specularColor);

color3 RDM_btdf(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, float extIOR, float intIOR, float roughness, color3 specularColor);

float RDM_Fresnel(float LdotH, float extIOR, float intIOR);

float schlick(float VdotN, float extIOR, float intIOR);
