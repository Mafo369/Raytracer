#pragma once

#include "scene.h"

color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                const Material *m, float uTex, float vTex, int face);

color3 RDM_brdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN,
                const Material *m, float extIOR, float intIOR );

color3 RDM_btdf(float LdotH, float NdotH, float VdotH, float LdotN,
                  float VdotN, const Material *m, float extIOR, float intIOR);

float RDM_Fresnel(float LdotH, float extIOR, float intIOR);

float schlick(float VdotN, float extIOR, float intIOR);
