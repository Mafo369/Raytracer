
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#pragma once

#include "rng.h"
#include "../defines.h"
#include <algorithm>
#include <memory>

const float Pi = M_PI;
const float InvPi = 1.f / Pi;
const float Inv2Pi = 1.f / (2.f * Pi);
const float Inv4Pi = 1.f / (4.f * Pi);
const float PiOver2 = Pi / 2.f;
const float PiOver4 = Pi / 4.f;

// Sampling Declarations
void StratifiedSample1D(float *samples, int nsamples, RNG &rng,
                        bool jitter = true);
void StratifiedSample2D(point2 *samples, int nx, int ny, RNG &rng,
                        bool jitter = true);
void LatinHypercube(float *samples, int nSamples, int nDim, RNG &rng);

template <typename Predicate>
int FindInterval(int size, const Predicate &pred) {
    int first = 0, len = size;
    while (len > 0) {
        int half = len >> 1, middle = first + half;
        // Bisect range based on value of _pred_ at _middle_
        if (pred(middle)) {
            first = middle + 1;
            len -= half + 1;
        } else
            len = half;
    }
    return clamp(first - 1, 0, size - 2);
}

struct Distribution1D {
    // Distribution1D Public Methods
    Distribution1D(const float *f, int n) : func(f, f + n), cdf(n + 1) {
        // Compute integral of step function at $x_i$
        cdf[0] = 0;
        for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

        // Transform step function integral into CDF
        funcInt = cdf[n];
        if (funcInt == 0) {
            for (int i = 1; i < n + 1; ++i) cdf[i] = float(i) / float(n);
        } else {
            for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
        }
    }
    int Count() const { return (int)func.size(); }
    float SampleContinuous(float u, float *pdf, int *off = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval((int)cdf.size(),
                                  [&](int index) { return cdf[index] <= u; });
        if (off) *off = offset;
        // Compute offset along CDF segment
        float du = u - cdf[offset];
        if ((cdf[offset + 1] - cdf[offset]) > 0) {
            du /= (cdf[offset + 1] - cdf[offset]);
        }

        // Compute PDF for sampled offset
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;

        // Return $x\in{}[0,1)$ corresponding to sample
        return (offset + du) / Count();
    }
    int SampleDiscrete(float u, float *pdf = nullptr,
                       float *uRemapped = nullptr) const {
        // Find surrounding CDF segments and _offset_
        int offset = FindInterval((int)cdf.size(),
                                  [&](int index) { return cdf[index] <= u; });
        if (pdf) *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
        if (uRemapped)
            *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
        return offset;
    }
    float DiscretePDF(int index) const {
        return func[index] / (funcInt * Count());
    }

    // Distribution1D Public Data
    std::vector<float> func, cdf;
    float funcInt;
};

point2 RejectionSampleDisk(RNG &rng);
vec3 UniformSampleHemisphere(const point2 &u);
float UniformHemispherePdf();
vec3 UniformSampleSphere(const point2 &u);
float UniformSpherePdf();
vec3 UniformSampleCone(const point2 &u, float thetamax);
vec3 UniformSampleCone(const point2 &u, float thetamax, const vec3 &x,
                           const vec3 &y, const vec3 &z);
float UniformConePdf(float thetamax);
point2 UniformSampleDisk(const point2 &u);
point2 ConcentricSampleDisk(const point2 &u);
point2 UniformSampleTriangle(const point2 &u);
class Distribution2D {
  public:
    // Distribution2D Public Methods
    Distribution2D(const float *data, int nu, int nv);
    point2 SampleContinuous(const point2 &u, float *pdf) const {
        float pdfs[2];
        int v;
        float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
        float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
        return point2(d0, d1);
    }
    float Pdf(const point2 &p) const {
        int iu = clamp(int(p[0] * pConditionalV[0]->Count()), 0,
                       pConditionalV[0]->Count() - 1);
        int iv =
            clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
        return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
    }

  private:
    // Distribution2D Private Data
    std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
    std::unique_ptr<Distribution1D> pMarginal;
};

// Sampling Inline Functions
template <typename T>
void Shuffle(T *samp, int count, int nDimensions, RNG &rng) {
    for (int i = 0; i < count; ++i) {
        int other = i + rng.UniformUInt32(count - i);
        for (int j = 0; j < nDimensions; ++j)
            std::swap(samp[nDimensions * i + j], samp[nDimensions * other + j]);
    }
}

inline vec3 CosineSampleHemisphere(const point2 &u) {
    point2 d = ConcentricSampleDisk(u);
    float z = std::sqrt(std::max((float)0, 1 - d.x * d.x - d.y * d.y));
    return vec3(d.x, d.y, z);
}

inline float CosineHemispherePdf(float cosTheta) { return cosTheta * InvPi; }

inline float BalanceHeuristic(int nf, float fPdf, int ng, float gPdf) {
    return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}

inline float PowerHeuristic(int nf, float fPdf, int ng, float gPdf) {
    float f = nf * fPdf, g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

//inline vec3 cosineWeightedSampling(vec3 normal) {
//    float phi = uniform01(engine) * 2.0 * M_PI;
//    float the = acos(1.0 - 2.0 * uniform01(engine)) / 2.0;
//
//    vec3 v0 = vec3(0, 1.0, 0);
//    if(dot(v0, normal) > 0.5 || dot(v0, normal) < -0.5)
//      v0 = vec3(0,0,1.0);
//    vec3 v1 = normalize(cross(v0, normal));
//    v0 = normalize(cross(v1, normal));
//
//    return normal * cos(the) + (v0 * cos(phi) + v1 * sin(phi)) * sin(the);
//}

inline vec3 random_dir(const vec3& normal) {
  float r1 = uniform01(engine);
  float r2 = uniform01(engine);

  float theta = 2.f * Pi * r1;
  float r = sqrt(1.f-r2);
  vec3 dir_rand_local = vec3(cos(theta)*r, sin(theta)*r, sqrt(r2));
  vec3 randV = vec3(uniform01(engine)-0.5, uniform01(engine)-0.5, uniform01(engine)-0.5);

  vec3 tangent1 = normalize(cross(normal, randV));
  vec3 tangent2 = cross(tangent1, normal);

  return dir_rand_local[2] * normal + dir_rand_local[0] * tangent1 + dir_rand_local[1] * tangent2;
  //float r1 = uniform01(engine);
	//float r2 = uniform01(engine);
	//float sr2 = sqrt(1. - r2);
	//vec3 direction_aleatoire_repere_local(cos(2 * M_PI*r1)*sr2, sin(2 * M_PI*r1)*sr2, sqrt(r2));
	///*vec3 aleatoire(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);	
	//vec3 tangent1 = cross(N, aleatoire); tangent1.normalize();*/
	//vec3 tangent1;
	//vec3 absN(abs(N[0]), abs(N[1]), abs(N[2]));
	//if (absN[0] <= absN[1] && absN[0]<=absN[2]) {
	//	tangent1 = vec3(0, -N[2], N[1]);
	//}else
	//	if (absN[1] <= absN[0] && absN[1]<=absN[2]) {
	//		tangent1 = vec3(-N[2], 0, N[0]);
	//	} else
	//		tangent1 = vec3(-N[1], N[0], 0);
  //tangent1 = normalize(tangent1);
	//vec3 tangent2 = cross(tangent1, N);

	//return direction_aleatoire_repere_local[2] * N + direction_aleatoire_repere_local[0] * tangent1 + direction_aleatoire_repere_local[1] * tangent2;
}

inline vec3 random_uniform() {
  float r1 = uniform01(engine);
  float r2 = uniform01(engine);
  vec3 result;
  result[0] = 2.f * cos(2.f * Pi * r1)*sqrt(r2*(1. - r2));
  result[1] = 2.f * sin(2.f * Pi * r1)*sqrt(r2*(1. - r2));
  result[2] = 1.f - 2.f * r2;
  return result;
}

inline vec3 uniformHemisphereDir(vec3 normal){
  float r1 = uniform01(engine);
  float r2 = uniform01(engine);
  float theta = 2.f * Pi * r1;
  float r = sqrt(1.f-(r2*r2));

  vec3 dir_rand_local = vec3(cos(theta)*r, sin(theta)*r, r2);
  vec3 randV = vec3(uniform01(engine)-0.5, uniform01(engine)-0.5, uniform01(engine)-0.5);

  vec3 tangent1 = normalize(cross(normal, randV));
  vec3 tangent2 = cross(tangent1, normal);

  return dir_rand_local[2] * normal + dir_rand_local[0] * tangent1 + dir_rand_local[1] * tangent2;
}
