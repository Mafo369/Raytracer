
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

#define NOMINMAX
#pragma once

#include "rng.h"
#include "../defines.h"
#include <vector>
#include <memory>

// Sampler Declarations
class Sampler {
  public:
    // Sampler Interface
    virtual ~Sampler();
    Sampler(int64_t samplesPerPixel);
    virtual void StartPixel(const point2 &p);
    virtual float Get1D() = 0;
    virtual point2 Get2D() = 0;
    point2 GetCameraSample(const point2 &pRaster);
    void Request1DArray(int n);
    void Request2DArray(int n);
    virtual int RoundCount(int n) const { return n; }
    const float *Get1DArray(int n);
    const point2 *Get2DArray(int n);
    virtual bool StartNextSample();
    virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
    virtual bool SetSampleNumber(int64_t sampleNum);
    int64_t CurrentSampleNumber() const { return currentPixelSampleIndex; }

    // Sampler Public Data
    const int64_t samplesPerPixel;

  protected:
    // Sampler Protected Data
    point2 currentPixel;
    int64_t currentPixelSampleIndex;
    std::vector<int> samples1DArraySizes, samples2DArraySizes;
    std::vector<std::vector<float>> sampleArray1D;
    std::vector<std::vector<point2>> sampleArray2D;

  private:
    // Sampler Private Data
    size_t array1DOffset, array2DOffset;
};

class PixelSampler : public Sampler {
  public:
    // PixelSampler Public Methods
    PixelSampler(int64_t samplesPerPixel, int nSampledDimensions);
    bool StartNextSample();
    bool SetSampleNumber(int64_t);
    float Get1D();
    point2 Get2D();

  protected:
    // PixelSampler Protected Data
    std::vector<std::vector<float>> samples1D;
    std::vector<std::vector<point2>> samples2D;
    int current1DDimension = 0, current2DDimension = 0;
    RNG rng;
};

class GlobalSampler : public Sampler {
  public:
    // GlobalSampler Public Methods
    bool StartNextSample();
    void StartPixel(const point2 &);
    bool SetSampleNumber(int64_t sampleNum);
    float Get1D();
    point2 Get2D();
    GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) {}
    virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
    virtual float SampleDimension(int64_t index, int dimension) const = 0;

  private:
    // GlobalSampler Private Data
    int dimension;
    int64_t intervalSampleIndex;
    static const int arrayStartDim = 5;
    int arrayEndDim;
};
