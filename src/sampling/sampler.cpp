
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

// core/sampler.cpp*
#include "sampler.h"
#include "sampling.h"
#include <iostream>

// Sampler Method Definitions
Sampler::~Sampler() {}

Sampler::Sampler( int64_t samplesPerPixel ) : samplesPerPixel( samplesPerPixel ) {}

CameraSample Sampler::GetCameraSample( const point2& pRaster ) {
    CameraSample sample;
    sample.xy = (point2)pRaster + Get2D();
    sample.t  = Get1D();
    sample.uv = Get2D();
    return sample;
}

void Sampler::StartPixel( const point2& p ) {
    currentPixel            = p;
    currentPixelSampleIndex = 0;
    // Reset array offsets for next pixel sample
    array1DOffset = array2DOffset = 0;
}

bool Sampler::StartNextSample() {
    // Reset array offsets for next pixel sample
    array1DOffset = array2DOffset = 0;
    return ++currentPixelSampleIndex < samplesPerPixel;
}

bool Sampler::SetSampleNumber( int64_t sampleNum ) {
    // Reset array offsets for next pixel sample
    array1DOffset = array2DOffset = 0;
    currentPixelSampleIndex       = sampleNum;
    return currentPixelSampleIndex < samplesPerPixel;
}

void Sampler::Request1DArray( int n ) {
    samples1DArraySizes.push_back( n );
    sampleArray1D.push_back( std::vector<float>( n * samplesPerPixel ) );
}

void Sampler::Request2DArray( int n ) {
    samples2DArraySizes.push_back( n );
    sampleArray2D.push_back( std::vector<point2>( n * samplesPerPixel ) );
}

const float* Sampler::Get1DArray( int n ) {
    if ( array1DOffset == sampleArray1D.size() ) return nullptr;
    return &sampleArray1D[array1DOffset++][currentPixelSampleIndex * n];
}

const point2* Sampler::Get2DArray( int n ) {
    if ( array2DOffset == sampleArray2D.size() ) return nullptr;
    return &sampleArray2D[array2DOffset++][currentPixelSampleIndex * n];
}

PixelSampler::PixelSampler( int64_t samplesPerPixel, int nSampledDimensions ) :
    Sampler( samplesPerPixel ) {
    for ( int i = 0; i < nSampledDimensions; ++i ) {
        samples1D.push_back( std::vector<float>( samplesPerPixel ) );
        samples2D.push_back( std::vector<point2>( samplesPerPixel ) );
    }
}

bool PixelSampler::StartNextSample() {
    current1DDimension = current2DDimension = 0;
    return Sampler::StartNextSample();
}

bool PixelSampler::SetSampleNumber( int64_t sampleNum ) {
    current1DDimension = current2DDimension = 0;
    return Sampler::SetSampleNumber( sampleNum );
}

float PixelSampler::Get1D() {
    if ( current1DDimension < samples1D.size() )
        return samples1D[current1DDimension++][currentPixelSampleIndex];
    else
        return rng.Uniformfloat();
}

point2 PixelSampler::Get2D() {
    if ( current2DDimension < samples2D.size() )
        return samples2D[current2DDimension++][currentPixelSampleIndex];
    else
        return point2( rng.Uniformfloat(), rng.Uniformfloat() );
}

void GlobalSampler::StartPixel( const point2& p ) {
    Sampler::StartPixel( p );
    dimension           = 0;
    intervalSampleIndex = GetIndexForSample( 0 );
    // Compute _arrayEndDim_ for dimensions used for array samples
    arrayEndDim = arrayStartDim + sampleArray1D.size() + 2 * sampleArray2D.size();

    // Compute 1D array samples for _GlobalSampler_
    for ( size_t i = 0; i < samples1DArraySizes.size(); ++i ) {
        int nSamples = samples1DArraySizes[i] * samplesPerPixel;
        for ( int j = 0; j < nSamples; ++j ) {
            int64_t index       = GetIndexForSample( j );
            sampleArray1D[i][j] = SampleDimension( index, arrayStartDim + i );
        }
    }

    // Compute 2D array samples for _GlobalSampler_
    int dim = arrayStartDim + samples1DArraySizes.size();
    for ( size_t i = 0; i < samples2DArraySizes.size(); ++i ) {
        int nSamples = samples2DArraySizes[i] * samplesPerPixel;
        for ( int j = 0; j < nSamples; ++j ) {
            int64_t idx           = GetIndexForSample( j );
            sampleArray2D[i][j].x = SampleDimension( idx, dim );
            sampleArray2D[i][j].y = SampleDimension( idx, dim + 1 );
        }
        dim += 2;
    }
}

bool GlobalSampler::StartNextSample() {
    dimension           = 0;
    intervalSampleIndex = GetIndexForSample( currentPixelSampleIndex + 1 );
    return Sampler::StartNextSample();
}

bool GlobalSampler::SetSampleNumber( int64_t sampleNum ) {
    dimension           = 0;
    intervalSampleIndex = GetIndexForSample( sampleNum );
    return Sampler::SetSampleNumber( sampleNum );
}

float GlobalSampler::Get1D() {
    if ( dimension >= arrayStartDim && dimension < arrayEndDim ) dimension = arrayEndDim;
    return SampleDimension( intervalSampleIndex, dimension++ );
}

point2 GlobalSampler::Get2D() {
    if ( dimension + 1 >= arrayStartDim && dimension < arrayEndDim ) dimension = arrayEndDim;
    point2 p( SampleDimension( intervalSampleIndex, dimension ),
              SampleDimension( intervalSampleIndex, dimension + 1 ) );
    dimension += 2;
    return p;
}
