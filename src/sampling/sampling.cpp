
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

#include "sampling.h"
#include <cmath>

// Sampling Function Definitions
void StratifiedSample1D( float* samp, int nSamples, RNG& rng, bool jitter ) {
    float invNSamples = (float)1 / nSamples;
    for ( int i = 0; i < nSamples; ++i ) {
        float delta = jitter ? rng.Uniformfloat() : 0.5f;
        samp[i]     = std::min( ( i + delta ) * invNSamples, OneMinusEpsilon );
    }
}

void StratifiedSample2D( point2* samp, int nx, int ny, RNG& rng, bool jitter ) {
    float dx = (float)1 / nx, dy = (float)1 / ny;
    for ( int y = 0; y < ny; ++y )
        for ( int x = 0; x < nx; ++x ) {
            float jx = jitter ? rng.Uniformfloat() : 0.5f;
            float jy = jitter ? rng.Uniformfloat() : 0.5f;
            samp->x  = std::min( ( x + jx ) * dx, OneMinusEpsilon );
            samp->y  = std::min( ( y + jy ) * dy, OneMinusEpsilon );
            ++samp;
        }
}

void LatinHypercube( float* samples, int nSamples, int nDim, RNG& rng ) {
    // Generate LHS samples along diagonal
    float invNSamples = (float)1 / nSamples;
    for ( int i = 0; i < nSamples; ++i )
        for ( int j = 0; j < nDim; ++j ) {
            float sj              = ( i + ( rng.Uniformfloat() ) ) * invNSamples;
            samples[nDim * i + j] = std::min( sj, OneMinusEpsilon );
        }

    // Permute LHS samples in each dimension
    for ( int i = 0; i < nDim; ++i ) {
        for ( int j = 0; j < nSamples; ++j ) {
            int other = j + rng.UniformUInt32( nSamples - j );
            std::swap( samples[nDim * j + i], samples[nDim * other + i] );
        }
    }
}

point2 RejectionSampleDisk( RNG& rng ) {
    point2 p;
    do {
        p.x = 1 - 2 * rng.Uniformfloat();
        p.y = 1 - 2 * rng.Uniformfloat();
    } while ( p.x * p.x + p.y * p.y > 1 );
    return p;
}

vec3 UniformSampleHemisphere( const point2& u ) {
    float z   = u[0];
    float r   = std::sqrt( std::max( (float)0, (float)1. - z * z ) );
    float phi = 2 * Pi * u[1];
    return vec3( r * std::cos( phi ), r * std::sin( phi ), z );
}

float UniformHemispherePdf() {
    return Inv2Pi;
}

vec3 UniformSampleSphere( const point2& u ) {
    float z   = 1 - 2 * u[0];
    float r   = std::sqrt( std::max( (float)0, (float)1 - z * z ) );
    float phi = 2 * Pi * u[1];
    return vec3( r * std::cos( phi ), r * std::sin( phi ), z );
}

float UniformSpherePdf() {
    return Inv4Pi;
}

point2 UniformSampleDisk( const point2& u ) {
    float r     = std::sqrt( u[0] );
    float theta = 2 * Pi * u[1];
    return point2( r * std::cos( theta ), r * std::sin( theta ) );
}

point2 ConcentricSampleDisk( const point2& u ) {
    // Map uniform random numbers to $[-1,1]^2$
    point2 uOffset = 2.f * u - vec2( 1, 1 );

    // Handle degeneracy at the origin
    if ( uOffset.x == 0 && uOffset.y == 0 ) return point2( 0, 0 );

    // Apply concentric mapping to point
    float theta, r;
    if ( std::abs( uOffset.x ) > std::abs( uOffset.y ) ) {
        r     = uOffset.x;
        theta = PiOver4 * ( uOffset.y / uOffset.x );
    }
    else {
        r     = uOffset.y;
        theta = PiOver2 - PiOver4 * ( uOffset.x / uOffset.y );
    }
    return r * point2( std::cos( theta ), std::sin( theta ) );
}

float UniformConePdf( float cosThetaMax ) {
    return 1 / ( 2 * Pi * ( 1 - cosThetaMax ) );
}

vec3 UniformSampleCone( const point2& u, float cosThetaMax ) {
    float cosTheta = ( (float)1 - u[0] ) + u[0] * cosThetaMax;
    float sinTheta = std::sqrt( (float)1 - cosTheta * cosTheta );
    float phi      = u[1] * 2 * Pi;
    return vec3( std::cos( phi ) * sinTheta, std::sin( phi ) * sinTheta, cosTheta );
}

vec3 UniformSampleCone( const point2& u,
                        float cosThetaMax,
                        const vec3& x,
                        const vec3& y,
                        const vec3& z ) {
    float cosTheta = std::lerp( u[0], cosThetaMax, 1.f );
    float sinTheta = std::sqrt( (float)1. - cosTheta * cosTheta );
    float phi      = u[1] * 2 * Pi;
    return std::cos( phi ) * sinTheta * x + std::sin( phi ) * sinTheta * y + cosTheta * z;
}

point2 UniformSampleTriangle( const point2& u ) {
    float su0 = std::sqrt( u[0] );
    return point2( 1 - su0, u[1] * su0 );
}

Distribution2D::Distribution2D( const float* func, int nu, int nv ) {
    pConditionalV.reserve( nv );
    for ( int v = 0; v < nv; ++v ) {
        // Compute conditional sampling distribution for $\tilde{v}$
        pConditionalV.emplace_back( new Distribution1D( &func[v * nu], nu ) );
    }
    // Compute marginal sampling distribution $p[\tilde{v}]$
    std::vector<float> marginalFunc;
    marginalFunc.reserve( nv );
    for ( int v = 0; v < nv; ++v )
        marginalFunc.push_back( pConditionalV[v]->funcInt );
    pMarginal.reset( new Distribution1D( &marginalFunc[0], nv ) );
}
