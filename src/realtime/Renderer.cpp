#include "Renderer.h"

#include "Walnut/Random.h"

#include "../integrator.h"

namespace Utils {

static uint32_t ConvertToRGBA( const glm::vec4& color ) {
    uint8_t r = (uint8_t)( color.r * 255.0f );
    uint8_t g = (uint8_t)( color.g * 255.0f );
    uint8_t b = (uint8_t)( color.b * 255.0f );
    uint8_t a = (uint8_t)( color.a * 255.0f );

    uint32_t result = ( a << 24 ) | ( b << 16 ) | ( g << 8 ) | r;
    return result;
}

} // namespace Utils

Renderer::Renderer() {}

void Renderer::OnResize( uint32_t width, uint32_t height ) {
    if ( m_FinalImage ) {
        // No resize necessary
        if ( m_FinalImage->GetWidth() == width && m_FinalImage->GetHeight() == height ) return;

        m_FinalImage->Resize( width, height );
        img = initImage( m_FinalImage->GetWidth(), m_FinalImage->GetHeight() );
        setCameraFOV( scene,
                      point3( 0, 0, 0 ),
                      vec3( 0, 0, 0 ),
                      vec3( 0, 0, 0 ),
                      60,
                      (float)m_FinalImage->GetWidth(),
                      (float)m_FinalImage->GetHeight(),
                      0.01 );
    }
    else {
        m_FinalImage = std::make_shared<Walnut::Image>( width, height, Walnut::ImageFormat::RGBA );
        setCameraFOV( scene,
                      point3( 0, 0, 0 ),
                      vec3( 0, 0, 0 ),
                      vec3( 0, 0, 0 ),
                      60,
                      (float)m_FinalImage->GetWidth(),
                      (float)m_FinalImage->GetHeight(),
                      0.01 );
    }

    delete[] m_ImageData;
    m_ImageData = new uint32_t[width * height];
    delete[] m_accumulationData;
    m_accumulationData = new glm::vec4[width * height];
}

void Renderer::Render( const CameraI& camera ) {
    if ( m_frameIndex == 1 )
        memset( m_accumulationData,
                0,
                m_FinalImage->GetHeight() * m_FinalImage->GetWidth() * sizeof( glm::vec4 ) );

    for ( uint32_t y = 0; y < m_FinalImage->GetHeight(); y++ ) {
#pragma omp parallel
        {
            std::unique_ptr<Sampler> tileSampler = sampler->Clone( time( NULL ) );
#pragma omp for schedule( dynamic )
            for ( uint32_t x = 0; x < m_FinalImage->GetWidth(); x++ ) {
                auto pixel = vec2( x, y );

                vec4 color( 0 );
                tileSampler->StartPixel( pixel );
                do {
                    vec3 dir = camera.GetRayDirections()[x + y * m_FinalImage->GetWidth()];
                    Ray ray  = Ray( camera.GetPosition(), normalize( dir ) );
                    ray.hasDifferentials = true;

                    Intersection intersection;
                    auto colorr = m_integrator->trace_ray(
                        scene, ray, tree, &intersection, tileSampler.get() );
                    color += vec4( colorr.r, colorr.g, colorr.b, 1 );
                } while ( tileSampler->StartNextSample() );

                color = color / (float)tileSampler->samplesPerPixel;

                m_accumulationData[x + y * m_FinalImage->GetWidth()] += color;

                glm::vec4 accumulatedColor = m_accumulationData[x + y * m_FinalImage->GetWidth()];
                accumulatedColor /= (float)m_frameIndex;

                // gamma-correction
                accumulatedColor.r = powf( accumulatedColor.r, 1.0f / 2.2 );
                accumulatedColor.g = powf( accumulatedColor.g, 1.0f / 2.2 );
                accumulatedColor.b = powf( accumulatedColor.b, 1.0f / 2.2 );

                accumulatedColor =
                    glm::clamp( accumulatedColor, glm::vec4( 0.0f ), glm::vec4( 1.0f ) );
                m_ImageData[x + y * m_FinalImage->GetWidth()] =
                    Utils::ConvertToRGBA( accumulatedColor );
            }
        }
    }

    m_FinalImage->SetData( m_ImageData );

    if ( m_settings.accumulate )
        m_frameIndex++;
    else
        m_frameIndex = 1;
}
