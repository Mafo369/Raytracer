#include "CPURenderer.h"
#include "../sampling/stratified.h"
#include <chrono>
#include "../Camera.h"
#include "../image.h"

void CPURenderer::init(std::string& imageName, glm::vec2 resolution, Scene* scene){
    m_imageName = imageName;
    m_scene = scene;
    m_renderImage = initImage( resolution.x, resolution.y );

    std::cout << "Output resolution: " << resolution.x << "x" << resolution.y << std::endl;

    m_kdTree = initKdTree( m_scene );
    m_sampler = new StratifiedSampler( 8, 8, true, 8 );
    std::cout << "Spp: " << m_sampler->samplesPerPixel << std::endl;

    m_tracer = new Pathtracer(5, m_sampler);
    m_tracer->preprocess(scene, m_sampler);
}

void CPURenderer::destroy(){
    freeImage(m_renderImage);
    m_renderImage = NULL;
}

void CPURenderer::run(){
    auto startTime = std::chrono::system_clock::now();

    for ( size_t j = 0; j < m_renderImage->height; j++ ) {
        if ( j != 0 ) printf( "\033[A\r" );
        float progress = (float)j / m_renderImage->height * 100.f;
        printf( "progress\t[" );
        int cpt = 0;
        for ( cpt = 0; cpt < progress; cpt += 5 )
            printf( "." );
        for ( ; cpt < 100; cpt += 5 )
            printf( " " );
        printf( "]\n" );
#pragma omp parallel
        {
            std::unique_ptr<Sampler> tileSampler = m_sampler->Clone( time( NULL ) );
#pragma omp for schedule( dynamic )
            for ( size_t i = 0; i < m_renderImage->width; i++ ) {
                color3 pixel_color( 0, 0, 0 );
                color3* ptr = getPixelPtr( m_renderImage, i, j );
                auto pixel  = vec2( i, j );

                tileSampler->StartPixel( pixel );
                do {
                    CameraSample cameraSample = tileSampler->GetCameraSample( pixel );
                    Ray rx = m_scene->cam->get_ray( cameraSample.xy.x,
                                         cameraSample.xy.y,
                                         cameraSample.uv.x,
                                         cameraSample.uv.y,
                                         vec2( int( i ), int( j ) ) );
                    if ( rx.hasDifferentials )
                        scaleDifferentials( &rx, 1.f / sqrt( tileSampler->samplesPerPixel ) );
                    Intersection intersection;
                    pixel_color += m_tracer->trace_ray( m_scene, &rx, m_kdTree, &intersection, tileSampler.get() );
                } while ( tileSampler->StartNextSample() );

                color3 avgColor = pixel_color / (float)tileSampler->samplesPerPixel;

                // gamma-correction
                avgColor.r = powf( avgColor.r, 1.0f / 2.2 );
                avgColor.g = powf( avgColor.g, 1.0f / 2.2 );
                avgColor.b = powf( avgColor.b, 1.0f / 2.2 );

                // if(avgColor.x < 1 || avgColor.y < 1 || avgColor.z < 1)
                //  std::cout << glm::to_string(avgColor) << std::endl;
                *ptr = avgColor;
            }
        }
    }
    auto stopTime = std::chrono::system_clock::now();
    std::cout
        << "Rendering took "
        << std::chrono::duration_cast<std::chrono::duration<double>>( stopTime - startTime ).count()
        << "s" << std::endl;

    saveImage( m_renderImage, m_imageName.data() );
}

void CPURenderer::scaleDifferentials( Ray* ray, float s ) {
    ray->dox = ray->orig + ( ray->dox - ray->orig ) * s;
    ray->doy = ray->orig + ( ray->doy - ray->orig ) * s;
    ray->ddx = ray->dir + ( ray->ddx - ray->dir ) * s;
    ray->ddy = ray->dir + ( ray->ddy - ray->dir ) * s;
}
