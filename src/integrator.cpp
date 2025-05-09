#include "integrator.h"
#include "LightDistribution.h"
#include "Medium.h"

// pathtracer
color3 trace_ray( Scene* scene,
                  Ray* ray,
                  KdTree* tree,
                  Intersection* intersection,
                  bool show_lights,
                  Sampler* sampler ) {
    color3 ret = color3( 0, 0, 0 );
    color3 throughput( 1.f );

    // if ( ray->depth > scene->depth ) return color3( 0.f );

    int depth;

    for ( depth = 0;; depth++ ) {

        if ( intersectKdTree( scene, tree, ray, intersection ) ) {
            intersection->hit = true;
            // Compute necessary differential information for texture filtering
            if ( ray->hasDifferentials ) intersection->computeDifferentials( ray );

            const Material& intersectedMaterial = scene->GetMaterial( intersection->materialIndex );

            if ( !isBlack( intersectedMaterial.m_emission ) && depth == 0 ) {
                ret += throughput *
                       static_cast<ShapeLight*>( scene->lights[0] )->L( *intersection, -ray->dir );
            }

            for ( auto& light : scene->lights ) {
                vec2 uLight      = sampler->Get2D();
                vec2 uScattering = sampler->Get2D();
                ret +=
                    throughput * directIllumination(
                                     scene, tree, ray, intersection, light, uScattering, uLight );
            }
            float pdf;
            vec3 wi;
            color3 bsdf =
                intersectedMaterial.sample( ray, intersection, sampler->Get2D(), &wi, &pdf );
            wi = normalize( wi );
            if ( isBlack( bsdf ) || pdf == 0.f ) break;

            throughput *= bsdf * abs( dot( wi, intersection->normal ) ) / pdf;

            ray = new Ray( intersection->position + acne_eps * wi, wi, 0, 10000, 0 );
        }
        else {
            for ( auto& env : scene->envLights )
                ret += throughput * env->Le( ray );
            break;
        }

        if ( depth >= 5 ) break;

        // if ( ray->tmax < 0 || !intersection->hit || ray->dir.z == 0 ) return ret;

        if ( scene->medium != nullptr )
            ret = ret * scene->medium->tr( *ray, scene->ysol ) +
                  scene->medium->sample( *ray, scene, tree, scene->ysol );
    }

    return ret;
}

color3 directIllumination( Scene* scene,
                           KdTree* tree,
                           Ray* ray,
                           Intersection* intersection,
                           Light* light,
                           const vec2& uScattering,
                           const vec2& uLight ) {
    color3 ret( 0 );

    vec3 wi;
    float lightPdf = 0, scatteringPdf = 0;
    bool vis;
    auto Li = light->sample_Li( scene, tree, *intersection, uLight, &wi, &lightPdf, &vis );
    const Material& intersectedMaterial = scene->GetMaterial( intersection->materialIndex );
    if ( lightPdf > 0 && !isBlack( Li ) ) {
        vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
        color3 bsdf = intersectedMaterial.eval( ray, intersection, wi, &scatteringPdf ) *
                      max( dot( normal, wi ), 0.f );
        if ( !isBlack( bsdf ) ) {
            if ( vis ) { Li = color3( 0 ); }
            if ( !isBlack( Li ) ) {
                float weight = PowerHeuristic( 1, lightPdf, 1, scatteringPdf );
                ret += bsdf * Li * weight / lightPdf;
            }
        }
    }

    color3 bsdf =
        intersectedMaterial.sample( ray, intersection, uScattering, &wi, &scatteringPdf ) *
        max( dot( intersection->normal, wi ), 0.f );
    if ( !isBlack( bsdf ) && scatteringPdf > 0 ) {
        float weight = 1;
        lightPdf     = light->pdf_Li( *intersection, wi );
        if ( lightPdf == 0 ) return ret;
        weight       = PowerHeuristic( 1, scatteringPdf, 1, lightPdf );
        auto ray_ref = Ray(
            intersection->position + ( wi * acne_eps ), normalize( wi ), 0, 10000, ray->depth + 1 );
        ray_ref.hasDifferentials = false;

        Intersection temp_intersection;
        auto foundIntersection = intersectKdTree( scene, tree, &ray_ref, &temp_intersection );
        color3 Li( 0.f );
        if ( foundIntersection ) {
            if ( !isBlack( scene->GetMaterial( temp_intersection.materialIndex ).m_emission ) )
                Li = static_cast<ShapeLight*>( light )->L( temp_intersection, -wi );
        }
        else
            Li = light->Le( &ray_ref );
        if ( !isBlack( Li ) ) ret += bsdf * Li * weight / scatteringPdf;
    }
    return ret;
}

Pathtracer::Pathtracer( int maxDepth, Sampler* sampler ) : Integrator( sampler, maxDepth ) {}

void Pathtracer::preprocess( Scene* scene, Sampler* sampler ) {
    m_lightDistrib = std::unique_ptr<LightDistribution> { new UniformLightDistribution( scene ) };
}

color3 Pathtracer::trace_ray( Scene* scene,
                              Ray& ray,
                              KdTree* tree,
                              Intersection* intersection,
                              Sampler* sampler ) {
    color3 ret = color3( 0, 0, 0 );
    color3 throughput( 1.f );

    // if ( ray->depth > scene->depth ) return color3( 0.f );

    int depth;

    for ( depth = 0;; depth++ ) {

        if ( intersectKdTree( scene, tree, &ray, intersection ) ) {
            intersection->hit = true;
            // Compute necessary differential information for texture filtering
            if ( ray.hasDifferentials ) intersection->computeDifferentials( &ray );
            const Material& intersectedMaterial = scene->GetMaterial( intersection->materialIndex );
            if ( !isBlack( intersectedMaterial.m_emission ) && depth == 0 ) {
                ret += throughput *
                       static_cast<ShapeLight*>( scene->lights[0] )->L( *intersection, -ray.dir );
            }

            int nLights = int( scene->lights.size() );
            if ( nLights > 0 ) {
                int lightNum;
                float lightPdf;
                const Distribution1D* lightDistrib =
                    m_lightDistrib->lookup( intersection->position );

                if ( lightDistrib ) {
                    lightNum = lightDistrib->SampleDiscrete( sampler->Get1D(), &lightPdf );
                }
                else {
                    lightNum = min( (int)( sampler->Get1D() * nLights ), nLights - 1 );
                    lightPdf = 1.f / nLights;
                }
                if ( lightPdf > 0 ) {
                    vec2 uLight      = sampler->Get2D();
                    vec2 uScattering = sampler->Get2D();
                    ret += throughput *
                           directIllumination( scene,
                                               tree,
                                               &ray,
                                               intersection,
                                               scene->lights[lightNum],
                                               uScattering,
                                               uLight ) /
                           lightPdf;
                }
            }

            float pdf;
            vec3 wi;
            color3 bsdf =
                intersectedMaterial.sample( &ray, intersection, sampler->Get2D(), &wi, &pdf );
            wi = normalize( wi );
            if ( isBlack( bsdf ) || pdf == 0.f ) break;
            if ( intersectedMaterial.m_MatType != TRANSPARENT )
                throughput *= bsdf * abs( dot( wi, intersection->normal ) ) / pdf;
            else {
                // std::cout << dot( wi, intersection->normal ) << std::endl;
            }

            ray = Ray( intersection->position + acne_eps * wi, wi, 0, 10000, 0 );
        }
        else {
            for ( auto& env : scene->envLights )
                ret += throughput * env->Le( &ray );
            break;
        }

        if ( depth >= maxDepth ) break;

        // if ( ray->tmax < 0 || !intersection->hit || ray->dir.z == 0 ) return ret;

        if ( scene->medium != nullptr )
            ret = ret * scene->medium->tr( ray, scene->ysol ) +
                  scene->medium->sample( ray, scene, tree, scene->ysol );
    }

    return ret;
}
