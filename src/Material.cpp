#include "Material.h"
#include "Light.h"
#include "bsdf.hpp"

inline float saturate( float x ) {
    return clamp( x, 0.0f, 1.0f );
}

bool computeBrdfData( BrdfData& data,
                      const vec3& v,
                      const vec3& wi,
                      const vec3& n,
                      const vec2& uv,
                      const color3& baseColor,
                      const float& metalness,
                      const float& roughness ) {
    vec3 vl = v + wi;
    vec3 h  = vl;
    h       = normalize( h );

    float LdotN = dot( wi, n );
    float VdotN = dot( v, n );

    // Backfacing
    if ( LdotN <= 0.0f || VdotN <= 0.0f ) return false;

    data.LdotN = min( max( 0.0001f, LdotN ), 1.0f );
    data.VdotN = min( max( 0.0001f, VdotN ), 1.0f );

    data.LdotH = saturate( dot( wi, h ) );
    data.NdotH = saturate( dot( n, h ) );
    data.VdotH = saturate( dot( v, h ) );
    data.uv    = uv;

    data.specularF0         = baseColorToSpecularF0( baseColor, metalness );
    data.diffuseReflectance = baseColorToDiffuseReflectance( baseColor, metalness );

    data.roughness = roughness;
    data.alpha     = data.roughness;
    data.alphaSq   = data.alpha * data.alpha;

    // Pre-calculate some more BRDF terms
    // data.F = color3(RDM_Fresnel( data.LdotH, 1.f, 2.0 ));
    data.F         = evalFresnel( data.specularF0, shadowedF90( data.specularF0 ), data.LdotH );
    data.metalness = metalness;

    return true;
};

color3 specularReflect( Ray* ray,
                        Intersection* intersection,
                        Scene* scene,
                        KdTree* tree,
                        Sampler* sampler ) {
    float pdf;
    vec3 wi;
    int type = 0;

    // BRDF
    color3 f = intersection->mat->sample_f(
        ray->dir, &wi, intersection->normal, point2( 0, 0 ) /*sampler->Get2D()*/, &pdf, type );

    if ( pdf > 0 && !isBlack( f ) ) {
        vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
        Ray ray_ref;
        rayInit( &ray_ref,
                 intersection->position + ( acne_eps * wi ),
                 wi,
                 ray->pixel,
                 0,
                 100000,
                 ray->depth + 1 );

        if ( ray->hasDifferentials ) {
            vec3 wo = ray->dir;
            vec3 dndx =
                intersection->dn[0] * intersection->dudx + intersection->dn[1] * intersection->dvdx;
            vec3 dndy =
                intersection->dn[0] * intersection->dudy + intersection->dn[1] * intersection->dvdy;
            vec3 dwodx = -ray->ddx - wo, dwody = -ray->ddy - wo;
            float dDNdx = dot( dwodx, normal ) + dot( wo, dndx );
            float dDNdy = dot( dwody, normal ) + dot( wo, dndy );
            vec3 ddx    = wi - dwodx + 2.f * vec3( dot( wo, normal ) * dndx + dDNdx * normal );
            vec3 ddy    = wi - dwody + 2.f * vec3( dot( wo, normal ) * dndy + dDNdy * normal );

            ray_ref.dox = ( intersection->position + ( acne_eps * ddx ) ) + intersection->dpdx;
            ray_ref.doy = ( intersection->position + ( acne_eps * ddy ) ) + intersection->dpdy;

            ray_ref.ddx = normalize( ddx );
            ray_ref.ddy = normalize( ddy );
        }

        Intersection temp_intersection;
        auto reflColor = trace_ray( scene, &ray_ref, tree, &temp_intersection );

        color3 Li = color3( 0 );

        // if ( !temp_intersection.hit ) {
        //    if ( scene->m_skyTexture != nullptr ) {
        //        vec3 dir = wi;
        //        float z  = asin( -dir.z ) / float( M_PI ) + 0.5;
        //        float x  = dir.x / ( abs( dir.x ) + abs( dir.y ) );
        //        float y  = dir.y / ( abs( dir.x ) + abs( dir.y ) );
        //        point3 p = point3( 0.5, 0.5, 0.0 ) +
        //                   z * ( x * point3( 0.5, 0.5, 0.0 ) + y * point3( -0.5, 0.5, 0.0 ) );
        //        // TODO: Multiply with intensity var
        //        color3 env = 0.7f * scene->m_skyTexture->value( p.x, p.y );
        //        Li         = env;
        //    }
        //}
        // else {
        Li = reflColor;
        //}

        return f * Li;
    }

    return color3( 0 );
}

color3 specularTransmission( Ray* ray,
                             Intersection* intersection,
                             Scene* scene,
                             KdTree* tree,
                             Sampler* sampler ) {
    float pdf;
    vec3 wi;
    int type = 1;

    // BTDF
    color3 f = intersection->mat->sample_f(
        -ray->dir, &wi, intersection->normal, point2( 0, 0 ) /*sampler->Get2D()*/, &pdf, type );

    if ( pdf > 0 && !isBlack( f ) ) {
        vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
        Ray ray_refr;
        rayInit( &ray_refr,
                 intersection->position + ( acne_eps * wi ),
                 wi,
                 ray->pixel,
                 0,
                 100000,
                 ray->depth + 1 );
        if ( ray->hasDifferentials ) {
            vec3 dndx =
                intersection->dn[0] * intersection->dudx + intersection->dn[1] * intersection->dvdx;
            vec3 dndy =
                intersection->dn[0] * intersection->dudy + intersection->dn[1] * intersection->dvdy;
            float eta = 1.f / intersection->mat->m_IOR;
            if ( !intersection->isOutside ) {
                eta  = 1.f / intersection->mat->m_IOR;
                dndx = -dndx;
                dndy = -dndy;
            }
            vec3 wo = -ray->dir;

            vec3 dwodx = -ray->ddx - wo, dwody = -ray->ddy - wo;
            float dDNdx = dot( dwodx, normal ) + dot( wo, dndx );
            float dDNdy = dot( dwody, normal ) + dot( wo, dndy );

            float mu = eta * dot( wo, normal ) - abs( dot( wi, normal ) );
            float dmudx =
                ( eta - ( eta * eta * dot( wo, normal ) ) / abs( dot( wi, normal ) ) ) * dDNdx;
            float dmudy =
                ( eta - ( eta * eta * dot( wo, normal ) ) / abs( dot( wi, normal ) ) ) * dDNdy;

            vec3 ddx = normalize( wi - eta * dwodx + vec3( mu * dndx + dmudx * normal ) );
            vec3 ddy = normalize( wi - eta * dwody + vec3( mu * dndy + dmudy * normal ) );

            ray_refr.dox = ( intersection->position + ( acne_eps * ddx ) ) + intersection->dpdx;
            ray_refr.doy = ( intersection->position + ( acne_eps * ddy ) ) + intersection->dpdy;
            ray_refr.ddx = ddx;
            ray_refr.ddy = ddy;
        }

        Intersection temp_inter;
        auto refrColor = trace_ray( scene, &ray_refr, tree, &temp_inter );

        color3 Li = color3( 0 );

        if ( temp_inter.hit ) {
            // if(!temp_inter.isOutside){
            //  refractionColor.r *= exp(-m_absorption.r * ray_refr.tmax);
            //  refractionColor.g *= exp(-m_absorption.g * ray_refr.tmax);
            //  refractionColor.b *= exp(-m_absorption.b * ray_refr.tmax);
            //}
            Li = refrColor;
        }
        else {
            // if ( scene->m_skyTexture != nullptr ) {
            //    vec3 dir = wi;
            //    float z  = asin( -dir.z ) / float( M_PI ) + 0.5;
            //    float x  = dir.x / ( abs( dir.x ) + abs( dir.y ) + 0.00001 );
            //    float y  = dir.y / ( abs( dir.x ) + abs( dir.y ) + 0.00001 );
            //    point3 p = point3( 0.5, 0.5, 0.0 ) +
            //               z * ( x * point3( 0.5, 0.5, 0.0 ) + y * point3( -0.5, 0.5, 0.0 ) );
            //    color3 env = 0.7f * scene->m_skyTexture->value( p.x, p.y );
            //    return env;
            //}
        }

        return f * Li;
    }

    return color3( 0 );
}
