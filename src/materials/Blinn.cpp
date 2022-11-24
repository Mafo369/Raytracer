#include "Blinn.h"
#include "../Light.h"
#include "../Object.h"
#include "../bsdf.hpp"

#include "../sampling/sampling.h"
#define SCATTER_SAMPLES 1
#define GI_SAMPLES 1

Blinn::Blinn() {
    m_IOR             = 1.0;
    m_specularColor   = color3( 0.7f );
    m_albedo    = color3( 0.5f );
    m_shininess       = 20.f;
    m_reflection      = vec3( 0, 0, 0 );
    m_refraction      = vec3( 0, 0, 0 );
    m_absorption      = vec3( 0, 0, 0 );
    m_reflectionGloss = 0;
    m_refractionGloss = 0;
    m_emission        = color3( 0 );
}

color3 Blinn::shade( Intersection* intersection, vec3 v, Light* light, float intensity ) {
    color3 ret   = color3( 0.f );
    vec3 n       = intersection->normal;
    auto lc      = light->getColor();
    auto samples = light->getSamples();

    if ( intersection->isOutside ) {
        for ( auto& sample : samples ) {
            vec3 l;
            if ( light->isDirectional() ) { l = -sample; }
            else {
                l = -normalize( intersection->position - sample );
            }
            float LdotN = dot( l, n );
            vec3 vl     = v + l;
            vec3 h      = vl;
            h           = normalize( h );
            float NdotH = dot( n, h );
            float s     = std::pow( NdotH, m_shininess );

            if ( LdotN >= 0 ) {
                if ( m_texture != nullptr ) {
                    vec3 duv[2];
                    duv[0]        = vec3( intersection->dudx, intersection->dvdx, 0 );
                    duv[1]        = vec3( intersection->dudy, intersection->dvdy, 0 );
                    vec3 texColor = m_texture->value( intersection->u, intersection->v, duv );
                    ret += lc * LdotN * ( texColor + s * m_specularColor );
                }
                else {
                    ret += lc * LdotN * ( m_albedo + s * m_specularColor );
                }
            }
        }
        ret = ( ret / float( samples.size() ) ) * intensity;
    }

    return ret;
}

Blinn::~Blinn() {}

color3 Blinn::textureColor( float u, float v, int face ) {
    if ( m_texture != nullptr ) return m_texture->value( u, v );
    return color3( 0.f );
}

color3 Blinn::sample_f( vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type ) {
    if ( type == 0 ) { // whitted reflection
        if ( isBlack( m_reflection ) && isBlack( m_refraction ) ) return color3( 0 );
        bool isOutside = dot( wo, normal ) < 0;
        if ( !isOutside ) normal = -normal;
        *wi  = normalize( reflect( wo, normal ) );
        *pdf = 1;

        if ( !isBlack( m_refraction ) ) {
            float n1, n2;
            if ( isOutside ) {
                n1 = 1.0;
                n2 = m_IOR;
            }
            else {
                n2 = 1.0;
                n1 = m_IOR;
            }
            wo = -wo;

            // calculate refraction ray direction
            float c1 = dot( normal, wo );
            float s1 = sqrt( 1.0 - c1 * c1 );
            float s2 = n1 / n2 * s1;
            float c2 = sqrt( 1.0 - s2 * s2 );

            if ( s2 * s2 > 1.0 ) return m_refraction;

            float r0 = ( n1 - n2 ) / ( n1 + n2 );
            r0 *= r0;
            float r;
            if ( n1 <= n2 )
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 );
            else
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 );

            return m_refraction * r;
        }

        // float absCosThetar = abs(dot(*wi, normal));
        // return m_diffuseColor / absCosThetar;
        return m_reflection;
    }
    else if ( type == 1 ) { // whitted refraction
        if ( isBlack( m_refraction ) ) return color3( 0 );
        float n1, n2;
        bool isOutside = dot( normal, wo ) > 0;
        if ( isOutside ) {
            n1 = 1.0;
            n2 = m_IOR;
        }
        else {
            n2     = 1.0;
            n1     = m_IOR;
            normal = -normal;
        }

        // calculate refraction ray direction
        float c1 = dot( normal, wo );
        float s1 = sqrt( 1.0 - c1 * c1 );
        float s2 = n1 / n2 * s1;
        float c2 = sqrt( 1.0 - s2 * s2 );
        vec3 p   = normalize( wo - c1 * normal );
        vec3 pt  = s2 * -p;
        vec3 nt  = c2 * -normal;

        // cannot refract
        if ( s2 * s2 > 1.0 ) return color3( 0 );

        float r0 = ( n1 - n2 ) / ( n1 + n2 );
        r0 *= r0;
        float r;
        if ( n1 <= n2 )
            r = r0 + ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 );
        else
            r = r0 + ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 );
        float fr = 1.0 - r;

        *wi  = normalize( pt + nt );
        *pdf = 1;
        // float absCosThetar = abs(dot(*wi, normal));
        // return m_diffuseColor / absCosThetar;
        return m_refraction * fr;
    }
    return color3( 0, 0, 0 );
}

color3 Blinn::ambientColor( Ray* ray, Intersection* intersection, color3 lightColor ) {
    if ( ray->hasDifferentials && m_texture != nullptr ) {
        vec3 duv[2];
        duv[0]        = vec3( intersection->dudx, intersection->dvdx, 0 );
        duv[1]        = vec3( intersection->dudy, intersection->dvdy, 0 );
        vec3 texColor = m_texture->value( intersection->u, intersection->v, duv );
        return texColor * lightColor;
    }
    return m_albedo * lightColor;
}

vec3 Phong_BRDF( const vec3& wi, const vec3& wo, const vec3& N, float phongExponent ) {
    vec3 reflected = reflect( wo, N );
    float d        = dot( reflected, wi );
    if ( d < 0 ) return vec3( 0 );
    vec3 lobe = pow( vec3( d ), vec3( phongExponent ) ) * ( vec3( phongExponent ) + vec3( 2 ) ) /
                ( 2.f * Pi );
    return lobe;
}

vec3 random_Phong( const vec3& R, float phong_exponent ) {
    float r1      = uniform01( engine );
    float r2      = uniform01( engine );
    float facteur = sqrt( 1 - std::pow( r2, 2. / ( phong_exponent + 1 ) ) );
    vec3 direction_aleatoire_repere_local( cos( 2 * M_PI * r1 ) * facteur,
                                           sin( 2 * M_PI * r1 ) * facteur,
                                           std::pow( r2, 1. / ( phong_exponent + 1 ) ) );
    // vec3 aleatoire(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);
    // vec3 tangent1 = cross(R, aleatoire); tangent1.normalize();
    vec3 tangent1;
    vec3 absR( abs( R[0] ), abs( R[1] ), abs( R[2] ) );
    if ( absR[0] <= absR[1] && absR[0] <= absR[2] ) { tangent1 = vec3( 0, -R[2], R[1] ); }
    else if ( absR[1] <= absR[0] && absR[1] <= absR[2] ) {
        tangent1 = vec3( -R[2], 0, R[0] );
    }
    else
        tangent1 = vec3( -R[1], R[0], 0 );
    tangent1 = normalize( tangent1 );

    vec3 tangent2 = cross( tangent1, R );

    return direction_aleatoire_repere_local[2] * R +
           direction_aleatoire_repere_local[0] * tangent1 +
           direction_aleatoire_repere_local[1] * tangent2;
}

color3 Blinn::scatterColor( Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection ) {
    auto ret = color3( 0.0 );

    if ( !isBlack( m_reflection ) ) {
        vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
        vec3 dirR   = reflect( ray->dir, normal );
        Ray new_ray;
        new_ray.dox = vec3( 0 );
        new_ray.doy = vec3( 0 );
        new_ray.ddx = vec3( 0 );
        new_ray.ddy = vec3( 0 );
        rayInit( &new_ray,
                 intersection->position + ( acne_eps * normal ),
                 normalize( dirR ),
                 ray->pixel,
                 0,
                 100000,
                 ray->depth + 1 );

        Intersection temp_inter;
        ret = trace_ray( scene, &new_ray, tree, &temp_inter );
    }
    else if ( !isBlack( m_refraction ) ) {
        Ray new_ray;
        new_ray.dox = vec3( 0 );
        new_ray.doy = vec3( 0 );
        new_ray.ddx = vec3( 0 );
        new_ray.ddy = vec3( 0 );
        float n1, n2;
        n1          = 1.0;
        n2          = m_IOR;
        vec3 normal = intersection->normal;
        if ( !intersection->isOutside ) {
            n1     = m_IOR;
            n2     = 1.0;
            normal = -normal;
        }
        vec3 v = -ray->dir;

        // calculate refraction ray direction
        float c1 = dot( normal, v );
        float s1 = sqrt( 1.0 - c1 * c1 );
        float s2 = n1 / n2 * s1;
        float c2 = sqrt( 1.0 - s2 * s2 );
        vec3 p   = normalize( v - c1 * normal );
        vec3 pt  = s2 * -p;
        vec3 nt  = c2 * -normal;

        vec3 refractDir = normalize( pt + nt );
        vec3 reflDir    = normalize( reflect( ray->dir, intersection->normal ) );

        if ( s2 * s2 <= 1.0 ) {
            float r0 = ( n1 - n2 ) / ( n1 + n2 );
            r0 *= r0;
            float r;
            if ( n1 <= n2 )
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 );
            else
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 );

            if ( uniform01( engine ) < r ) {
                rayInit( &new_ray,
                         intersection->position + ( 0.01f * reflDir ),
                         normalize( reflDir ),
                         ray->pixel,
                         0,
                         1000000,
                         ray->depth + 1 );
            }
            else {
                rayInit( &new_ray,
                         intersection->position - ( 0.001f * normal ),
                         normalize( refractDir ),
                         ray->pixel,
                         0,
                         100000,
                         ray->depth + 1 );
            }
        }
        else {
            rayInit( &new_ray,
                     intersection->position + ( 0.001f * normal ),
                     normalize( reflDir ),
                     ray->pixel,
                     0,
                     100000,
                     ray->depth + 1 );
        }
        Intersection refInter;
        ret += trace_ray( scene, &new_ray, tree, &refInter );
    }
    else {
        float phongExponent = m_shininess;
        // auto sphereL = scene->lights[0];
        auto sphereL = scene->objects[scene->objects.size() - 1];
        vec3 axePO   = normalize( intersection->position - sphereL->geom.sphere.center );
        vec3 dirA    = random_dir( axePO );
        vec3 ptA     = dirA * sphereL->geom.sphere.radius + sphereL->geom.sphere.center;
        vec3 wi      = normalize( ptA - intersection->position );

        float d_light2 = distance( ptA, intersection->position );
        d_light2 *= d_light2;

        vec3 Np = dirA;

        Intersection temp_inter;
        Ray rayS;
        vec3 origin = intersection->position + ( acne_eps * wi );
        rayInit( &rayS, origin, wi, vec2( 0, 0 ), 0.f, sqrt( d_light2 ) );
        rayS.shadow = true;
        rayS.dox    = vec3( 0.f );
        rayS.doy    = vec3( 0.f );
        rayS.ddx    = vec3( 0.f );
        rayS.ddy    = vec3( 0.f );
        if ( intersectKdTree( scene, tree, &rayS, &temp_inter ) &&
             rayS.tmax * rayS.tmax < d_light2 * 0.99 ) { //|| dot(intersection->normal, wi) < 0){
            ret += color3( 0 );
        }
        else {
            color3 brdf = intersection->mat->m_albedo / Pi +
                          m_specularColor *
                              Phong_BRDF( wi, ray->dir, intersection->normal, phongExponent );
            float J   = 1. * dot( Np, -wi ) / d_light2;
            float pdf = dot( axePO, dirA ) /
                        ( Pi * sphereL->geom.sphere.radius * sphereL->geom.sphere.radius );
            if ( pdf > 0 )
                ret += sphereL->mat->m_emission / sqr( sphereL->geom.sphere.radius ) *
                       max( 0.f, dot( intersection->normal, wi ) ) * J * brdf / pdf;
        }

        // Indirect illumination
        color3 iColor      = color3( 0, 0, 0 );
        auto specularColor = m_specularColor;
        float proba        = 1. - ( specularColor.r + specularColor.g + specularColor.b ) / 3;
        vec3 R             = normalize( reflect( ray->dir, intersection->normal ) );
        bool sample_diffuse;
        if ( uniform01( engine ) < proba ) {
            sample_diffuse = true;
            dirA           = random_dir( intersection->normal );
        }
        else {
            sample_diffuse = false;
            dirA           = random_Phong( R, phongExponent );
        }
        if ( dot( dirA, intersection->normal ) < 0 ) return ret;
        if ( dot( dirA, R ) < 0 ) return ret;

        Ray ray_ref;
        rayInit( &ray_ref,
                 intersection->position + ( acne_eps * dirA ),
                 normalize( dirA ),
                 ray->pixel,
                 0,
                 100000,
                 ray->depth + 1 );
        ray_ref.dox = vec3( 0 );
        ray_ref.doy = vec3( 0 );
        ray_ref.ddx = vec3( 0 );
        ray_ref.ddy = vec3( 0 );

        float pdfPhong =
            ( phongExponent + 1.f ) / ( 2.f * Pi ) * pow( dot( R, dirA ), phongExponent );
        float pdf = proba * dot( intersection->normal, dirA ) / ( Pi ) + ( 1.f - proba ) * pdfPhong;
        if ( pdf <= 0 ) return ret;

        Intersection temp_intersection;
        auto reflColor = trace_ray( scene, &ray_ref, tree, &temp_intersection, false );

        if ( sample_diffuse )
            iColor += reflColor * intersection->mat->m_albedo *
                      ( dot( intersection->normal, dirA ) / ( Pi ) / pdf );
        else
            iColor += reflColor *
                      ( dot( intersection->normal, dirA ) *
                        Phong_BRDF( dirA, ray->dir, intersection->normal, phongExponent ) / pdf ) *
                      m_specularColor;

        ret += iColor;
    }

    return ret;
}

color3 Blinn::reflectionColor( Scene* scene,
                               KdTree* tree,
                               Ray* ray,
                               Intersection* intersection,
                               color3& color,
                               vec3 normal ) {
    vec3 r = normalize( reflect( ray->dir, normal ) );
    Ray ray_ref;
    rayInit( &ray_ref,
             intersection->position + ( acne_eps * r ),
             r,
             ray->pixel,
             0,
             100000,
             ray->depth + 1 );

    if ( ray->hasDifferentials ) {
        vec3 wo = ray->dir;
        vec3 wi = r;
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

    if ( !temp_intersection.hit ) {
        if ( scene->m_skyTexture != nullptr ) {
            vec3 dir = r;
            float z  = asin( -dir.z ) / float( M_PI ) + 0.5;
            float x  = dir.x / ( abs( dir.x ) + abs( dir.y ) );
            float y  = dir.y / ( abs( dir.x ) + abs( dir.y ) );
            point3 p = point3( 0.5, 0.5, 0.0 ) +
                       z * ( x * point3( 0.5, 0.5, 0.0 ) + y * point3( -0.5, 0.5, 0.0 ) );
            // TODO: Multiply with intensity var
            color3 env = 0.7f * scene->m_skyTexture->value( p.x, p.y );
            color += m_reflection * env;
        }
    }
    else if ( intersection->isOutside ) {
        color += ( reflColor * m_reflection );
    }
    return reflColor;
}

color3 Blinn::refractionColor( Scene* scene,
                               KdTree* tree,
                               Ray* ray,
                               Intersection* intersection,
                               color3 reflectionShade,
                               vec3 normal ) {
    float n1, n2;
    color3 refractionShade = color3( 0, 0, 0 );
    if ( intersection->isOutside ) {
        n1 = 1.0;
        n2 = m_IOR;
    }
    else {
        n2 = 1.0;
        n1 = m_IOR;
    }
    vec3 v = -ray->dir;

    // calculate refraction ray direction
    float c1 = dot( normal, v );
    float s1 = sqrt( 1.0 - c1 * c1 );
    float s2 = n1 / n2 * s1;
    float c2 = sqrt( 1.0 - s2 * s2 );
    vec3 p   = normalize( v - c1 * normal );
    vec3 pt  = s2 * -p;
    vec3 nt  = c2 * -normal;

    vec3 refractDir = normalize( pt + nt );

    if ( s2 * s2 <= 1.0 ) {
        Ray ray_refr;
        rayInit( &ray_refr,
                 intersection->position + ( acne_eps * refractDir ),
                 refractDir,
                 ray->pixel,
                 0,
                 100000,
                 ray->depth + 1 );
        if ( ray->hasDifferentials ) {
            vec3 dndx =
                intersection->dn[0] * intersection->dudx + intersection->dn[1] * intersection->dvdx;
            vec3 dndy =
                intersection->dn[0] * intersection->dudy + intersection->dn[1] * intersection->dvdy;
            float eta = 1.f / m_IOR;
            if ( !intersection->isOutside ) {
                eta  = 1.f / m_IOR;
                dndx = -dndx;
                dndy = -dndy;
            }
            vec3 wo = -ray->dir;
            vec3 wi = refractDir;

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
        color3 refrColor = trace_ray( scene, &ray_refr, tree, &temp_inter );

        if ( temp_inter.hit ) {
            // Schlick's approximation for transmittance vs. reflectance
            float r0 = ( n1 - n2 ) / ( n1 + n2 );
            r0 *= r0;
            float r;
            if ( n1 <= n2 )
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 );
            else
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 );
            float t = 1.0 - r;

            // compute total refraction color
            color3 refractionColor = m_refraction * ( t * refrColor + r * reflectionShade );

            if ( !temp_inter.isOutside ) {
                refractionColor.r *= exp( -m_absorption.r * ray_refr.tmax );
                refractionColor.g *= exp( -m_absorption.g * ray_refr.tmax );
                refractionColor.b *= exp( -m_absorption.b * ray_refr.tmax );
            }
            refractionShade += refractionColor;
        }
        else {
            if ( scene->m_skyTexture != nullptr ) {
                vec3 dir = refractDir;
                float z  = asin( -dir.z ) / float( M_PI ) + 0.5;
                float x  = dir.x / ( abs( dir.x ) + abs( dir.y ) + 0.00001 );
                float y  = dir.y / ( abs( dir.x ) + abs( dir.y ) + 0.00001 );
                point3 p = point3( 0.5, 0.5, 0.0 ) +
                           z * ( x * point3( 0.5, 0.5, 0.0 ) + y * point3( -0.5, 0.5, 0.0 ) );
                color3 env = 0.7f * scene->m_skyTexture->value( p.x, p.y );
                refractionShade += env;
            }
        }
    }
    else {
        refractionShade += m_refraction * reflectionShade;
    }
    return refractionShade;
}
