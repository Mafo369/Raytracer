#include "Material.h"
#include "Light.h"
#include "bsdf.hpp"
#include "defines.h"
#include "integrator.h"
#include "sampling/sampling.h"

#define GI 1
#define SCATTER_SAMPLES 1
#define GI_SAMPLES 1

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
    else if ( absR[1] <= absR[0] && absR[1] <= absR[2] ) { tangent1 = vec3( -R[2], 0, R[0] ); }
    else
        tangent1 = vec3( -R[1], R[0], 0 );
    tangent1 = normalize( tangent1 );

    vec3 tangent2 = cross( tangent1, R );

    return direction_aleatoire_repere_local[2] * R +
           direction_aleatoire_repere_local[0] * tangent1 +
           direction_aleatoire_repere_local[1] * tangent2;
}

vec3 random_Cook( vec3 r ) {
    return vec3( 0 );
}
inline float rsqrt( float x ) {
    return inversesqrt( x );
}

// Samples a microfacet normal for the Beckmann distribution using walter's method.
// Source: "Microfacet Models for Refraction through Rough Surfaces" by Walter et al.
// PDF is 'D * NdotH'
vec3 sampleBeckmannWalter( vec3 Vlocal, vec2 alpha2D, vec2 u ) {
    float alpha = dot( alpha2D, vec2( 0.5f, 0.5f ) );

    // Equations (28) and (29) from Walter's paper for Beckmann distribution
    float tanThetaSquared = -( alpha * alpha ) * log( 1.0f - u.x );
    float phi             = ( 2.f * Pi ) * u.y;

    // Calculate cosTheta and sinTheta needed for conversion to H vector
    float cosTheta = rsqrt( 1.0f + tanThetaSquared );
    float sinTheta = sqrt( 1.0f - cosTheta * cosTheta );

    // Convert sampled spherical coordinates to H vector
    return normalize( vec3( sinTheta * cos( phi ), sinTheta * sin( phi ), cosTheta ) );
}

// Calculates rotation quaternion from input vector to the vector (0, 0, 1)
// Input vector must be normalized!
vec4 getRotationToZAxis( vec3 input ) {

    // Handle special case when input is exact or near opposite of (0, 0, 1)
    if ( input.z < -0.99999f ) return vec4( 1.0f, 0.0f, 0.0f, 0.0f );

    return normalize( vec4( input.y, -input.x, 0.0f, 1.0f + input.z ) );
}

// Samples a direction within a hemisphere oriented along +Z axis with a cosine-weighted
// distribution Source: "Sampling Transformations Zoo" in Ray Tracing Gems by Shirley et al.
vec3 sampleHemisphereCook( vec2 u, float& pdf ) {

    float a = sqrt( u.x );
    float b = ( 2 * Pi ) * u.y;

    vec3 result = vec3( a * cos( b ), a * sin( b ), sqrt( 1.0f - u.x ) );

    pdf = result.z * ( InvPi );

    return result;
}

vec3 sampleHemisphereCook( vec2 u ) {
    float pdf;
    return sampleHemisphereCook( u, pdf );
}

// Calculates rotation quaternion from vector (0, 0, 1) to the input vector
// Input vector must be normalized!
vec4 getRotationFromZAxis( vec3 input ) {

    // Handle special case when input is exact or near opposite of (0, 0, 1)
    if ( input.z < -0.99999f ) return vec4( 1.0f, 0.0f, 0.0f, 0.0f );

    return normalize( vec4( -input.y, input.x, 0.0f, 1.0f + input.z ) );
}

// Returns the quaternion with inverted rotation
vec4 invertRotation( vec4 q ) {
    return vec4( -q.x, -q.y, -q.z, q.w );
}

// Optimized point rotation using quaternion
// Source: https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
vec3 rotatePoint( vec4 q, vec3 v ) {
    const vec3 qAxis = vec3( q.x, q.y, q.z );
    return 2.0f * dot( qAxis, v ) * qAxis + ( q.w * q.w - dot( qAxis, qAxis ) ) * v +
           2.0f * q.w * cross( qAxis, v );
}

float Smith_G1_Beckmann_Walter( float a ) {
    if ( a < 1.6f ) {
        return ( ( 3.535f + 2.181f * a ) * a ) / ( 1.0f + ( 2.276f + 2.577f * a ) * a );
    }
    else { return 1.0f; }
}
// Function to calculate 'a' parameter for lambda functions needed in Smith G term
// This is a version for shape invariant (isotropic) NDFs
// Note: makse sure NdotS is not negative
float Smith_G_a( float alpha, float NdotS ) {
    return NdotS / ( max( 0.00001f, alpha ) * sqrt( 1.0f - min( 0.99999f, NdotS * NdotS ) ) );
}

float Smith_G1_Beckmann_Walter( float alpha, float NdotS, float alphaSquared, float NdotSSquared ) {
    return Smith_G1_Beckmann_Walter( Smith_G_a( alpha, NdotS ) );
}

// Smith G2 term (masking-shadowing function)
// Separable version assuming independent (uncorrelated) masking and shadowing, uses G1 functions
// for selected NDF
float Smith_G2_Separable( float alpha, float NdotL, float NdotV ) {
    float aL = Smith_G_a( alpha, NdotL );
    float aV = Smith_G_a( alpha, NdotV );
    return Smith_G1_Beckmann_Walter( aL ) * Smith_G1_Beckmann_Walter( aV );
}

// Weight for the reflection ray sampled from Beckmann distribution using Walter's method
float specularSampleWeightBeckmannWalter( float alpha,
                                          float alphaSquared,
                                          float NdotL,
                                          float NdotV,
                                          float HdotL,
                                          float NdotH ) {
    return ( HdotL * Smith_G2_Separable( alpha, NdotL, NdotV ) ) / ( NdotV * NdotH );
}

// Samples a reflection ray from the rough surface using selected microfacet distribution and
// sampling method Resulting weight includes multiplication by cosine (NdotL) term
vec3 sampleSpecularMicrofacet( vec3 Vlocal,
                               float alpha,
                               float alphaSquared,
                               vec3 specularF0,
                               vec2 u,
                               vec3& weight ) {

    // Sample a microfacet normal (H) in local space
    vec3 Hlocal;
    if ( alpha == 0.0f ) {
        // Fast path for zero roughness (perfect reflection), also prevents NaNs appearing due to
        // divisions by zeroes
        Hlocal = vec3( 0.0f, 0.0f, 1.0f );
    }
    else {
        // For non-zero roughness, this calls VNDF sampling for GG-X distribution or Walter's
        // sampling for Beckmann distribution
        Hlocal = sampleBeckmannWalter( Vlocal, vec2( alpha, alpha ), u );
    }

    // Reflect view direction to obtain light vector
    vec3 Llocal = reflect( -Vlocal, Hlocal );

    // Note: HdotL is same as HdotV here
    // Clamp dot products here to small value to prevent numerical instability. Assume that rays
    // incident from below the hemisphere have been filtered
    float HdotL       = max( 0.00001f, min( 1.0f, dot( Hlocal, Llocal ) ) );
    const vec3 Nlocal = vec3( 0.0f, 0.0f, 1.0f );
    float NdotL       = max( 0.00001f, min( 1.0f, dot( Nlocal, Llocal ) ) );
    float NdotV       = max( 0.00001f, min( 1.0f, dot( Nlocal, Vlocal ) ) );
    float NdotH       = max( 0.00001f, min( 1.0f, dot( Nlocal, Hlocal ) ) );
    vec3 F            = evalFresnel( specularF0, shadowedF90( specularF0 ), HdotL );

    // Calculate weight of the sample specific for selected sampling method
    // (this is microfacet BRDF divided by PDF of sampling method - notice how most terms cancel
    // out)
    weight =
        F * specularSampleWeightBeckmannWalter( alpha, alphaSquared, NdotL, NdotV, HdotL, NdotH );

    return Llocal;
}

void createCoordinateSystem( const vec3& N, vec3& Nt, vec3& Nb ) {
    if ( std::fabs( N.x ) > std::fabs( N.y ) )
        Nt = vec3( N.z, 0, -N.x ) / sqrtf( N.x * N.x + N.z * N.z );
    else
        Nt = vec3( 0, -N.z, N.y ) / sqrtf( N.y * N.y + N.z * N.z );
    Nb = cross( N, Nt );
}

vec3 uniformSampleHemisphere( const float& r1, const float& r2 ) {
    // cos(theta) = r1 = y
    // cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
    float sinTheta = sqrtf( 1 - r1 * r1 );
    float phi      = 2 * M_PI * r2;
    float x        = sinTheta * cosf( phi );
    float z        = sinTheta * sinf( phi );
    return vec3( x, r1, z );
}

float Beckmann_D( float alphaSquared, float NdotH ) {
    float cos2Theta   = NdotH * NdotH;
    float numerator   = exp( ( cos2Theta - 1.0f ) / ( alphaSquared * cos2Theta ) );
    float denominator = Pi * alphaSquared * cos2Theta * cos2Theta;
    return numerator / denominator;
}

// PDF of sampling a reflection vector L using 'sampleBeckmannWalter'.
// Note that PDF of sampling given microfacet normal is (D * NdotH). Remaining terms (1.0f / (4.0f *
// LdotH)) are specific for reflection case, and come from multiplying PDF by jacobian of reflection
// operator
float sampleBeckmannWalterReflectionPdf( float alpha,
                                         float alphaSquared,
                                         float NdotH,
                                         float NdotV,
                                         float LdotH ) {
    NdotH = max( 0.00001f, NdotH );
    LdotH = max( 0.00001f, LdotH );
    return Beckmann_D( max( 0.00001f, alphaSquared ), NdotH ) * NdotH / ( 4.0f * LdotH );
}

Material::Material( size_t UID, MaterialModel matModel, MatType matType ) :
    m_MatModel( matModel ), m_UID( UID ) {
    if ( matModel == MaterialModel::COOK_TORRANCE ) {
        m_IOR       = 1.0;
        m_roughness = 0.1;
        m_metalness = 0.001f;
        m_albedo    = color3( 1.f );
        m_MatType   = matType;
    }
    else if ( matModel == MaterialModel::BLINN ) {
        m_IOR             = 1.0;
        m_specularColor   = color3( 0.7f );
        m_albedo          = color3( 0.5f );
        m_shininess       = 20.f;
        m_reflection      = vec3( 0, 0, 0 );
        m_refraction      = vec3( 0, 0, 0 );
        m_absorption      = vec3( 0, 0, 0 );
        m_reflectionGloss = 0;
        m_refractionGloss = 0;
        m_emission        = color3( 0 );
    }
}

color3 Material::f( const vec3& wo, const vec3& wi, const vec3& n ) {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        float LdotN = abs( dot( wi, n ) );
        vec3 vl     = wo + wi;
        vec3 h      = vl;
        if ( h.x == 0 && h.y == 0 && h.z == 0 ) return color3( 0 );
        h           = normalize( h );
        float LdotH = dot( wi, h );
        float NdotH = dot( n, h );
        float VdotH = dot( wo, h );
        float VdotN = abs( dot( wo, n ) );
        return color3( 0 );
        // return RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, m_texture,
        //           m_diffuseColor, m_specularColor, m_roughness, m_IOR, 0, 0, -1) * LdotN;
    }
    else if ( m_MatModel == MaterialModel::BLINN ) { return color3( 0.0 ); }
    return color3( 0.0 );
}

color3 Material::sample_f( vec3 wo, vec3* wi, vec3 normal, const point2& u, float* pdf, int type ) {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        if ( type == 0 ) {
            bool isOutside = dot( wo, normal ) < 0;
            normal         = isOutside ? normal : -normal;
            *wi            = reflect( wo, normal );
            *pdf           = 1;

            vec3 v      = wo * -1.0f;
            vec3 h      = v + *wi;
            h           = h / length( h );
            float LdotH = dot( *wi, h );

            float fr =
                isOutside ? RDM_Fresnel( LdotH, 1.f, m_IOR ) : RDM_Fresnel( LdotH, m_IOR, 1.f );

            return color3( 0 );
            // return fr * m_specularColor;
        }
        else if ( type == 1 ) {
            if ( m_MatType == TRANSPARENT ) {
                // REFRACTION + REFLECTION

                bool isOutside = dot( -wo, normal ) < 0;
                normal         = isOutside ? normal : -normal;

                vec3 r = reflect( -wo, normal );

                vec3 h      = wo + r;
                h           = h / length( h );
                float LdotH = dot( r, h );

                float refractionRatio = isOutside ? ( 1.f / m_IOR ) : m_IOR;
                float fr =
                    isOutside ? RDM_Fresnel( LdotH, 1.f, m_IOR ) : RDM_Fresnel( LdotH, m_IOR, 1.f );

                vec3 unit_direction = normalize( -wo );

                if ( fr < 1.0f ) {
                    *wi  = refract( unit_direction, normal, refractionRatio );
                    *pdf = 1;
                    return vec3( 1.0f - fr );
                }
            }
        }
        else if ( type == 2 ) {
            vec3 wh = Beckmann_Sample_wh( wo, u, m_roughness, normal );
            *wi     = normalize( reflect( wo, wh ) );

            float LdotH = dot( *wi, wh );
            float LdotN = dot( *wi, normal );
            *pdf =
                BeckmannPdf( wo, wh, normal, m_roughness, LdotH, LdotN ) / ( 4.f * dot( wo, wh ) );

            return f( wo, *wi, normal );
        }
        return color3( 0 );
    }
    else if ( m_MatModel == MaterialModel::BLINN ) {
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
                    r = r0 + ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) *
                                 ( 1 - c1 );
                else
                    r = r0 + ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) *
                                 ( 1 - c2 );

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
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 );
            else
                r = r0 +
                    ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 );
            float fr = 1.0 - r;

            *wi  = normalize( pt + nt );
            *pdf = 1;
            // float absCosThetar = abs(dot(*wi, normal));
            // return m_diffuseColor / absCosThetar;
            return m_refraction * fr;
        }
        return color3( 0, 0, 0 );
    }
    return color3( 0.0 );
}

float Material::pdf( const vec3& wo, const vec3& wi, const vec3& n ) {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        // if(!SameHemisphere(wo, wi)) return 0;
        vec3 wh     = normalize( wo + wi );
        float LdotH = dot( wi, wh );
        float LdotN = dot( wi, n );
        return BeckmannPdf( wo, wh, n, m_roughness, LdotH, LdotN ) / ( 4.f * dot( wo, wh ) );
    }
    else if ( m_MatModel == MaterialModel::BLINN ) { return 0.0; }
    return 0.0;
}

color3 Material::shade( Intersection* intersection, vec3 v, Light* light, float intensity ) {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        color3 ret           = color3( 0.f );
        vec3 n               = intersection->normal;
        vec3 intersectionPos = intersection->position;
        bool outside         = intersection->isOutside;
        auto lc              = light->getColor();
        auto samples         = light->getSamples();

        //! \todo compute bsdf, return the shaded color taking into account the
        //! lightcolor

        if ( m_MatType == TRANSPARENT ) {
            for ( auto& sample : samples ) {
                vec3 lp     = sample - intersectionPos;
                vec3 l      = lp / length( lp );
                float LdotN = abs( dot( l, n ) );
                if ( LdotN == 0.0f ) continue;
                float VdotN = abs( dot( v, n ) );
                float extIOR, intIOR;
                if ( outside ) {
                    extIOR = 1.f;
                    intIOR = m_IOR;
                }
                else {
                    extIOR = m_IOR;
                    intIOR = 1.f;
                }

                // REFLECTION
                vec3 hr     = ( v + l );
                hr          = hr / length( hr );
                float LdotH = abs( dot( l, hr ) );
                float NdotH = abs( dot( n, hr ) );
                float VdotH = abs( dot( v, hr ) );
                // auto brdfColor = RDM_brdf(
                //    LdotH, NdotH, VdotH, LdotN, VdotN, extIOR, intIOR, m_roughness,
                //    m_specularColor );

                // REFRACTION
                hr    = -extIOR * l - intIOR * v;
                hr    = hr / length( hr );
                LdotH = abs( dot( l, hr ) );
                NdotH = abs( dot( n, hr ) );
                VdotH = abs( dot( v, hr ) );
                // auto btdfColor = RDM_btdf(
                //    LdotH, NdotH, VdotH, LdotN, VdotN, extIOR, intIOR, m_roughness,
                //    m_specularColor );

                //// BSDF
                // ret += ( brdfColor + btdfColor ) * LdotN;
            }
            ret = lc * ( ret / float( samples.size() ) ) * intensity;
        }
        else {
            for ( auto& sample : samples ) {
                vec3 lp     = sample - intersectionPos;
                vec3 l      = lp / length( lp );
                float LdotN = abs( dot( l, n ) );
                vec3 vl     = v + l;
                vec3 h      = vl;
                if ( h.x == 0 && h.y == 0 && h.z == 0 ) continue;
                h           = normalize( h );
                float LdotH = dot( l, h );
                float NdotH = dot( n, h );
                float VdotH = dot( v, h );
                float VdotN = abs( dot( v, n ) );
                // ret += lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, m_texture,
                //       m_diffuseColor, m_specularColor, m_roughness, m_IOR, uTex, vTex, face) *
                //       LdotN;
            }
            ret = ( ret / float( samples.size() ) ) * intensity;
        }

        return ret;
    }
    else if ( m_MatModel == MaterialModel::BLINN ) {
        color3 ret   = color3( 0.f );
        vec3 n       = intersection->normal;
        auto lc      = light->getColor();
        auto samples = light->getSamples();

        if ( intersection->isOutside ) {
            for ( auto& sample : samples ) {
                vec3 l;
                if ( light->isDirectional() ) { l = -sample; }
                else { l = -normalize( intersection->position - sample ); }
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
                    else { ret += lc * LdotN * ( m_albedo + s * m_specularColor ); }
                }
            }
            ret = ( ret / float( samples.size() ) ) * intensity;
        }

        return ret;
    }
    return color3( 0.0 );
}

color3 Material::textureColor( float u, float v, int face ) {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) { return m_texture->value( u, v, face ); }
    else if ( m_MatModel == MaterialModel::BLINN ) {
        if ( m_texture != nullptr ) return m_texture->value( u, v );
        return color3( 0.f );
    }
    return color3( 0.0 );
}

color3 Material::ambientColor( Ray* ray, Intersection* intersection, color3 lightColor ) {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        if ( m_texture != nullptr )
            return m_texture->value( intersection->u, intersection->v ) * lightColor;
        return m_albedo * lightColor;
    }
    else if ( m_MatModel == MaterialModel::BLINN ) {
        if ( ray->hasDifferentials && m_texture != nullptr ) {
            vec3 duv[2];
            duv[0]        = vec3( intersection->dudx, intersection->dvdx, 0 );
            duv[1]        = vec3( intersection->dudy, intersection->dvdy, 0 );
            vec3 texColor = m_texture->value( intersection->u, intersection->v, duv );
            return texColor * lightColor;
        }
        return m_albedo * lightColor;
    }
    return color3( 0.0 );
}

color3
Material::scratchAPixelScatter( Ray* ray, Scene* scene, KdTree* tree, Intersection* intersection ) {
    color3 ret( 0 );
    vec3 Nt, Nb;
    createCoordinateSystem( intersection->normal, Nt, Nb );
    float r1    = uniform01( engine );
    float r2    = uniform01( engine );
    vec3 sample = uniformSampleHemisphere( r1, r2 );
    vec3 wi( sample.x * Nb.x + sample.y * intersection->normal.x + sample.z * Nt.x,
             sample.x * Nb.y + sample.y * intersection->normal.y + sample.z * Nt.y,
             sample.x * Nb.z + sample.y * intersection->normal.z + sample.z * Nt.z );

    color3 directL = color3( 0 );

    Intersection temp_inter;
    vec3 origin = intersection->position + ( acne_eps * wi );
    Ray rayS    = Ray( origin, wi, 0.f, 10000 );
    rayS.shadow = true;
    rayS.dox    = vec3( 0.f );
    rayS.doy    = vec3( 0.f );
    rayS.ddx    = vec3( 0.f );
    rayS.ddy    = vec3( 0.f );
    if ( intersectKdTree(
             scene, tree, &rayS, &temp_inter ) ) { //|| dot(intersection->normal, wi) < 0){
        ret += color3( 0 );
    }
    else {
        color3 lc = scene->sky->getRadiance( rayS );
        directL   = r1 * lc;
    }
    color3 indirectL = color3( 0 );
#ifdef GI
    float r1I    = uniform01( engine );
    float r2I    = uniform01( engine );
    vec3 sampleI = uniformSampleHemisphere( r1I, r2I );
    vec3 wiI( sampleI.x * Nb.x + sampleI.y * intersection->normal.x + sampleI.z * Nt.x,
              sampleI.x * Nb.y + sampleI.y * intersection->normal.y + sampleI.z * Nt.y,
              sampleI.x * Nb.z + sampleI.y * intersection->normal.z + sampleI.z * Nt.z );

    Ray ray_ref = Ray(
        intersection->position + ( wiI * acne_eps ), normalize( wiI ), 0, 10000, ray->depth + 1 );
    ray_ref.hasDifferentials = true;

    Intersection temp_intersection;
    auto reflColor = trace_ray( scene, &ray_ref, tree, &temp_intersection );
    // float pdf = 1.f / (2.f * Pi);
    // indirectL = (r1I * reflColor / pdf);
    indirectL = r1I * reflColor;
#endif

    ret = m_albedo * ( directL + indirectL );
    // ret = ((m_albedo / Pi) * directL / (1.f / (2.f * Pi))) + indirectL * (m_albedo / Pi);
    // ret /= 2.f;
    return ret;
}

color3 Material::myScatter( Ray* ray, Scene* scene, KdTree* tree, Intersection* intersection ) {
    color3 ret( 0 );
    // auto sphereL = scene->objects[scene->objects.size() - 1];
    ////vec3 axePO   = normalize( intersection->position - sphereL->geom.sphere.center );
    // vec3 axePO = intersection->normal;
    // vec3 dirA    = random_dir( axePO );
    ////vec3 ptA     = dirA * sphereL->geom.sphere.radius + sphereL->geom.sphere.center;
    ////vec3 wi      = normalize( ptA - intersection->position );
    // vec3 wi = normalize(dirA);

    // vec3 Np        = dirA;
    ////float d_light2 = distance( ptA, intersection->position );
    ////d_light2 *= d_light2;

    // vec3 Nt, Nb;
    // createCoordinateSystem(intersection->normal, Nt, Nb);
    // float r1 = uniform01(engine);
    // float r2 = uniform01(engine);
    // vec3 sample = uniformSampleHemisphere(r1, r2);
    // vec3 wi(
    // sample.x * Nb.x + sample.y * intersection->normal.x + sample.z * Nt.x,
    // sample.x * Nb.y + sample.y * intersection->normal.y + sample.z * Nt.y,
    // sample.x * Nb.z + sample.y * intersection->normal.z + sample.z * Nt.z);

    // vec3 wi = random_dir(intersection->normal);
    //
    vec3 V = -ray->dir;
    // Ignore incident ray coming from "below" the hemisphere
    if ( !intersection->isOutside ) return ret;

    // Transform view direction into local space of our sampling routines
    // (local space is oriented so that its positive Z axis points along the shading normal)
    vec4 qRotationToZ   = getRotationToZAxis( intersection->normal );
    vec3 Vlocal         = rotatePoint( qRotationToZ, V );
    const float3 Nlocal = float3( 0.0f, 0.0f, 1.0f );

    vec3 rayDirectionLocal = vec3( 0.0f, 0.0f, 0.0f );
    color3 sampleWeight;

    if ( m_MatType == DIFFUSE ) {
        BrdfData data;
        // Sample diffuse ray using cosine-weighted hemisphere sampling
        rayDirectionLocal =
            sampleHemisphereCook( vec2( uniform01( engine ), uniform01( engine ) ) );
        color3 albedo( 0.f );
        if ( m_texture != nullptr && ray->hasDifferentials ) {
            vec3 duv[2];
            duv[0] = vec3( intersection->dudx, intersection->dvdx, 0 );
            duv[1] = vec3( intersection->dudy, intersection->dvdy, 0 );
            albedo = m_texture->value( intersection->u, intersection->v, duv );
        }
        else { albedo = m_albedo; }
        if ( computeBrdfData( data,
                              Vlocal,
                              rayDirectionLocal,
                              Nlocal,
                              vec2( intersection->u, intersection->v ),
                              albedo,
                              m_metalness,
                              m_roughness ) ) {
            // Function 'diffuseTerm' is predivided by PDF of sampling the cosine weighted
            // hemisphere
            sampleWeight = data.diffuseReflectance * InvPi;

#if COMBINE_BRDFS_WITH_FRESNEL
            // Sample a half-vector of specular BRDF. Note that we're reusing random variable
            // 'u' here, but correctly it should be an new independent random number
            vec2 u         = vec2( uniform01( engine ), uniform01( engine ) );
            vec3 Hspecular = sampleBeckmannWalter( Vlocal, vec2( data.alpha ), u );

            // Clamp HdotL to small value to prevent numerical instability. Assume that rays
            // incident from below the hemisphere have been filtered
            float VdotH = max( 0.00001f, min( 1.0f, dot( Vlocal, Hspecular ) ) );
            sampleWeight *=
                ( vec3( 1.0f, 1.0f, 1.0f ) -
                  evalFresnel( data.specularF0, shadowedF90( data.specularF0 ), VdotH ) );
#endif
        }
    }
    else if ( m_MatType == SPECULAR ) {
        vec2 u      = vec2( uniform01( engine ), uniform01( engine ) );
        float alpha = m_roughness * m_roughness;
        color3 albedo( 0.f );
        if ( m_texture != nullptr && ray->hasDifferentials ) {
            vec3 duv[2];
            duv[0] = vec3( intersection->dudx, intersection->dvdx, 0 );
            duv[1] = vec3( intersection->dudy, intersection->dvdy, 0 );
            albedo = m_texture->value( intersection->u, intersection->v, duv );
        }
        else { albedo = m_albedo; }
        auto specularF0 = baseColorToSpecularF0( albedo, m_metalness );
        rayDirectionLocal =
            sampleSpecularMicrofacet( Vlocal, alpha, alpha * alpha, specularF0, u, sampleWeight );
    }

    // Prevent tracing direction with no contribution
    if ( luminance( sampleWeight ) == 0.0f ) return ret;

    // Transform sampled direction Llocal back to V vector space
    auto rayDirection =
        normalize( rotatePoint( invertRotation( qRotationToZ ), rayDirectionLocal ) );

    // Prevent tracing direction "under" the hemisphere (behind the triangle)
    if ( dot( intersection->normal, rayDirection ) <= 0.0f ) return ret;

    Ray ray_ref              = Ray( intersection->position + ( rayDirection * acne_eps ),
                       normalize( rayDirection ),
                       0,
                       10000,
                       ray->depth + 1 );
    ray_ref.hasDifferentials = true;
    Intersection temp_intersection;
    auto reflColor = trace_ray( scene, &ray_ref, tree, &temp_intersection, true );
    ret += reflColor * sampleWeight;
    return ret;
}

color3 Material::scatterColor( Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection ) {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        auto ret = color3( 0.0 );
        if ( m_MatType == TRANSPARENT ) {
            // REFRACTION + REFLECTION
            vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;

            float LdotH = dot( -ray->dir, normal );
            float f     = intersection->isOutside ? RDM_Fresnel( LdotH, 1.f, m_IOR )
                                                  : RDM_Fresnel( LdotH, m_IOR, 1.f );

            Ray new_ray;
            if ( f < 1.0f ) {
                if ( uniform01( engine ) < f ) { // reflect
                    vec3 r  = reflect( ray->dir, normal );
                    new_ray = Ray(
                        intersection->position + ( acne_eps * r ), r, 0, 100000, ray->depth + 1 );
                }
                else { // refract
                    float refractionRatio = intersection->isOutside ? ( 1.f / m_IOR ) : m_IOR;
                    vec3 refr             = refract( ray->dir, normal, refractionRatio );
                    new_ray               = Ray( intersection->position + ( acne_eps * refr ),
                                   refr,
                                   0,
                                   100000,
                                   ray->depth + 1 );
                }
            }
            else {
                vec3 r = reflect( ray->dir, normal );
                new_ray =
                    Ray( intersection->position + ( acne_eps * r ), r, 0, 100000, ray->depth + 1 );
            }
            Intersection new_inter;
            color3 color = trace_ray( scene, &new_ray, tree, &new_inter );
            ret += color;
        }
        else {
            ret += myScatter( ray, scene, tree, intersection );
            // ret += scratchAPixelScatter(ray, scene, tree, intersection);
        }

        return ret;
    }
    else if ( m_MatModel == MaterialModel::BLINN ) {

        auto ret = color3( 0.0 );

        if ( !isBlack( m_reflection ) ) {
            vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
            vec3 dirR   = reflect( ray->dir, normal );
            Ray new_ray;
            new_ray.dox = vec3( 0 );
            new_ray.doy = vec3( 0 );
            new_ray.ddx = vec3( 0 );
            new_ray.ddy = vec3( 0 );
            new_ray     = Ray( intersection->position + ( acne_eps * normal ),
                           normalize( dirR ),
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
                    r = r0 + ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) *
                                 ( 1 - c1 );
                else
                    r = r0 + ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) *
                                 ( 1 - c2 );

                if ( uniform01( engine ) < r ) {
                    new_ray = Ray( intersection->position + ( 0.01f * reflDir ),
                                   normalize( reflDir ),
                                   0,
                                   1000000,
                                   ray->depth + 1 );
                }
                else {
                    new_ray = Ray( intersection->position - ( 0.001f * normal ),
                                   normalize( refractDir ),
                                   0,
                                   100000,
                                   ray->depth + 1 );
                }
            }
            else {
                new_ray = Ray( intersection->position + ( 0.001f * normal ),
                               normalize( reflDir ),
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
            vec3 origin = intersection->position + ( acne_eps * wi );
            Ray rayS    = Ray( origin, wi, 0.f, sqrt( d_light2 ) );
            rayS.shadow = true;
            rayS.dox    = vec3( 0.f );
            rayS.doy    = vec3( 0.f );
            rayS.ddx    = vec3( 0.f );
            rayS.ddy    = vec3( 0.f );
            if ( intersectKdTree( scene, tree, &rayS, &temp_inter ) &&
                 rayS.tmax * rayS.tmax <
                     d_light2 * 0.99 ) { //|| dot(intersection->normal, wi) < 0){
                ret += color3( 0 );
            }
            else {
                color3 brdf = scene->GetMaterial( intersection->materialIndex ).m_albedo / Pi +
                              m_specularColor *
                                  Phong_BRDF( wi, ray->dir, intersection->normal, phongExponent );
                float J   = 1. * dot( Np, -wi ) / d_light2;
                float pdf = dot( axePO, dirA ) /
                            ( Pi * sphereL->geom.sphere.radius * sphereL->geom.sphere.radius );
                if ( pdf > 0 )
                    ret += scene->GetMaterial( sphereL->m_MaterialIndex ).m_emission /
                           sqr( sphereL->geom.sphere.radius ) *
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

            Ray ray_ref = Ray( intersection->position + ( acne_eps * dirA ),
                               normalize( dirA ),
                               0,
                               100000,
                               ray->depth + 1 );
            ray_ref.dox = vec3( 0 );
            ray_ref.doy = vec3( 0 );
            ray_ref.ddx = vec3( 0 );
            ray_ref.ddy = vec3( 0 );

            float pdfPhong =
                ( phongExponent + 1.f ) / ( 2.f * Pi ) * pow( dot( R, dirA ), phongExponent );
            float pdf =
                proba * dot( intersection->normal, dirA ) / ( Pi ) + ( 1.f - proba ) * pdfPhong;
            if ( pdf <= 0 ) return ret;

            Intersection temp_intersection;
            auto reflColor = trace_ray( scene, &ray_ref, tree, &temp_intersection, false );

            if ( sample_diffuse )
                iColor += reflColor * scene->GetMaterial( intersection->materialIndex ).m_albedo *
                          ( dot( intersection->normal, dirA ) / ( Pi ) / pdf );
            else
                iColor +=
                    reflColor *
                    ( dot( intersection->normal, dirA ) *
                      Phong_BRDF( dirA, ray->dir, intersection->normal, phongExponent ) / pdf ) *
                    m_specularColor;

            ret += iColor;
        }

        return ret;
    }
    return color3( 0.0 );
}

color3
Material::eval( Ray* ray, Intersection* intersection, const vec3& wi, float* scatteringPdf ) const {
    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        if ( m_MatType == TRANSPARENT ) return color3( 0 );
        BrdfData data;
        color3 bsdf( 0 );
        color3 albedo( 0 );
        if ( m_texture != nullptr && ray->hasDifferentials ) {
            vec3 duv[2];
            duv[0] = vec3( intersection->dudx, intersection->dvdx, 0 );
            duv[1] = vec3( intersection->dudy, intersection->dvdy, 0 );
            albedo = m_texture->value( intersection->u, intersection->v, duv );
        }
        else { albedo = m_albedo; }
        if ( computeBrdfData( data,
                              -ray->dir,
                              wi,
                              intersection->normal,
                              vec2( intersection->u, intersection->v ),
                              albedo,
                              m_metalness,
                              m_roughness ) ) {
            bsdf           = RDM_bsdf( data, nullptr, -1 );
            *scatteringPdf = sampleBeckmannWalterReflectionPdf(
                data.alpha, data.alphaSq, data.NdotH, data.VdotN, data.LdotH );
        }
        return bsdf;
    }
    else if ( m_MatModel == MaterialModel::BLINN ) { return color3( 0.0 ); }
    return color3( 0.0 );
}

color3 Material::sample( Ray* ray,
                         Intersection* intersection,
                         const vec2& uScattering,
                         vec3* wi,
                         float* scatteringPdf ) const {

    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) {
        if ( m_MatType == TRANSPARENT ) {
            // REFRACTION + REFLECTION
            vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;

            float LdotH = dot( -ray->dir, normal );
            float f     = intersection->isOutside ? RDM_Fresnel( LdotH, 1.f, m_IOR )
                                                  : RDM_Fresnel( LdotH, m_IOR, 1.f );

            Ray new_ray;
            if ( f < 1.0f ) {
                if ( uniform01( engine ) < f ) { // reflect
                    *wi = reflect( ray->dir, normal );
                }
                else { // refract
                    float refractionRatio = intersection->isOutside ? ( 1.f / m_IOR ) : m_IOR;
                    *wi                   = refract( ray->dir, normal, refractionRatio );
                }
            }
            else { *wi = reflect( ray->dir, normal ); }
            *scatteringPdf = 1.f;
            return color3( 1 );
        }
        vec3 V = -ray->dir;
        vec2 u( uniform01( engine ), uniform01( engine ) );
        color3 bsdf( 0 );
        // Ignore incident ray coming from "below" the hemisphere
        if ( !intersection->isOutside ) return color3( 0 );

        // Transform view direction into local space of our sampling routines
        // (local space is oriented so that its positive Z axis points along the shading normal)
        vec4 qRotationToZ = getRotationToZAxis( intersection->normal );
        vec3 Vlocal       = rotatePoint( qRotationToZ, V );

        float alpha = m_roughness;

        // Sample a microfacet normal (H) in local space
        vec3 Hlocal;
        if ( alpha == 0.0f ) {
            // Fast path for zero roughness (perfect reflection), also prevents NaNs appearing due
            // to divisions by zeroes
            Hlocal = vec3( 0.0f, 0.0f, 1.0f );
        }
        else {
            // For non-zero roughness, this calls VNDF sampling for GG-X distribution or Walter's
            // sampling for Beckmann distribution
            Hlocal = sampleBeckmannWalter( Vlocal, vec2( alpha, alpha ), u );
        }

        // Reflect view direction to obtain light vector
        vec3 Llocal = reflect( -Vlocal, Hlocal );
        // Transform sampled direction Llocal back to V vector space
        *wi = normalize( rotatePoint( invertRotation( qRotationToZ ), Llocal ) );

        BrdfData data;
        color3 albedo( 0.f );
        if ( m_texture != nullptr && ray->hasDifferentials ) {
            vec3 duv[2];
            duv[0] = vec3( intersection->dudx, intersection->dvdx, 0 );
            duv[1] = vec3( intersection->dudy, intersection->dvdy, 0 );
            albedo = m_texture->value( intersection->u, intersection->v, duv );
        }
        else { albedo = m_albedo; }
        if ( computeBrdfData( data,
                              -ray->dir,
                              *wi,
                              intersection->normal,
                              vec2( intersection->u, intersection->v ),
                              albedo,
                              m_metalness,
                              m_roughness ) ) {
            bsdf           = RDM_bsdf( data, nullptr, -1 );
            *scatteringPdf = sampleBeckmannWalterReflectionPdf(
                data.alpha, data.alphaSq, data.NdotH, data.VdotN, data.LdotH );
        }
        return bsdf;
    }
    else if ( m_MatModel == MaterialModel::BLINN ) { return color3( 0.0 ); }
    return color3( 0.0 );
}

color3 Material::refractionColor( Scene* scene,
                                  KdTree* tree,
                                  Ray* ray,
                                  Intersection* intersection,
                                  color3 reflectionShade,
                                  vec3 normal ) {

    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) { return color3( 0.0 ); }
    else if ( m_MatModel == MaterialModel::BLINN ) {
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
            Ray ray_refr = Ray( intersection->position + ( acne_eps * refractDir ),
                                refractDir,
                                0,
                                100000,
                                ray->depth + 1 );
            if ( ray->hasDifferentials ) {
                vec3 dndx = intersection->dn[0] * intersection->dudx +
                            intersection->dn[1] * intersection->dvdx;
                vec3 dndy = intersection->dn[0] * intersection->dudy +
                            intersection->dn[1] * intersection->dvdy;
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
                    r = r0 + ( 1.0 - r0 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) * ( 1 - c1 ) *
                                 ( 1 - c1 );
                else
                    r = r0 + ( 1.0 - r0 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) * ( 1 - c2 ) *
                                 ( 1 - c2 );
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
                // if ( scene->m_skyTexture != nullptr ) {
                //    vec3 dir = refractDir;
                //    float z  = asin( -dir.z ) / float( M_PI ) + 0.5;
                //    float x  = dir.x / ( abs( dir.x ) + abs( dir.y ) + 0.00001 );
                //    float y  = dir.y / ( abs( dir.x ) + abs( dir.y ) + 0.00001 );
                //    point3 p = point3( 0.5, 0.5, 0.0 ) +
                //               z * ( x * point3( 0.5, 0.5, 0.0 ) + y * point3( -0.5, 0.5, 0.0 ) );
                //    color3 env = 0.7f * scene->m_skyTexture->value( p.x, p.y );
                //    refractionShade += env;
                //}
            }
        }
        else { refractionShade += m_refraction * reflectionShade; }
        return refractionShade;
    }
    return color3( 0.0 );
}

color3 Material::reflectionColor( Scene* scene,
                                  KdTree* tree,
                                  Ray* ray,
                                  Intersection* intersection,
                                  color3& color,
                                  vec3 normal ) {

    if ( m_MatModel == MaterialModel::COOK_TORRANCE ) { return color3( 0.0 ); }
    else if ( m_MatModel == MaterialModel::BLINN ) {
        vec3 r = normalize( reflect( ray->dir, normal ) );
        Ray ray_ref =
            Ray( intersection->position + ( acne_eps * r ), r, 0, 100000, ray->depth + 1 );

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

        // if ( !temp_intersection.hit ) {
        //    if ( scene->m_skyTexture != nullptr ) {
        //        vec3 dir = r;
        //        float z  = asin( -dir.z ) / float( M_PI ) + 0.5;
        //        float x  = dir.x / ( abs( dir.x ) + abs( dir.y ) );
        //        float y  = dir.y / ( abs( dir.x ) + abs( dir.y ) );
        //        point3 p = point3( 0.5, 0.5, 0.0 ) +
        //                   z * ( x * point3( 0.5, 0.5, 0.0 ) + y * point3( -0.5, 0.5, 0.0 ) );
        //        // TODO: Multiply with intensity var
        //        color3 env = 0.7f * scene->m_skyTexture->value( p.x, p.y );
        //        color += m_reflection * env;
        //    }
        //}
        // else if ( intersection->isOutside ) {
        if ( intersection->isOutside ) { color += ( reflColor * m_reflection ); }
        return reflColor;
    }
    return color3( 0.0 );
}
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
    color3 f = scene->GetMaterial( intersection->materialIndex )
                   .sample_f( ray->dir,
                              &wi,
                              intersection->normal,
                              point2( 0, 0 ) /*sampler->Get2D()*/,
                              &pdf,
                              type );

    if ( pdf > 0 && !isBlack( f ) ) {
        vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
        Ray ray_ref =
            Ray( intersection->position + ( acne_eps * wi ), wi, 0, 100000, ray->depth + 1 );

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
    color3 f = scene->GetMaterial( intersection->materialIndex )
                   .sample_f( -ray->dir,
                              &wi,
                              intersection->normal,
                              point2( 0, 0 ) /*sampler->Get2D()*/,
                              &pdf,
                              type );

    if ( pdf > 0 && !isBlack( f ) ) {
        vec3 normal = intersection->isOutside ? intersection->normal : -intersection->normal;
        Ray ray_refr =
            Ray( intersection->position + ( acne_eps * wi ), wi, 0, 100000, ray->depth + 1 );
        if ( ray->hasDifferentials ) {
            vec3 dndx =
                intersection->dn[0] * intersection->dudx + intersection->dn[1] * intersection->dvdx;
            vec3 dndy =
                intersection->dn[0] * intersection->dudy + intersection->dn[1] * intersection->dvdy;
            float eta = 1.f / scene->GetMaterial( intersection->materialIndex ).m_IOR;
            if ( !intersection->isOutside ) {
                eta  = 1.f / scene->GetMaterial( intersection->materialIndex ).m_IOR;
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
