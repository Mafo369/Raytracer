#include "defines.h"
#include "raytracer.h"
#include "Light.h"
#include "sampling/sampling.h"
#include "Sky.h"

// from PBR
color3 directIllumination(Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection, Light* light) {
    color3 ret(0);

    vec3 wi;
    float lightPdf = 0, scatteringPdf = 0;
    bool vis;
    vec2 u(uniform01(engine), uniform01(engine));
    auto Li = light->sample_Li(scene, tree, *intersection, u, &wi, &lightPdf, &vis);
    if(lightPdf > 0 && !isBlack(Li)){
        color3 bsdf = intersection->mat->eval(ray, intersection, wi, &scatteringPdf) * max(dot(intersection->normal, wi), 0.f);
        if(!isBlack(bsdf)){
          if(vis){
            Li = color3(0);
          }
          if(!isBlack(Li)){
              float weight = PowerHeuristic(1, lightPdf, 1, scatteringPdf);
              ret += bsdf * Li * weight / lightPdf;
          }
        }
    }

    color3 bsdf = intersection->mat->sample(ray, intersection, &wi, &scatteringPdf) * max(dot(intersection->normal, wi), 0.f);
    if(!isBlack(bsdf) && scatteringPdf > 0){
        float weight = 1;
        lightPdf = light->pdf_Li(*intersection, wi);
        if(lightPdf == 0) return ret;
        weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
        Ray ray_ref;
        ray_ref.hasDifferentials = true;
        rayInit( &ray_ref,
                 intersection->position + ( wi * acne_eps ),
                 normalize( wi ),
                 ray->pixel,
                 0,
                 10000,
                 ray->depth + 1 );

        Intersection temp_intersection;
        auto foundIntersection = intersectKdTree(scene, tree, &ray_ref, &temp_intersection);
        color3 Li(0.f);
        if(foundIntersection){
          if(!isBlack(temp_intersection.mat->m_emission))
            Li = static_cast<ShapeLight*>(light)->L(temp_intersection, -wi);
        }
        else
          Li = scene->sky->getRadiance(*ray);
        if(!isBlack(Li)) ret += bsdf * Li * weight / scatteringPdf;
    }
    return ret;
}
