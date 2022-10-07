#pragma once

#include "defines.h"

typedef struct image_s {
    size_t width;
    size_t height;
    color3 *data;
} RenderImage;

RenderImage *loadPng(const char *filename); //Not used

color3 *getPixelPtr(RenderImage *img, size_t x, size_t y);
RenderImage *initImage(size_t width, size_t height);
void freeImage(RenderImage *img);
void saveImage(RenderImage *img, char *basename);


