#ifndef __IMAGE_H__
#define __IMAGE_H__

#include "defines.h"

typedef struct image_s {
    size_t width;
    size_t height;
    color3 *data;
} Image;

color3 *getPixelPtr(Image *img, size_t x, size_t y);
Image *initImage(size_t width, size_t height);
void freeImage(Image *img);
void saveImage(Image *img, char *basename);


#endif
