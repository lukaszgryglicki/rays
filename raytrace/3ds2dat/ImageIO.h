#ifndef _IMAGE_IO_H__MDBMA_
#define _IMAGE_IO_H__MDBMA_

typedef struct _BMPTag
{
 char ident[2];
 int fsize;
 int dummy;
 int offset;
 int dummy2;
 int bm_x;
 int bm_y;
 short planes;
 short bpp;
 int compress;
 int nbytes;
 int no_matter[4];
} BMPTag;

typedef struct _TGATag
{
 unsigned char imageTypeCode;
 short int imageWidth;
 short int imageHeight;
 unsigned char bitCount;
 unsigned char* imageData;
} TGATag;

bool LoadBMP(char* fn, unsigned char*& image, unsigned long& width, unsigned long& height);
bool LoadTGA(char* fn, unsigned char*& image, unsigned long& width, unsigned long& height);

#endif

