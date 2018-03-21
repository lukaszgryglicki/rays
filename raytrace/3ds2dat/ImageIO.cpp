#include "ImageIO.h"
#include <stdio.h>

int init_bmp(BMPTag* b)
{
 int i;
 if (!b) return 0;
 b->ident[0]='B';
 b->ident[1]='M';
 b->fsize=0;
 b->dummy=0;
 b->offset=sizeof(BMPTag);
 b->bm_x=b->bm_y=0x20;
 b->dummy2=40;
 b->bpp=0x18;
 b->planes=1;
 b->compress=0;
 b->nbytes=4*32*32;
 for (i=0;i<4;i++) b->no_matter[i]=0;
 return 1;
}

bool LoadBMP(char* fn, unsigned char*& image, unsigned long& width, unsigned long& height)
{
 FILE* plik;
 char b,m;
 int i,j;
 BMPTag bm_handle;
 plik = fopen(fn, "rb");
 if (!plik) return false;
 if (!init_bmp(&bm_handle)) return false;
 i = fscanf(plik,"%c%c",&b,&m);
 if (i != 2) return 0;
 if (b != 'B' || m != 'M') return false;
 fread(&bm_handle.fsize,4,1,plik);
 fread(&bm_handle.dummy,4,1,plik);
 fread(&bm_handle.offset,4,1,plik);
 fread(&bm_handle.dummy2,4,1,plik);
 fread(&bm_handle.bm_x,4,1,plik);
 fread(&bm_handle.bm_y,4,1,plik);
 fread(&bm_handle.planes,2,1,plik);
 fread(&bm_handle.bpp,2,1,plik);
 if (bm_handle.bpp != 24) return false;
 fseek(plik,bm_handle.offset,SEEK_SET);
 width  = bm_handle.bm_x;
 height = bm_handle.bm_y;
 image = new unsigned char[3*width*height+1];
 //printf("%dx%d\n", width, height);
 fread(image, 1, 3*width*height, plik);
 /*for (i=0;i<height;i++)  
 for (j=0;j<width;j++)
    {
     fscanf(plik,"%c%c%c", &b,&g,&r);
     //set_color(t, j, i, r, g, b, 235);
	 image[3*(width * i + j)  ] = b;
	 image[3*(width * i + j)+1] = g;
	 image[3*(width * i + j)+2] = r;
	 image[3*(width * i + j)+3] = 255;
    }*/
 fclose(plik);
 printf("Image alignment is not supported yet.\n");
 return false;
}

bool LoadTGA(char* fn, unsigned char*& image, unsigned long& width, unsigned long& height)
{
 FILE* filePtr;
 TGATag tgaFile;
 unsigned char ucharBad;
 short int sintBad;
 long imageSize;
 int colorMode;
 long imageIdx;
 unsigned char colorSwap;
 filePtr = fopen(fn, "rb");
 if (!filePtr) return false;
 fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
 fread(&ucharBad, sizeof(unsigned char), 1, filePtr);

 fread(&tgaFile.imageTypeCode, sizeof(unsigned char), 1, filePtr);
 if (tgaFile.imageTypeCode != 2 && tgaFile.imageTypeCode != 3)
   {
	printf("TGA is probably compressed, and that's why it is not supported.\n");
    fclose(filePtr);
    return false;
   }
 fread(&sintBad, sizeof(short int), 1, filePtr);
 fread(&sintBad, sizeof(short int), 1, filePtr);
 fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
 fread(&sintBad, sizeof(short int), 1, filePtr);
 fread(&sintBad, sizeof(short int), 1, filePtr);

 fread(&tgaFile.imageWidth, sizeof(short int), 1, filePtr);
 fread(&tgaFile.imageHeight, sizeof(short int), 1, filePtr);

 fread(&tgaFile.bitCount, sizeof(unsigned char), 1, filePtr);
 fread(&ucharBad, sizeof(unsigned char), 1, filePtr);

 colorMode = tgaFile.bitCount / 8;
 imageSize = tgaFile.imageWidth * tgaFile.imageHeight * colorMode;

 tgaFile.imageData = new unsigned char[imageSize];
 fread(tgaFile.imageData, sizeof(unsigned char), imageSize, filePtr);

 for (imageIdx = 0; imageIdx < imageSize; imageIdx += colorMode)
   {
	   colorSwap = tgaFile.imageData[imageIdx];
	   tgaFile.imageData[imageIdx] = tgaFile.imageData[imageIdx + 2];
	   tgaFile.imageData[imageIdx + 2] = colorSwap;
   }
 fclose(filePtr);
 return true;
}

