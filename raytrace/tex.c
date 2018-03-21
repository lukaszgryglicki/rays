#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

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

void init_bmp(BMPTag* b)
{
 int i;
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
 b->nbytes=3*32*32;
 for (i=0;i<4;i++) b->no_matter[i]=0;
}

void error(char* fmt, ...)
{
 va_list lst;
 va_start(lst,fmt);
 printf("Error: \t");
 vprintf(fmt,lst);
 printf("\n\tClient halted\n");
 fflush(stdout);
 va_end(lst);
 exit(1);
}


void write_bmp(char* fn, float fff)
{
 FILE* plik;
 BMPTag bm_handle;
 int x,y,i,j,r,g,b;
 float dx,dy,d,angle;
 if (!(plik = fopen(fn, "w"))) error("cannot write to file: %s", fn);
 init_bmp(&bm_handle);
 fprintf(plik,"%c%c",'B', 'M');
 x = 512;
 y = 512;
 bm_handle.bm_y = x;
 bm_handle.bm_x = y;
 printf("Dimnesions: (%dx%d)\n", bm_handle.bm_x, bm_handle.bm_y);
 bm_handle.fsize = sizeof(BMPTag)+(bm_handle.bm_y*bm_handle.bm_x*3);
 fwrite(&bm_handle.fsize,4,1,plik);
 fwrite(&bm_handle.dummy,4,1,plik);
 bm_handle.offset=sizeof(BMPTag);
 bm_handle.planes=1;
 bm_handle.bpp=24;
 bm_handle.nbytes = x * y * 3;
 fwrite(&bm_handle.offset,4,1,plik);
 fwrite(&bm_handle.dummy2,4,1,plik);
 fwrite(&bm_handle.bm_x,4,1,plik);
 fwrite(&bm_handle.bm_y,4,1,plik);
 fwrite(&bm_handle.planes,2,1,plik);
 fwrite(&bm_handle.bpp,2,1,plik);
 fwrite(&bm_handle.compress,4,1,plik);
 fwrite(&bm_handle.nbytes,4,1,plik);
 for (i=0;i<4;i++)  fwrite(&bm_handle.no_matter[i],4,1,plik);
 fseek(plik,bm_handle.offset,SEEK_SET);
  for (i=0;i<bm_handle.bm_y;i++)  
    {
     for (j=0;j<bm_handle.bm_x;j++) 
       {
        dx = (float)(i - x/2.) / (float)(x/2.);
        dy = (float)(j - y/2.) / (float)(y/2.);
	
        d = sqrt(dx*dx + dy*dy);
	if (d >= 1.) d = 1 - (d - 1.);

	d = pow(d, fff);
	
        r = (int)(d * 255.);

	dx = pow(fabs(dx), 1. - 1./fff);
	dy = pow(fabs(dy), 1. - 1./fff);

	if (dx < 0.000001) dx = 0.000001;
	if (dy < 0.000001) dy = 0.000001;
	
	angle = atan(dy/dx);
	angle += M_PI / 2.;
	angle /= M_PI;

	angle = pow(angle, fff);
	   
	g = (int)(angle * 255.);
	
	angle = atan(dx/dy);
	angle += M_PI / 2.;
	angle /= M_PI;

	angle = pow(angle, fff);
	   
	b = (int)(angle * 255.);

        fprintf(plik,"%c%c%c", b, g, r);
       }
    }
 fclose(plik);
 printf("Writen bitmap: %s\n", fn);
}

int main(int lb, char** par)
{
 if (lb < 3) { printf("want 2 args: texname and pow_factor\n"); return 1; }
 write_bmp(par[1], atof(par[2]));
 return 0;
}

