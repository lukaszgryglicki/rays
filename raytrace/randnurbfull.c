#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double d_rand()
{
 return (double)rand() / (double)RAND_MAX;
}

int i_rand(int i)
{
 return rand() % i;
}

void rnurbs(int num)
{
 int dimx,dimy;
 int nx, ny;
 int i,j,k;
 double x,y,z,tr,tg,tb;
 double sr,sg,sb,cr,cg,cb;
 double fr,fg,fb,sf,nd,tx,ty;
 double w,big,fact;
 printf("nNURBS: %d\n", num);
 for (i=0;i<num;i++)
  {
   printf("NURBS: %d\n",i);
   printf("  {\n");
   dimx = 2 + i_rand(4);
   dimy = 2 + i_rand(4);
   nx = dimx + 2 + i_rand(5);
   ny = dimy + 2 + i_rand(5);
   printf("DimU: %d\n", dimx);
   printf("DimV: %d\n", dimy);
   printf("nptsU: %d\n", nx);
   printf("nptsV: %d\n", ny);
   printf("divU: %d\n", 2*nx+i_rand(2*nx));
   printf("divV: %d\n", 2*ny+i_rand(2*nx));
   printf("ITex: %d\n", i_rand(2));
   printf("Faces: %d\n", i_rand(2)+1);
   printf("Tid: %d\n", i_rand(34));
   printf("userNodes: 0\n");
   printf("userKnots: 0\n");
   big = 600.;
   for (j=0;j<nx;j++)
     {
      printf("Line=%d\n", j);
      for (k=0;k<ny;k++)
        {
	 fact = 1.;
	 if (j == 0 || j == nx-1) fact *= .5;
	 if (j == 1 || j == nx-2) fact *= .8;
	 if (k == 0 || k == ny-1) fact *= .5;
	 if (k == 1 || k == ny-2) fact *= .8;
	 x = d_rand() * big - big/2.;
	 y = d_rand() * big - big/2.;
	 z = d_rand() * big - big/2.;
	 x *= fact;
	 y *= fact;
	 z *= fact;
	 tr = d_rand();
	 tg = d_rand();
	 tb = d_rand();
	 sr = d_rand();
	 sg = d_rand();
	 sb = d_rand();
	 cr = d_rand();
	 cg = d_rand();
	 cb = d_rand();
	 fr = 1. + d_rand()/3.;
	 fg = 1. + d_rand()/3.;
	 fb = 1. + d_rand()/3.;
	 sf = 100. + d_rand() * 200.;
	 nd = d_rand() * d_rand() * d_rand() * 6.;
	 tx = (double)j / (double)(nx-1);
	 ty = (double)k / (double)(ny-1);
	 tx += d_rand() * 0.2 - 0.1;
	 ty += d_rand() * 0.2 - 0.1;
	 w = d_rand() * 4.;
	 if (tx < 0.) tx  = 0.;
	 if (ty < 0.) ty  = 0.;
	 if (tx > 1.) tx  = 1.;
	 if (ty > 1.) ty  = 1.;
         printf("{v=(%f,%f,%f),t=(%f,%f,%f),s=(%f,%f,%f),c=(%f,%f,%f),"
         "f=(%f,%f,%f),sf=%f,nd=%f%%,tc=(%f,%f),w=%f}\n"
	 , x, y, z, tr, tg, tb, sr, sg, sb, cr, cg, cb, fr, fg, fb, 
	 sf, nd, tx, ty, w);
        }
     }
   printf("}\n");
  }
}

int main(int lb, char** par)
{
 time_t tm;
 time(&tm);
 srand((int)tm);
 if (lb != 2) { printf("%s numNurbs\n", par[0]); return 1; }
 else rnurbs(atoi(par[1]));
 return 0;
}
