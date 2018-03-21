#include <math.h>
#include <stdlib.h>

#define REAL double
#define PI 3.14159265

int rec, xres, yres, nsteps, sidx;
double xr, yr, zr, xt, yt, zt, xs, ys, zs;
double xr0, yr0, zr0, xt0, yt0, zt0, xs0, ys0, zs0;
char fname[1024], outdir[128];

void rotatex(REAL** m, REAL ang)
{
 ang /= 180./PI;
 m[1][1] = cos(ang);
 m[2][1] = sin(ang);
 m[1][2] = -sin(ang);
 m[2][2] = cos(ang);
}


void rotatey(REAL** m, REAL ang)
{
 ang /= 180./PI;
 m[0][0] = cos(ang);
 m[2][0] = -sin(ang);
 m[0][2] = sin(ang);
 m[2][2] = cos(ang);
}


void rotatez(REAL** m, REAL ang)
{
 ang /= 180./PI;
 m[0][0] = cos(ang);
 m[1][0] = sin(ang);
 m[0][1] = -sin(ang);
 m[1][1] = cos(ang);
}


void translatex(REAL** m, REAL arg)
{
 m[0][3] = arg;
}


void translatey(REAL** m, REAL arg)
{
 m[1][3] = arg;
}


void translatez(REAL** m, REAL arg)
{
 m[2][3] = arg;
}


void translate(REAL** m, REAL x, REAL y, REAL z)
{
 translatex(m, x);
 translatey(m, y);
 translatez(m, z);
}


void scalex(REAL** m, REAL arg)
{
 m[0][0] = arg;
}


void scaley(REAL** m, REAL arg)
{
 m[1][1] = arg;
}


void scalez(REAL** m, REAL arg)
{
 m[2][2] = arg;
}


void scale(REAL** m, REAL x, REAL y, REAL z)
{
 scalex(m, x);
 scaley(m, y);
 scalez(m, z);
}

REAL* vector(int siz)
{
 if (siz <= 0) return NULL;
 return (REAL*)malloc(siz*sizeof(REAL));
}


REAL** matrix(int siz)
{
 REAL** mem;
 int i;
 if (siz <= 0) return NULL;
 mem = (REAL**)malloc(siz*sizeof(REAL*));
 for (i=0;i<siz;i++) mem[i] = vector(siz);
 return mem;
}


REAL** matrix2(int siz1, int siz2)
{
 REAL** mem;
 int i;
 if (siz1 <= 0 || siz2 <= 0) return NULL;
 mem = (REAL**)malloc(siz1*sizeof(REAL*));
 for (i=0;i<siz1;i++) mem[i] = vector(siz2);
 return mem;
}


REAL*** matrix3(int siz1, int siz2, int siz3)
{
 REAL*** mem;
 int i;
 if (siz1 <= 0 || siz2 <= 0 || siz3 <= 0) return NULL;
 mem = (REAL***)malloc(siz1*sizeof(REAL**));
 for (i=0;i<siz1;i++) mem[i] = matrix2(siz2, siz3);
 return mem;
}


void I_matrix(REAL** dst, int siz)
{
 int i,j;
 for (i=0;i<siz;i++) for (j=0;j<siz;j++)
   {
    if (i == j) dst[i][j] = 1.;
    else dst[i][j] = 0.;
   }
}


REAL** matrix_mul(REAL** m, REAL** n, int ma, int mb, int na, int nb)
{
 REAL** dst;
 int k,j,i;
 if (ma <= 0 || mb  <= 0 || na <= 0 || nb <=0 || mb != na) return NULL;
 if (!m || !n) return NULL;
 dst = (REAL**)malloc(ma*sizeof(REAL*));
 if (!dst) return NULL;
 for (i=0;i<ma;i++) dst[i] = (REAL*)malloc(nb*sizeof(REAL));
 for (i=0;i<ma;i++)
 for (j=0;j<nb;j++)
    {
     dst[i][j] = 0.0;
     for (k=0;k<mb;k++) dst[i][j] += m[i][k] * n[k][j];
    }
 return dst;
}


void matrix_mul_vector(REAL* dst, REAL** m, REAL* v, int len)
{
 int i,j;
 for (i=0;i<len;i++)
    {
     dst[i] = 0.0;
     for (j=0;j<len;j++) dst[i] += v[j] * m[i][j];
    }
}


void free_matrix(REAL** m, int siz)
{
 int i;
 if (siz <= 0) return;
 for (i=0;i<siz;i++) free(m[i]);
 free(m);
}


void free_matrix3(REAL*** m, int siz1, int siz2)
{
 int i,j;
 if (siz1 <= 0 || siz2 <= 0) return;
 for (i=0;i<siz1;i++)
 for (j=0;j<siz2;j++) free(m[i][j]);
 for (i=0;i<siz1;i++) free(m[i]);
 free(m);
}

void m_ownmatrix(REAL*** m, REAL** mat)
{
 REAL** newm;
 newm = matrix_mul(*m, mat, 4, 4, 4, 4);
 free_matrix(*m, 4);
 *m = newm;
}


void m_translate(REAL*** m, REAL x, REAL y, REAL z)
{
 REAL** trans;
 REAL** newm;
 trans = matrix(4);
 I_matrix(trans, 4);
 translate(trans, x, y, z);
 newm = matrix_mul(*m, trans, 4, 4, 4, 4);
 free_matrix(*m, 4);
 free_matrix(trans, 4);
 *m = newm;
}


void m_scale(REAL*** m, REAL x, REAL y, REAL z)
{
 REAL** scal;
 REAL** newm;
 scal = matrix(4);
 I_matrix(scal, 4);
 scale(scal, x, y, z);
 newm = matrix_mul(*m, scal, 4, 4, 4, 4);
 free_matrix(*m, 4);
 free_matrix(scal, 4);
 *m = newm;
}


void m_rotatex(REAL*** m, REAL ang)
{
 REAL** rot;
 REAL** newm;
 rot = matrix(4);
 I_matrix(rot, 4);
 rotatex(rot, ang);
 newm = matrix_mul(*m, rot, 4, 4, 4, 4);
 free_matrix(*m, 4);
 free_matrix(rot, 4);
 *m = newm;
}


void m_rotatey(REAL*** m, REAL ang)
{
 REAL** rot;
 REAL** newm;
 rot = matrix(4);
 I_matrix(rot, 4);
 rotatey(rot, ang);
 newm = matrix_mul(*m, rot, 4, 4, 4, 4);
 free_matrix(*m, 4);
 free_matrix(rot, 4);
 *m = newm;
}


void m_rotatez(REAL***m, REAL ang)
{
 REAL** rot;
 REAL** newm;
 rot = matrix(4);
 I_matrix(rot, 4);
 rotatez(rot, ang);
 newm = matrix_mul(*m, rot, 4, 4, 4, 4);
 free_matrix(*m, 4);
 free_matrix(rot, 4);
 *m = newm;
}

void runRT(REAL** m, int i)
{
    char cmd[2048];
    sprintf(cmd, "./rays.fast -u -b 1000000 -r %d -x %d -y %d -q 95 -J -i %s -o %s/%08d.bmp -5 \"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\"",
	 rec, xres, yres, fname, outdir, i, m[0][0], m[0][1], m[0][2], m[0][3], m[1][0], m[1][1], m[1][2], m[1][3], m[2][0], m[2][1], m[2][2], m[2][3], m[3][0], m[3][1], m[3][2], m[3][3]);
    printf("%s\n", cmd);
 system(cmd);
    sprintf(cmd, "rm -f %s/%08d.bmp", outdir, i);
 system(cmd);
}

void rotator()
{
 int i;
 double angle;
 double** m;
 double _xr, _yr, _zr, _xt, _yt, _zt, _xs, _ys, _zs;
 _xr = xr0;
 _yr = yr0;
 _zr = zr0;
 _xt = xt0;
 _yt = yt0;
 _zt = zt0;
 _xs = xs0;
 _ys = ys0;
 _zs = zs0;
 m = matrix(4);
 i = 0;
 for (i=0;i<nsteps;i++)
 {
  _xr += xr;
  _yr += yr;
  _zr += zr;
  _xt += xt;
  _yt += yt;
  _zt += zt;
  _xs *= xs;
  _ys *= ys;
  _zs *= zs;
  if (i >= sidx)
  {
	  I_matrix(m, 4);
	  m_translate(&m, _xt, _yt, _zt);
	  m_rotatez(&m, _zr);
	  m_rotatey(&m, _yr);
	  m_rotatex(&m, _xr);
	  m_scale(&m, _xs, _ys, _zs);
 	 runRT(m, i);
  }
 }
 free_matrix(m, 4);
}


int main(int lb, char** par)
{
 printf("args: filename rec xres yres xr yr zr xt yt zt xs ys zs nsteps outdir sidx xr0 yr0 zr0 xt0 yt0 zt0 xs0 ys0 zs0 \n");
 if (lb < 26) return 1;
 strcpy(fname, par[1]);
 sscanf(par[2], "%d", &rec);
 sscanf(par[3], "%d", &xres);
 sscanf(par[4], "%d", &yres);
 sscanf(par[5], "%lf", &xr);
 sscanf(par[6], "%lf", &yr);
 sscanf(par[7], "%lf", &zr);
 sscanf(par[8], "%lf", &xt);
 sscanf(par[9], "%lf", &yt);
 sscanf(par[10], "%lf", &zt);
 sscanf(par[11], "%lf", &xs);
 sscanf(par[12], "%lf", &ys);
 sscanf(par[13], "%lf", &zs);
 sscanf(par[14], "%d", &nsteps);
 strcpy(outdir, par[15]);
 sscanf(par[16], "%d", &sidx);
 sscanf(par[17], "%lf", &xr0);
 sscanf(par[18], "%lf", &yr0);
 sscanf(par[19], "%lf", &zr0);
 sscanf(par[20], "%lf", &xt0);
 sscanf(par[21], "%lf", &yt0);
 sscanf(par[22], "%lf", &zt0);
 sscanf(par[23], "%lf", &xs0);
 sscanf(par[24], "%lf", &ys0);
 sscanf(par[25], "%lf", &zs0);
 rotator();
 return 0;
}
