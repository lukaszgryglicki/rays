/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define COORD_LX 0
#define COORD_LY 1
#define COORD_LZ 2
#define COORD_TR 3
#define COORD_TG 4
#define COORD_TB 5
#define COORD_SR 6
#define COORD_SG 7
#define COORD_SB 8
#define COORD_CR 9
#define COORD_CG 10
#define COORD_CB 11
#define COORD_FR 12
#define COORD_FG 13
#define COORD_FB 14
#define COORD_SF 15
#define COORD_ND 16
#define COORD_TX 17
#define COORD_TY 18
#define S 19

typedef struct _NURBS
{
 int n1,n2;
 int p1,p2;
 int m1,m2;
 int d1,d2;
 int ntri;
 int idx;
 int itex;
 int faces;
 int tid;
 double*** D;
 double*** normals;
 double*** values;
 double** w;
 double* t1;
 double* t2;
 double* knot1;
 double* knot2;
} NURBS;

int lt;

void panic(char* str)
{
 printf("\npanic: %s\n", str);
 exit(1);
}

double basis_bspline(double* knot, double val, int dim, int off, int I, int J)
{
 register double up;
 register double down;
 double factor1;
 double factor2;
 register volatile int Ioff, Joff;
 factor1 = factor2 = 0.0;
 Ioff = I+off;
 Joff = J+off;
 if (dim == 0)
   {
    if (val >= knot[Ioff] && val <= knot[Joff]) return 1.0;
    else return 0.0;
   }
 up = (val-knot[Ioff]) * basis_bspline(knot, val, dim-1, off, I, J-1);
 down = knot[Joff-1] - knot[Ioff];
 if (down==0.0) factor1 = 0.0;
 else factor1 = up/down;
 up = (knot[Joff]-val) * basis_bspline(knot, val, dim-1, off, I+1, J);
 down = knot[Joff] - knot[Ioff+1];
 if (down==0.0) factor2 += 0.0;
 else factor2 = up/down;
 return factor1+factor2;
}

double bsp(double* knot, double t, int dim, int i)
{
 return basis_bspline(knot, t, dim, dim, i-dim, i+1);
}

double nurbs(NURBS* nurb, double u, double v, int xyz)
{
 register int i,j;
 double up,down;
 double factor;
 double b1,b2;
 down = 0.0;
 up = 0.0;
 for (i=0;i<nurb->n1;i++)
 {
 b1 = bsp(nurb->knot1, u, nurb->p1, i);
 for (j=0;j<nurb->n2;j++)
   {
    b2 = bsp(nurb->knot2, v, nurb->p2, j);
    factor = nurb->D[i][j][xyz]*b1*b2*nurb->w[i][j];
    up += factor;
    down += b1*b2*nurb->w[i][j];
   }
 }
 if (down==0.0) return 0.0;
 return up/down;
}


double* vector(int siz)
{
 if (siz <= 0) return NULL;
 return (double*)malloc(siz*sizeof(double));
}


double** matrix2(int siz1, int siz2)
{
 double** mem;
 int i;
 if (siz1 <= 0 || siz2 <= 0) return NULL;
 mem = (double**)malloc(siz1*sizeof(double*));
 for (i=0;i<siz1;i++) mem[i] = vector(siz2);
 return mem;
}


double*** matrix3(int siz1, int siz2, int siz3)
{
 double*** mem;
 int i;
 if (siz1 <= 0 || siz2 <= 0 || siz3 <= 0) return NULL;
 mem = (double***)malloc(siz1*sizeof(double**));
 for (i=0;i<siz1;i++) mem[i] = matrix2(siz2, siz3);
 return mem;
}

void read_nurbs(FILE* f, NURBS* n)
{
 int idx,nr;
 int uK, uN;
 int i,j,k;
 nr = fscanf(f, "NURBS: %d\n", &idx);
 if (nr != 1 || idx < 0) panic("cannot read nurbs header");
 fscanf(f, "{\n");
 nr = fscanf(f, "DimU: %d\n", &n->p1);
 if (nr != 1) panic("cannot read DimU");
 nr = fscanf(f, "DimV: %d\n", &n->p2);
 if (nr != 1) panic("cannot read DimV");
 nr = fscanf(f, "nptsU: %d\n", &n->n1);
 if (nr != 1) panic("cannot read nptsU");
 nr = fscanf(f, "nptsV: %d\n", &n->n2);
 if (nr != 1) panic("cannot read nptsV");
 nr = fscanf(f, "divU: %d\n", &n->d1);
 if (nr != 1) panic("cannot read divU");
 nr = fscanf(f, "divV: %d\n", &n->d2);
 nr = fscanf(f, "ITex: %d\n", &n->itex);
 if (nr != 1) panic("cannot read ITex");
 nr = fscanf(f, "Faces: %d\n", &n->faces);
 if (nr != 1) panic("cannot read Faces");
 nr = fscanf(f, "Tid: %d\n", &n->tid);
 if (nr != 1) panic("cannot read Tid");
 if (nr != 1) panic("cannot read divV");
 if (n->n1 <= n->p1) panic("not enough U points");
 if (n->n2 <= n->p2) panic("not enough V points");
 if (n->n1 < 2 || n->n2 < 2) panic("not enough points");
 if (n->p1 < 1 || n->p2 < 1) panic("too low degree");
 if (n->p1 >= 10 || n->p2 >= 10) panic("too high degree");
 if (n->n1 > 32 || n->n2 > 32) panic("too high points number");
 nr = fscanf(f, "userNodes: %d\n", &uN);
 if (nr != 1) panic("cannot read uN");
 nr = fscanf(f, "userKnots: %d\n", &uK);
 if (nr != 1) panic("cannot read uK");
 n->t1 = (double*)malloc(n->n1*sizeof(double));
 n->t2 = (double*)malloc(n->n2*sizeof(double));
 if (uN)
   {
    fscanf(f, "UNodes:");
    for (i=0;i<n->n1;i++) 
      {
       nr = fscanf(f, " %lf", &n->t1[i]);
       if (nr != 1) panic("cannot read Unode");
      }
    fscanf(f, "\n");
    fscanf(f, "VNodes:");
    for (i=0;i<n->n2;i++) 
      {
       nr = fscanf(f, " %lf", &n->t2[i]);
       if (nr != 1) panic("cannot read Vnode");
      }
    fscanf(f, "\n");
   }
 else
   {
    for (i=0;i<n->n1;i++) n->t1[i] = (double)i / (double)(n->n1 - 1);
    for (i=0;i<n->n2;i++) n->t2[i] = (double)i / (double)(n->n2 - 1);
   }
 /*printf("NodesU:");
 for (i=0;i<n->n1;i++) printf(" %f", n->t1[i]);
 printf("\n");
 printf("NodesV:");
 for (i=0;i<n->n2;i++) printf(" %f", n->t2[i]);
 printf("\n");*/
 n->m1 = n->n1 + n->p1;
 n->m2 = n->n2 + n->p2;
 n->knot1 = (double*)malloc((n->m1+1)*sizeof(double));
 n->knot2 = (double*)malloc((n->m2+1)*sizeof(double));
 if (uK)
   {
    fscanf(f, "UKnots:");
    for (i=0;i<=n->m1;i++) 
      {
       nr = fscanf(f, " %lf", &n->knot1[i]);
       if (nr != 1) panic("cannot read Uknot");
      }
    fscanf(f, "\n");
    fscanf(f, "VKnots:");
    for (i=0;i<=n->m2;i++) 
      {
       nr = fscanf(f, " %lf", &n->knot2[i]);
       if (nr != 1) panic("cannot read Vknot");
      }
    fscanf(f, "\n");
   }
 else
   {
    for (i=0;i<=n->p1;i++) n->knot1[i] = 0.;
    for (i=0;i<=n->p2;i++) n->knot2[i] = 0.;
    for (j=1;j<n->n1-n->p1;j++)
      {
	n->knot1[j+n->p1] = 0.;
	for (i=j;i<j+n->p1;i++) n->knot1[j+n->p1] += n->t1[i];
	n->knot1[j+n->p1] /= (double)n->p1;
      }
    for (j=1;j<n->n2-n->p2;j++)
      {
	n->knot2[j+n->p2] = 0.;
	for (i=j;i<j+n->p2;i++) n->knot2[j+n->p2] += n->t2[i];
	n->knot2[j+n->p2] /= (double)n->p2;
      }
    for (i=n->m1-n->p1;i<=n->m1;i++) n->knot1[i] = 1.;
    for (i=n->m2-n->p2;i<=n->m2;i++) n->knot2[i] = 1.;
   }
 /*printf("KnotsU:");
 for (i=0;i<=n->m1;i++) printf(" %f", n->knot1[i]);
 printf("\n");
 printf("KnotsV:");
 for (i=0;i<=n->m2;i++) printf(" %f", n->knot2[i]);
 printf("\n");*/
 n->D = matrix3(n->n1, n->n2, S);
 n->w = matrix2(n->n1, n->n2);
 for (i=0;i<n->n1;i++)
   {
    nr = fscanf(f, "Line=%d\n", &k);
/*    printf(" i = %d\n", i);*/
    if (nr != 1 || k != i) panic("No line definition");
    for (j=0;j<n->n2;j++)
      {
       nr = fscanf(f, 
       "{v=(%lf,%lf,%lf),t=(%lf,%lf,%lf),s=(%lf,%lf,%lf),c=(%lf,%lf,%lf),"
       "f=(%lf,%lf,%lf),sf=%lf,nd=%lf%%,tc=(%lf,%lf),w=%lf}\n",&n->D[i][j][0],&n->D[i][j][1],&n->D[i][j][2],
        &n->D[i][j][3],&n->D[i][j][4],&n->D[i][j][5],&n->D[i][j][6],&n->D[i][j][7],&n->D[i][j][8],
        &n->D[i][j][9],&n->D[i][j][10],&n->D[i][j][11],&n->D[i][j][12],&n->D[i][j][13],&n->D[i][j][14],
	&n->D[i][j][15],&n->D[i][j][16],&n->D[i][j][17],&n->D[i][j][18], &n->w[i][j]);
        if (nr != S+1) { printf("data scanf error at: (%d,%d)\n", i,j); panic("error"); }
      }
     fscanf(f, "\n");
    }
 n->ntri = 2*n->d1*n->d2;
 fscanf(f, "}\n");
}

double b0(double* knot, double t, int i)
{
/* printf("checking: [%d,%d] %f in [%f,%f]\n", i,i+1,t,knot[i],knot[i+1]);*/
 if (t >= knot[i] && t <= knot[i+1]) return 1.;
 else return 0;
}

double** pre_basis(double* knot, double t, int dim, int idx)
{
 int i,crd;
 double** prep;
 double up;
 double down;
 double f1, f2;
/* printf("t = %f, dim = %d, idx = %d\n", t, dim, idx);*/
 prep = (double**)malloc((dim+1)*sizeof(double*));
 for (i=0;i<=dim;i++) prep[i] = (double*)malloc(((dim-i)+1)*sizeof(double));
 for (crd=0;crd<=dim;crd++)
  {
   if (crd == 0)
     {
      for (i=idx;i<=idx+dim;i++)
        {
         prep[crd][i-idx] = b0(knot, t, i);
/*         printf("prep[%d][%d] = %f\n", crd, i-idx, prep[crd][i-idx]);*/
        }
     }
   else
     {
      for (i=idx;i<=idx+dim-crd;i++)
        {
         up = (t - knot[i]) * prep[crd-1][i-idx];
	 down = knot[i+crd] - knot[i];
	 if (down != 0.) f1 = up/down;
	 else f1 = 0.;
	 up = (knot[i+crd+1] - t) * prep[crd-1][i-idx+1];
	 down = knot[i+crd+1] - knot[i+1];
	 if (down != 0.) f2 = up/down;
	 else f2 = 0.;
	 prep[crd][i-idx] = f1 + f2;
/*         printf("prep[%d][%d] = %f\n", crd, i-idx, prep[crd][i-idx]);*/
        }
     }
  }
 return prep;
}

double*** precompute_nurbs(double* knot, double t, int dim, int npts)
{
 int i,crd,idx;
 double*** prep;
 double up;
 double down;
 double f1, f2;
/* printf("t = %f, dim = %d, idx = %d\n", t, dim, idx);*/
 prep = (double***)malloc(npts*sizeof(double**));
 for (idx=0;idx<npts;idx++)
 {
  prep[idx] = (double**)malloc((dim+1)*sizeof(double*));
  for (i=0;i<=dim;i++) prep[idx][i] = (double*)malloc(((dim-i)+1)*sizeof(double));
  for (crd=0;crd<=dim;crd++)
   {
    if (crd == 0)
      {
       for (i=idx;i<=idx+dim;i++)
         {
          prep[idx][crd][i-idx] = b0(knot, t, i);
         }
      }
    else
      {
       for (i=idx;i<=idx+dim-crd;i++)
         {
          up = (t - knot[i]) * prep[idx][crd-1][i-idx];
	  down = knot[i+crd] - knot[i];
	  if (down != 0.) f1 = up/down;
	  else f1 = 0.;
	  up = (knot[i+crd+1] - t) * prep[idx][crd-1][i-idx+1];
	  down = knot[i+crd+1] - knot[i+1];
	  if (down != 0.) f2 = up/down;
	  else f2 = 0.;
	  prep[idx][crd][i-idx] = f1 + f2;
        }
     }
  }
 }
 return prep;
}

void prenurbs_free(double**** ptr, int dim, int npts)
{
 int i,j;
 for (i=0;i<npts;i++)
   for (j=0;j<=dim;j++) free((*ptr)[i][j]);
 for (i=0;i<npts;i++) free((*ptr)[i]);
 free(*ptr);
 *ptr = NULL;
}

double fastnurbs(NURBS* nurb, double u, double v, int xyz)
{
 register int i,j;
 double up,down;
 double factor;
 double b1,b2;
 double*** pre1;
 double*** pre2;
 pre1 = precompute_nurbs(nurb->knot1, u, nurb->p1, nurb->n1);
 pre2 = precompute_nurbs(nurb->knot2, v, nurb->p2, nurb->n2);
 down = 0.0;
 up = 0.0;
 for (i=0;i<nurb->n1;i++)
   {
    b1 = pre1[i][nurb->p1][0];
    for (j=0;j<nurb->n2;j++)
       {
        b2 = pre2[j][nurb->p2][0];
        factor = nurb->D[i][j][xyz]*b1*b2*nurb->w[i][j];
        up += factor;
        down += b1*b2*nurb->w[i][j];
       }
    }
 prenurbs_free(&pre1, nurb->p1, nurb->n1);
 prenurbs_free(&pre2, nurb->p2, nurb->n2);
 if (down != 0.) return up/down;
 else return 0.;
}

void fastnurbs_array(NURBS* nurb, double u, double v, double* tab, int siz)
{
 register int i,j,k;
 double up,down;
 double factor;
 double b1,b2;
 double*** pre1;
 double*** pre2;
 pre1 = precompute_nurbs(nurb->knot1, u, nurb->p1, nurb->n1);
 pre2 = precompute_nurbs(nurb->knot2, v, nurb->p2, nurb->n2);
 for (k=0;k<siz;k++)
 {
  down = 0.0;
  up = 0.0;
   for (i=0;i<nurb->n1;i++)
   {
    b1 = pre1[i][nurb->p1][0];
    for (j=0;j<nurb->n2;j++)
       {
        b2 = pre2[j][nurb->p2][0];
        factor = nurb->D[i][j][k]*b1*b2*nurb->w[i][j];
        up += factor;
        down += b1*b2*nurb->w[i][j];
       }
    }
  if (down != 0.) tab[k] = up/down; 
  else tab[k] = 0.0;
/*  printf("tab[%d] = %f\n", k, tab[k]);*/
 }
 prenurbs_free(&pre1, nurb->p1, nurb->n1);
 prenurbs_free(&pre2, nurb->p2, nurb->n2);
}

void calc_normal(NURBS* nu, double u, double v, double* x, double* y, double* z)
{
 double x1,x2,x3;
 double y1,y2,y3;
 double z1,z2,z3;
 double nx,ny,nz;
 double dx1,dx2;
 double dy1,dy2;
 double dz1,dz2;
 double epsu;
 double epsv;
 double len;
/* double** pre;*/
 
 epsu = 3e-4;
 epsv = 3e-4;

 if (1 - 2.*epsu <= u) epsu *= -1.;
 if (1 - 2.*epsv <= v) epsv *= -1.;

 
/* pre = pre_basis(nu->knot1, 0.345, 3, 1);*/
/* len = bsp(nu->knot1, 0.345, 3, 1);*/
/* printf("Should equal: (%f == %f)\n", len, pre[3][0]);*/
 
 
 
 x1 = fastnurbs(nu, u, v, COORD_LX);
 x2 = fastnurbs(nu, u+epsu, v, COORD_LX);
 x3 = fastnurbs(nu, u, v+epsv, COORD_LX);
 y1 = fastnurbs(nu, u, v, COORD_LY);
 y2 = fastnurbs(nu, u+epsu, v, COORD_LY);
 y3 = fastnurbs(nu, u, v+epsv, COORD_LY);
 z1 = fastnurbs(nu, u, v, COORD_LZ);
 z2 = fastnurbs(nu, u+epsu, v, COORD_LZ);
 z3 = fastnurbs(nu, u, v+epsv, COORD_LZ);
 
 dx1 = x2-x1;
 dx2 = x3-x1;
 dy1 = y2-y1;
 dy2 = y3-y1;
 dz1 = z2-z1;
 dz2 = z3-z1;
 
 nx = dy1*dz2 - dz1*dy2;
 ny = dz1*dx2 - dx1*dz2;
 nz = dx1*dy2 - dy1*dx2;
 
 len = sqrt(nx*nx + ny*ny + nz*nz);
 nx /= len;
 ny /= len;
 nz /= len;

 if (epsu * epsv > 0.) 
   {
    nx *= -1.;
    ny *= -1.;
    nz *= -1.;
   }

 *x = nx;
 *y = ny;
 *z = nz;
}

void write_nurbs(FILE* f, NURBS* n, int idx)
{
 double u,U,v,V; 
/* double*** vals;*/
 int i,j, k;
 n->values = matrix3(n->d1+1, n->d2+1, S);
/* vals     = matrix3(n->d1+1, n->d2+1, S);*/
 n->normals = matrix3(n->d1+1, n->d2+1, 3);
 for (i=0;i<=n->d1;i++)
   {
    u = (double)i/(double)n->d1;
    for (j=0;j<=n->d2;j++)
     {
      v = (double)j/(double)n->d2;
      calc_normal(n, u, v, &n->normals[i][j][0], &n->normals[i][j][1], &n->normals[i][j][2]);
/*      for (k=0;k<S;k++) vals[i][j][k] = nurbs(n, u, v, k);*/
      fastnurbs_array(n, u, v, n->values[i][j], S);
      /*for (k=0;k<S;k++) if (vals[i][j][k] != n->values[i][j][k]) 
        { 
	 printf("error! i,j,k=%d,%d,%d (%f != %f)\n", i,j,k,vals[i][j][k], n->values[i][j][k]); 
	}*/
/*      printf("(%f,%f,%f)\n", n->values[i][j][0], n->values[i][j][1], n->values[i][j][2]);*/
     }
   }
 k = idx + n->idx;
 for (i=0;i<n->d1;i++)
  {
   u = (double)i/(double)n->d1;
   U = (double)(i+1)/(double)n->d1;
   for (j=0;j<n->d2;j++)
    {
     v = (double)j/(double)n->d2;
     V = (double)(j+1)/(double)n->d2;
     
       fprintf(f,"Triangle: %d\n", k);
       fprintf(f,"{\n");
       fprintf(f,"a: Vertex: (%f,%f,%f)\n",
	 n->values[i][j][COORD_LX], n->values[i][j][COORD_LY], n->values[i][j][COORD_LZ]);
       fprintf(f,"b: Vertex: (%f,%f,%f)\n",
	 n->values[i+1][j+1][COORD_LX], n->values[i+1][j+1][COORD_LY], n->values[i+1][j+1][COORD_LZ]);
       fprintf(f,"c: Vertex: (%f,%f,%f)\n",
	 n->values[i+1][j][COORD_LX], n->values[i+1][j][COORD_LY], n->values[i+1][j][COORD_LZ]);
       if (!n->itex)
         {
          fprintf(f,"texA: TexCoord: (%f,%f)\n", u, v);
          fprintf(f,"texB: TexCoord: (%f,%f)\n", U, V);
          fprintf(f,"texC: TexCoord: (%f,%f)\n", U, v);
         }
       else
         {
          fprintf(f,"texA: TexCoord: (%f,%f)\n", 
		  n->values[i][j][COORD_TX], n->values[i][j][COORD_TY]);
          fprintf(f,"texB: TexCoord: (%f,%f)\n", 
		  n->values[i+1][j+1][COORD_TX], n->values[i+1][j+1][COORD_TY]);
          fprintf(f,"texC: TexCoord: (%f,%f)\n", 
		  n->values[i+1][j][COORD_TX], n->values[i+1][j][COORD_TY]);
         }
       fprintf(f,"na: Vector: (%f,%f,%f)\n",
	       n->normals[i][j][COORD_LX], n->normals[i][j][COORD_LY], n->normals[i][j][COORD_LZ]);
       fprintf(f,"nb: Vector: (%f,%f,%f)\n",
	       n->normals[i+1][j+1][COORD_LX], n->normals[i+1][j+1][COORD_LY], n->normals[i+1][j+1][COORD_LZ]);
       fprintf(f,"nc: Vector: (%f,%f,%f)\n",
	       n->normals[i+1][j][COORD_LX], n->normals[i+1][j][COORD_LY], n->normals[i+1][j][COORD_LZ]);
       fprintf(f,"transparencyA: RGB: (%f,%f,%f)\n", 
	       n->values[i][j][COORD_TR], n->values[i][j][COORD_TG], n->values[i][j][COORD_TB]);
       fprintf(f,"specularA: RGB: (%f,%f,%f)\n", 
	       n->values[i][j][COORD_SR], n->values[i][j][COORD_SG], n->values[i][j][COORD_SB]);
       fprintf(f,"diffuseA: RGB: (%f,%f,%f)\n", 
	       n->values[i][j][COORD_CR], n->values[i][j][COORD_CG], n->values[i][j][COORD_CB]);

       fprintf(f,"transparencyB: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j+1][COORD_TR], n->values[i+1][j+1][COORD_TG], n->values[i+1][j+1][COORD_TB]);
       fprintf(f,"specularB: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j+1][COORD_SR], n->values[i+1][j+1][COORD_SG], n->values[i+1][j+1][COORD_SB]);
       fprintf(f,"diffuseB: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j+1][COORD_CR], n->values[i+1][j+1][COORD_CG], n->values[i+1][j+1][COORD_CB]);
       
       fprintf(f,"transparencyC: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j][COORD_TR], n->values[i+1][j][COORD_TG], n->values[i+1][j][COORD_TB]);
       fprintf(f,"specularC: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j][COORD_SR], n->values[i+1][j][COORD_SG], n->values[i+1][j][COORD_SB]);
       fprintf(f,"diffuseC: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j][COORD_CR], n->values[i+1][j][COORD_CG], n->values[i+1][j][COORD_CB]);
       fprintf(f,"transparencyFactR: (1,%f)\n", n->values[i][j][COORD_FR]);
       fprintf(f,"transparencyFactG: (1,%f)\n", n->values[i][j][COORD_FG]);
       fprintf(f,"transparencyFactB: (1,%f)\n", n->values[i][j][COORD_FB]);
       fprintf(f,"normalDist: %f%%\n", n->values[i][j][COORD_ND]);
       fprintf(f,"specularFact: %f\n", n->values[i][j][COORD_SF]);
       fprintf(f,"faces: %d\n", n->faces);
       fprintf(f,"texture: %d\n", n->tid);
       fprintf(f,"}\n");
       k++;
       fprintf(f,"Triangle: %d\n", k);
       fprintf(f,"{\n");
       fprintf(f,"a: Vertex: (%f,%f,%f)\n",
	 n->values[i][j][COORD_LX], n->values[i][j][COORD_LY], n->values[i][j][COORD_LZ]);
       fprintf(f,"b: Vertex: (%f,%f,%f)\n",
	 n->values[i][j+1][COORD_LX], n->values[i][j+1][COORD_LY], n->values[i][j+1][COORD_LZ]);
       fprintf(f,"c: Vertex: (%f,%f,%f)\n",
	 n->values[i+1][j+1][COORD_LX], n->values[i+1][j+1][COORD_LY], n->values[i+1][j+1][COORD_LZ]);
       if (!n->itex)
         {
          fprintf(f,"texA: TexCoord: (%f,%f)\n", u, v);
          fprintf(f,"texB: TexCoord: (%f,%f)\n", u, V);
          fprintf(f,"texC: TexCoord: (%f,%f)\n", U, V);
	 }
       else
         {
          fprintf(f,"texA: TexCoord: (%f,%f)\n", 
		  n->values[i][j][COORD_TX], n->values[i][j][COORD_TY]);
          fprintf(f,"texB: TexCoord: (%f,%f)\n", 
		  n->values[i][j+1][COORD_TX], n->values[i][j+1][COORD_TY]);
          fprintf(f,"texC: TexCoord: (%f,%f)\n", 
		  n->values[i+1][j+1][COORD_TX], n->values[i+1][j+1][COORD_TY]);
         }
       fprintf(f,"na: Vector: (%f,%f,%f)\n",
	       n->normals[i][j][COORD_LX], n->normals[i][j][COORD_LY], n->normals[i][j][COORD_LZ]);
       fprintf(f,"nb: Vector: (%f,%f,%f)\n",
	       n->normals[i][j+1][COORD_LX], n->normals[i][j+1][COORD_LY], n->normals[i][j+1][COORD_LZ]);
       fprintf(f,"nc: Vector: (%f,%f,%f)\n",
	       n->normals[i+1][j+1][COORD_LX], n->normals[i+1][j+1][COORD_LY], n->normals[i+1][j+1][COORD_LZ]);
       fprintf(f,"transparencyA: RGB: (%f,%f,%f)\n", 
	       n->values[i][j][COORD_TR], n->values[i][j][COORD_TG], n->values[i][j][COORD_TB]);
       fprintf(f,"specularA: RGB: (%f,%f,%f)\n", 
	       n->values[i][j][COORD_SR], n->values[i][j][COORD_SG], n->values[i][j][COORD_SB]);
       fprintf(f,"diffuseA: RGB: (%f,%f,%f)\n", 
	       n->values[i][j][COORD_CR], n->values[i][j][COORD_CG], n->values[i][j][COORD_CB]);

       fprintf(f,"transparencyB: RGB: (%f,%f,%f)\n", 
	       n->values[i][j+1][COORD_TR], n->values[i][j+1][COORD_TG], n->values[i][j+1][COORD_TB]);
       fprintf(f,"specularB: RGB: (%f,%f,%f)\n", 
	       n->values[i][j+1][COORD_SR], n->values[i][j+1][COORD_SG], n->values[i][j+1][COORD_SB]);
       fprintf(f,"diffuseB: RGB: (%f,%f,%f)\n", 
	       n->values[i][j+1][COORD_CR], n->values[i][j+1][COORD_CG], n->values[i][j+1][COORD_CB]);
       
       fprintf(f,"transparencyC: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j+1][COORD_TR], n->values[i+1][j+1][COORD_TG], n->values[i+1][j+1][COORD_TB]);
       fprintf(f,"specularC: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j+1][COORD_SR], n->values[i+1][j+1][COORD_SG], n->values[i+1][j+1][COORD_SB]);
       fprintf(f,"diffuseC: RGB: (%f,%f,%f)\n", 
	       n->values[i+1][j+1][COORD_CR], n->values[i+1][j+1][COORD_CG], n->values[i+1][j+1][COORD_CB]);
       fprintf(f,"transparencyFactR: (1,%f)\n", n->values[i][j][COORD_FR]);
       fprintf(f,"transparencyFactG: (1,%f)\n", n->values[i][j][COORD_FG]);
       fprintf(f,"transparencyFactB: (1,%f)\n", n->values[i][j][COORD_FB]);
       fprintf(f,"normalDist: %f%%\n", n->values[i][j][COORD_ND]);
       fprintf(f,"specularFact: %f\n", n->values[i][j][COORD_SF]);
       fprintf(f,"faces: %d\n", n->faces);
       fprintf(f,"texture: %d\n", n->tid);
       fprintf(f,"}\n");
       k++;
    }
  }
}

void write_header(FILE* f, int num)
{
fprintf(f, "Screen: (800,600)\n");
fprintf(f, "Backup: 200\n");
fprintf(f, "Observer: Vertex: (0,0,-450)\n");
fprintf(f, "Light: Vertex: (50,50,250)\n");	
fprintf(f, "TexDirectory: textures\n");
fprintf(f, "NumTextures: 33\n");
fprintf(f, "nTriangles: %d\n",num);
}

void nurbs2dat(char* istr, char* ostr, int idx)
{ 
 FILE *in, *out;
 int n,nr,i;
 int ntris;
 NURBS* nurbs;
 in = fopen(istr, "r");
 if (!in) { printf("cannot open: '%s'\n", istr); exit(1); }
 out = fopen(ostr, "w");
 if (!out) { fclose(in); printf("cannot write to: '%s'\n", ostr); exit(1); }
 nr = fscanf(in, "nNURBS: %d\n", &n);
 if (nr != 1) panic("cannot read nNURBS");
 if (n <= 0) panic("bad nNURBS value");
 /* start */
 nurbs = (NURBS*)malloc(n*sizeof(NURBS));
 for (i=0;i<n;i++) read_nurbs(in, &nurbs[i]);
 ntris = 0;
 for (i=0;i<n;i++) 
   {
    nurbs[i].idx = ntris;
    ntris += nurbs[i].ntri;
/*    printf("NURBS: %d, nTris: %d, IDX: %d\n", i, nurbs[i].ntri, nurbs[i].idx);*/
   }
 if (lt)
   {
    write_header(out, idx+ntris);
    fprintf(out, "#all nurbs\n");
    fprintf(out, "ListTransform: [%d,%d]\n{\n}\n", idx, idx + ntris - 1);
    for (i=0;i<n;i++)
      {
       fprintf(out, "#nurbs %d\n", i);
       fprintf(out, "ListTransform: [%d,%d]\n{\n}\n", nurbs[i].idx+idx, nurbs[i].idx+idx+nurbs[i].ntri-1);
      }
   }
 for (i=0;i<n;i++) write_nurbs(out, &nurbs[i], idx);
 /* end */
 fclose(out);
 fclose(in);
 printf("Done.\n");
}


int main(int lb, char** par)
{
 if (lb != 8) { printf("%s infile.nurbs outfile.dat idx lt\n", par[0]); exit(1); }
 lt = atoi(par[4]);
 nurbs2dat(par[1], par[2], atoi(par[3]));
 return 0;
}

