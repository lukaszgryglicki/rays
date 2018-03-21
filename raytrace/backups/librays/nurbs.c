/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <GL/glut.h>
#include "libfparser.h"
#define PI  3.1415926F
#define PI2 6.2831853F
/*#define DEBUG*/
int wx,wy;
int cx,cy;
int n,p,m,s,si,sj;
int solid = 1;
#ifdef DEBUG
long long int bcnt;
#endif
GLfloat angX, angY, angZ;
GLfloat*** D;
GLfloat*** P;
GLfloat*** pts;
GLfloat**** tricubic = NULL;
GLfloat**** tri_mesh = NULL;
GLfloat*    tricubic_knot = NULL;
GLfloat*** Q;
GLfloat** N;
GLfloat** NI;
GLfloat* t;
GLfloat* knot;
GLfloat* ctlpts;
GLUnurbsObj *theNurb;
GLfloat scalef = 20.;
GLfloat precision;
GLfloat step = .1;
GLfloat grid = .5;
time_t t1,t2;
int fps;
int epsU = 4;
int epsV = 4;
int epsUprev = -1;
int epsVprev = -1;
int want_d = 0;
int want_p = 0;
int want_q = 0;
int auto_move = 1;
int want_nurb = 1;
int moving = 0;
int deform = 0;
int uvt = 0;
int moving_tri = 0;
int u_par = 0;
int v_par = 0;
int t_par = 0;
int meshD = 64;
int mapping = 0;
int want_trictl = 1;
enum {
    MENU_AUTO_MOVE,
    MENU_D,
    MENU_Q,
    MENU_P,
    MENU_SPLUS,
    MENU_SMINUS,
    MENU_UPLUS,
    MENU_UMINUS,
    MENU_VPLUS,
    MENU_VMINUS,
    MENU_NORMAL,
    MENU_WRITE_BSP,
    MENU_LOWER_DEGREE,
    MENU_HIGHER_DEGREE,
    MENU_QUIT
};
void regenerate_pts();
GLfloat*** matrix3(int siz1, int siz2, int siz3);
void free_matrix3(GLfloat*** m, int siz1, int siz2);
GLfloat**** matrix4(int siz1, int siz2, int siz3, int siz4);
void free_matrix4(GLfloat**** m, int siz1, int siz2, int siz3);
GLfloat** matrix(int siz);
void free_matrix(GLfloat** m, int siz);
GLfloat** invert_matrix(GLfloat** srcC, int dim);
void construct_Q();
void construct_P();
void compute_nu();
GLfloat tricube(GLfloat u, GLfloat v, GLfloat t, int xyz);
#ifdef DEBUG
void print_matrix3(GLfloat*** m, int siz1, int siz2, int siz3, char* name);
#endif

GLfloat tri_basis_bspline(GLfloat val, int dim, int off, int I, int J)
{
 register GLfloat up;
 register GLfloat down;
 GLfloat factor1;
 GLfloat factor2;
 register volatile int Ioff, Joff;
 factor1 = factor2 = 0.0;
 Ioff = I+off;
 Joff = J+off;
#ifdef DEBUG
 bcnt++;
#endif
 if (dim == 0)
   {
    if (val >= tricubic_knot[Ioff] && val <= tricubic_knot[Joff]) return 1.0;
    else return 0.0;
   }
 up = (val-tricubic_knot[Ioff]) * tri_basis_bspline(val, dim-1, off, I, J-1);
 down = tricubic_knot[Joff-1] - tricubic_knot[Ioff];
 if (down==0.0) factor1 = 0.0;
 else factor1 = up/down;
 up = (tricubic_knot[Joff]-val) * tri_basis_bspline(val, dim-1, off, I+1, J);
 down = tricubic_knot[Joff] - tricubic_knot[Ioff+1];
 if (down==0.0) factor2 += 0.0;
 else factor2 = up/down;
 return factor1+factor2;
}

GLfloat basis_bspline(GLfloat val, int dim, int off, int I, int J)
{
 register GLfloat up;
 register GLfloat down;
 GLfloat factor1;
 GLfloat factor2;
 register volatile int Ioff, Joff;
 factor1 = factor2 = 0.0;
 Ioff = I+off;
 Joff = J+off;
#ifdef DEBUG
 bcnt++;
#endif
 if (dim == 0)
   {
    if (val >= knot[Ioff] && val <= knot[Joff]) return 1.0;
    else return 0.0;
   }
 up = (val-knot[Ioff]) * basis_bspline(val, dim-1, off, I, J-1);
 down = knot[Joff-1] - knot[Ioff];
 if (down==0.0) factor1 = 0.0;
 else factor1 = up/down;
 up = (knot[Joff]-val) * basis_bspline(val, dim-1, off, I+1, J);
 down = knot[Joff] - knot[Ioff+1];
 if (down==0.0) factor2 += 0.0;
 else factor2 = up/down;
 return factor1+factor2;
}


GLfloat bsp(GLfloat t, int dim, int i)
{
 return basis_bspline(t, dim, dim, i-dim, i+1);
}

GLfloat tri_bsp(GLfloat t, int dim, int i)
{
 return tri_basis_bspline(t, dim, dim, i-dim, i+1);
}

void compute_trimesh()
{
 GLint i,j,k,l;
 GLfloat I,J,K;
 if (!deform) return;
 for (i=0;i<=meshD;i+=16)
   {
    I = (GLfloat)i / (GLfloat)meshD;
    for (j=0;j<=meshD;j+=16)
      {
       J = (GLfloat)j / (GLfloat)meshD;
       for (k=0;k<=meshD;k++)
         {
           K = (GLfloat)k / (GLfloat)meshD;
           for (l=0;l<s;l++)  tri_mesh[i][j][k][l] = tricube(I, J, K, l);
         }
      }
   }
 for (i=0;i<=meshD;i+=16)
   {
    I = (GLfloat)i / (GLfloat)meshD;
    for (j=0;j<=meshD;j++)
      {
       J = (GLfloat)j / (GLfloat)meshD;
       for (k=0;k<=meshD;k+=16)
         {
           K = (GLfloat)k / (GLfloat)meshD;
           for (l=0;l<s;l++)  tri_mesh[i][j][k][l] = tricube(I, J, K, l);
         }
      }
   }
 for (i=0;i<=meshD;i++)
   {
    I = (GLfloat)i / (GLfloat)meshD;
    for (j=0;j<=meshD;j+=16)
      {
       J = (GLfloat)j / (GLfloat)meshD;
       for (k=0;k<=meshD;k+=16)
         {
           K = (GLfloat)k / (GLfloat)meshD;
           for (l=0;l<s;l++)  tri_mesh[i][j][k][l] = tricube(I, J, K, l);
         }
      }
   }
}

void init_tricubic_engine()
{
 GLint i,j,k;
 tricubic_knot = (GLfloat*)malloc((4+3+1)*sizeof(GLfloat));
 for (i=0;i<4;i++) tricubic_knot[i] = 0.;
 for (i=4;i<8;i++) tricubic_knot[i] = 1.;
 tricubic = matrix4(4, 4, 4, s);
 tri_mesh = matrix4(meshD+1, meshD+1, meshD+1, s);
 for (i=0;i<4;i++)
 for (j=0;j<4;j++)
 for (k=0;k<4;k++)
   {
    tricubic[i][j][k][0] = (GLfloat)i / 3.;
    tricubic[i][j][k][1] = (GLfloat)j / 3.;
    tricubic[i][j][k][2] = (GLfloat)k / 3.;
   }
 compute_trimesh();
 printf("tricubic_init done.\n");
}


void help()
{
 printf("usage: {basename} d|n (l[w|o] file| f[w|o] dim npts def| r[w|o] dim npts | u[w|o] dim npts def1 def2 def3)\n");
 printf("d|n stands for use deform tricubic bezier box or not\n");
 printf("(l/r)button: menu\n");
 printf("1/2/3: toggle display of: D,Q,P\n");
 printf("4: auto_move toogle\n");
 printf("5/6/7/8: curve smoothness\n");
 printf("+/-: scaling\n");
 printf("p: select point\n");
 printf("z: tricubic noise by 3%%\n");
 printf("w/s/a/d/e/x: move selected point\n");
 printf("i/k/j/l/o/m: rotations\n");
 printf("\\: change mapping of tricube: u->v->t\n");
 printf("/:toggle display mode nurbs/tricube\n");
 printf(".:manually recompute trimesh\n");
 printf(",:toggle want_trictl\n");
 printf("y: start moving mapping tricubic (toggle)\n");
 printf("u: move to next tricubic point (according to \\)\n");
 printf("n/b: next/prev interpol point U.\n");
 printf("c/v: next/prev interpol point V.\n");
 printf("0/9: increase/decrease degree\n");
 printf("r: write BSpline surface\n");
 printf("`: write raylab.dat for RT\n");
 printf("only surface is deformed, control net remains unchanged.\n");
 printf("to see how tricubic bezier box deforms space E^3\n");
}


void free_points()
{
 if (pts)
  {
   if (epsU > 0 && epsV >0) free_matrix3(pts, epsU+1, epsV+1);
   pts = NULL;
  }
 if (ctlpts) free(ctlpts);
 ctlpts = NULL;
}


int modify(GLint* mod, GLint* modp, GLfloat c)
{
 if (solid)
   {
    if (c < 1) precision *= 1.5;
    else precision /= 1.5;
    gluNurbsProperty(theNurb, GLU_SAMPLING_TOLERANCE, precision);
   }
 if ((*mod) <= 4 && c < 1.) return 0;
 free_points();
 *modp = *mod;
 *mod  = (GLint)(*mod * c);
 return 1;
}


void next_point()
{
 sj++;
 if (sj == n+1) { si++; sj = 0; }
 if (si == n+1) si = 0;
}


void prev_point()
{
 sj--;
 if (sj == -1) { si--; sj = n; }
 if (si == -1) si = n;
}


void right_point()
{
 si++;
 if (si == n+1) { sj++; si = 0; }
 if (sj == n+1) sj = 0;
}


void left_point()
{
 si--;
 if (si == -1) { sj--; si = n; }
 if (sj == -1) sj = n;
}


void re_interpolate()
{
 free_matrix3(Q, n+1, n+1);
 construct_Q();
#ifdef DEBUG
 print_matrix3(Q, n+1, n+1, s, "Q");
#endif
 free_matrix3(P, n+1, n+1);
 construct_P();
#ifdef DEBUG
 print_matrix3(P, n+1, n+1, s, "P");
#endif
 pts = NULL;
 fflush(stdout);
 free_points();
 regenerate_pts();
}

void display_tricube()
{
 int i,j,k;
 if (!deform) return;
 for (i=0;i<4;i++)
 { 
  printf("\t***\n");
  for (j=0;j<4;j++)
    {
     for (k=0;k<4;k++) printf("(%1.3f,%1.3f,%1.3f)\t ", tricubic[i][j][k][0], tricubic[i][j][k][1], tricubic[i][j][k][2]);
     printf("\n");
    }
  printf("\n");
 }
}


void moving_point()
{
 if (moving)
   {
    re_interpolate();
    moving = 0;
   }
 else moving = 1;
}


void move_point(GLfloat x, GLfloat y, GLfloat z)
{
/* GLfloat gX, gY, gZ, tmp;*/
 if (moving)
   {
    D[si][sj][0] += x;
    D[si][sj][1] += y;
    D[si][sj][2] += z;
   }
 if (deform && moving_tri)
   {
    tricubic[u_par][v_par][t_par][0] += x/5.;
    tricubic[u_par][v_par][t_par][1] += y/5.;
    tricubic[u_par][v_par][t_par][2] += z/5.;
    printf("new tricubic for: (%d,%d,%d) is (%1.3f,%1.3f,%1.3f)\n", u_par, v_par, t_par,
		    tricubic[u_par][v_par][t_par][0],
		    tricubic[u_par][v_par][t_par][1],
		    tricubic[u_par][v_par][t_par][2]);
   }
 /*gX = (PI*angX)/180.;
 gY = (PI*angY)/180.;
 gZ = (PI*angZ)/180.;
 tmp  = y*cos(gX) - z*sin(gX);
 z    = y*sin(gX) + z*cos(gX);
 y    = tmp;
 tmp  = x*cos(gY) - z*sin(gY);
 z    = x*sin(gY) + z*cos(gY);
 x    = tmp;
 tmp  = x*cos(gZ) - y*sin(gZ);
 y    = x*sin(gZ) + y*cos(gZ);
 x    = tmp;*/
}


void write_bspline()
{
 int i,j;
 FILE* plik;
 char fname[32];
 i = 1;
 while (1)
 {
  sprintf(fname, "bspline_surface%d", i);
  plik = fopen(fname, "r");
  if (plik) { fclose(plik); i++; }
  else break;
 }
 plik = fopen(fname, "w");
 if (!plik) { printf("cannot write to: %s\n", fname); return; }
 printf("writing bspline surface definition to: %s\n", fname);
 fprintf(plik, "%d,%d,%d\n", n,p,s);
 for (i=0;i<=n;i++)
 for (j=0;j<=n;j++) fprintf(plik, "(%f,%f,%f)\n", D[i][j][0], D[i][j][1], D[i][j][2]);
 fclose(plik);
}


void full_re_interpolate()
{
 int i,j;
 if (knot) free(knot);
 knot = (GLfloat*)malloc((m+1)*sizeof(GLfloat));
 for (i=0;i<=p;i++) knot[i] = 0.;
 for (j=1;j<=n-p;j++)
    {
	knot[j+p] = 0.;
	for (i=j;i<=j+p-1;i++) knot[j+p] += t[i];
	knot[j+p] /= (GLfloat)p;
    }
 for (i=m-p;i<=m;i++) knot[i] = 1.;
 if (N) free_matrix(N, n+1);
 if (NI) free_matrix(NI, n+1);
 N = matrix(n+1);
 compute_nu();
 NI = invert_matrix(N, n+1);
 re_interpolate();
}


void lower_deg()
{
 p--;
 if (p < 1) { p = 1; return; }
 m = n+p+1;
 full_re_interpolate();
}


void higher_deg()
{
 p++;
 if (p > n) { p = n; return; }
 m = n+p+1;
 full_re_interpolate();
}

GLfloat nurbs(GLfloat u, GLfloat v, int xyz);

void calc_normal(double u, double v, double* x, double* y, double* z)
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
 
 epsu = 3e-4;
 epsv = 3e-4;

 if (1 - 2.*epsu <= u) epsu *= -1.;
 if (1 - 2.*epsv <= v) epsv *= -1.;
 
 x1 = nurbs(u, v, 0);
 x2 = nurbs(u+epsu, v, 0);
 x3 = nurbs(u, v+epsv, 0);
 y1 = nurbs(u, v, 1);
 y2 = nurbs(u+epsu, v, 1);
 y3 = nurbs(u, v+epsv, 1);
 z1 = nurbs(u, v, 2);
 z2 = nurbs(u+epsu, v, 2);
 z3 = nurbs(u, v+epsv, 2);
 
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

void write_raylab()
{
 FILE* plik;
 GLdouble x1,y1,z1;
 GLdouble x2,y2,z2;
 GLdouble x3,y3,z3;
 GLdouble U,V;
 GLint i,j,k;
 GLdouble u,v, s;
 GLdouble tr,tg,tb;
 GLdouble cr,cg,cb;
 GLdouble sr,sg,sb;
 GLdouble sf, tfr,tfg, tfb;
 GLdouble n1x,n1y,n1z;
 GLdouble n2x,n2y,n2z;
 GLdouble n3x,n3y,n3z;
 GLint faces;
 GLint tid;
 tr = tg = 0;
 tb = 1.;
 sr = sg = sb = .7;
 cr = 1;
 cg = .8;
 cb = 0.;
 sf = 255.;
 tfr = 1.4;
 tfg = 1.5;
 tfb = 1.6;
 faces = 1.;
 tid = 10;
 plik = fopen("raylab.dat", "w");
 if (!plik) return;
/* fprintf(plik, "%d\n", 2*epsU*epsV);*/
 k = 12;
 s = 220.;
  for (i=0;i<epsU;i++)
  {
   u = (GLfloat)i/(GLfloat)epsU;
   U = (GLfloat)(i+1)/(GLfloat)epsU;
   for (j=0;j<epsV;j++)
    {
     v = (GLfloat)j/(GLfloat)epsV;
     V = (GLfloat)(j+1)/(GLfloat)epsV;
     x1 = s*pts[i][j][0]; 
     y1 = s*pts[i][j][1]; 
     z1 = s*pts[i][j][2];
     x2 = s*pts[i+1][j+1][0];
     y2 = s*pts[i+1][j+1][1];
     z2 = s*pts[i+1][j+1][2];
     x3 = s*pts[i+1][j][0];
     y3 = s*pts[i+1][j][1];
     z3 = s*pts[i+1][j][2];

     calc_normal(u, v, &n1x, &n1y, &n1z);
     calc_normal(U, V, &n2x, &n2y, &n2z);
     calc_normal(U, v, &n3x, &n3y, &n3z);
     
     /*fprintf(plik, "%f %f %f\n", s*pts[i][j][0], s*pts[i][j][1], s*pts[i][j][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i+1][j+1][0], s*pts[i+1][j+1][1], s*pts[i+1][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i+1][j][0], s*pts[i+1][j][1], s*pts[i+1][j][2]);
 
     fprintf(plik, "%f %f %f\n", s*pts[i][j+1][0], s*pts[i][j+1][1], s*pts[i][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i+1][j+1][0], s*pts[i+1][j+1][1], s*pts[i+1][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i][j][0], s*pts[i][j][1], s*pts[i][j][2]);*/
       fprintf(plik,"Triangle: %d\n", k);
       fprintf(plik,"{\n");
       fprintf(plik,"a: Vertex: (%f,%f,%f)\n",x1,z1-350.,y1);
       fprintf(plik,"b: Vertex: (%f,%f,%f)\n",x2,z2-350.,y2);
       fprintf(plik,"c: Vertex: (%f,%f,%f)\n",x3,z3-350.,y3);
       fprintf(plik,"texA: TexCoord: (%f,%f)\n", u, v);
       fprintf(plik,"texB: TexCoord: (%f,%f)\n", U, V);
       fprintf(plik,"texC: TexCoord: (%f,%f)\n", U, v);
       fprintf(plik,"na: Vector: (%f,%f,%f)\n",n1x,n1z,n1y);
       fprintf(plik,"nb: Vector: (%f,%f,%f)\n",n2x,n2z,n2y);
       fprintf(plik,"nc: Vector: (%f,%f,%f)\n",n3x,n3z,n3y);
       fprintf(plik,"transparency: RGB: (%f,%f,%f)\n", tr, tg, tb);
       fprintf(plik,"specular: RGB: (%f,%f,%f)\n", sr, sg, sb);
       fprintf(plik,"diffuse: RGB: (%f,%f,%f)\n", cr, cg, cb);
       fprintf(plik,"surface: ABCD: (0,0,0,0)\n");
       fprintf(plik,"transparencyFactR: (1,%f)\n", tfr);
       fprintf(plik,"transparencyFactG: (1,%f)\n", tfg);
       fprintf(plik,"transparencyFactB: (1,%f)\n", tfb);
       fprintf(plik,"specularFact: %f\n", sf);
       fprintf(plik,"faces: %d\n",faces);
       fprintf(plik,"texture: %d\n", tid);
       fprintf(plik,"}\n");
     k++;
     x1 = s*pts[i][j][0];
     y1 = s*pts[i][j][1];
     z1 = s*pts[i][j][2];
     x2 = s*pts[i][j+1][0]; 
     y2 = s*pts[i][j+1][1]; 
     z2 = s*pts[i][j+1][2];
     x3 = s*pts[i+1][j+1][0];
     y3 = s*pts[i+1][j+1][1];
     z3 = s*pts[i+1][j+1][2];
     
     calc_normal(u, v, &n1x, &n1y, &n1z);
     calc_normal(u, V, &n2x, &n2y, &n2z);
     calc_normal(U, V, &n3x, &n3y, &n3z);
     
 
     /*fprintf(plik, "%f %f %f\n", s*pts[i][j+1][0], s*pts[i][j+1][1], s*pts[i][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i+1][j+1][0], s*pts[i+1][j+1][1], s*pts[i+1][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i][j][0], s*pts[i][j][1], s*pts[i][j][2]);*/
       fprintf(plik,"Triangle: %d\n", k);
       fprintf(plik,"{\n");
       fprintf(plik,"a: Vertex: (%f,%f,%f)\n",x1,z1-350.,y1);
       fprintf(plik,"b: Vertex: (%f,%f,%f)\n",x2,z2-350.,y2);
       fprintf(plik,"c: Vertex: (%f,%f,%f)\n",x3,z3-350.,y3);
       fprintf(plik,"texA: TexCoord: (%f,%f)\n", u, v);
       fprintf(plik,"texB: TexCoord: (%f,%f)\n", u, V);
       fprintf(plik,"texC: TexCoord: (%f,%f)\n", U, V);
       fprintf(plik,"na: Vector: (%f,%f,%f)\n",n1x,n1z,n1y);
       fprintf(plik,"nb: Vector: (%f,%f,%f)\n",n2x,n2z,n2y);
       fprintf(plik,"nc: Vector: (%f,%f,%f)\n",n3x,n3z,n3y);
       fprintf(plik,"transparency: RGB: (%f,%f,%f)\n", tr, tg, tb);
       fprintf(plik,"specular: RGB: (%f,%f,%f)\n", sr, sg, sb);
       fprintf(plik,"diffuse: RGB: (%f,%f,%f)\n", cr, cg, cb);
       fprintf(plik,"surface: ABCD: (0,0,0,0)\n");
       fprintf(plik,"transparencyFactR: (1,%f)\n", tfr);
       fprintf(plik,"transparencyFactG: (1,%f)\n", tfg);
       fprintf(plik,"transparencyFactB: (1,%f)\n", tfb);
       fprintf(plik,"specularFact: %f\n", sf);
       fprintf(plik,"faces: %d\n",faces);
       fprintf(plik,"texture: %d\n", tid);
       fprintf(plik,"}\n");
     k++;
    }
  }
 printf("writen %d triangles\n", k);
}

void write_triangulation()
{
 FILE* plik;
 GLint i,j,k;
 GLdouble u,v, s;
 plik = fopen("triangulation", "w");
 if (!plik) return;
 fprintf(plik, "%d\n", 2*epsU*epsV);
 k = 0;
 s = .1;
  for (i=0;i<epsU;i++)
  {
   u = (GLfloat)i/(GLfloat)epsU;
   for (j=0;j<epsV;j++)
    {
     v = (GLfloat)j/(GLfloat)epsV;
     
     fprintf(plik, "%f %f %f\n", s*pts[i][j][0], s*pts[i][j][1], s*pts[i][j][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i+1][j+1][0], s*pts[i+1][j+1][1], s*pts[i+1][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i+1][j][0], s*pts[i+1][j][1], s*pts[i+1][j][2]);
 
     fprintf(plik, "%f %f %f\n", s*pts[i][j+1][0], s*pts[i][j+1][1], s*pts[i][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i+1][j+1][0], s*pts[i+1][j+1][1], s*pts[i+1][j+1][2]);
     fprintf(plik, "%f %f %f\n", s*pts[i][j][0], s*pts[i][j][1], s*pts[i][j][2]);
     k++;
    }
  }
 printf("writen %d triangles\n", k);
}

void mapping_uvt()
{
 if (!deform) return;
 uvt++;
 if (uvt == 3) uvt = 0;
 printf("current_configuration of tricube: (%d,%d,%d) moving: %d, uvt=%d\n", u_par, v_par, t_par, moving_tri, uvt);
}

void moving_tricube()
{
 if (!deform) return;
 if (moving_tri)
   {
    re_interpolate();
    moving_tri = 0;
    display_tricube();
   }
 else moving_tri = 1;
 printf("current_configuration of tricube: (%d,%d,%d) moving: %d, uvt=%d\n", u_par, v_par, t_par, moving_tri, uvt);
}

void move_tripoint()
{
 if (!deform) return;
 if (uvt == 0) u_par ++;
 if (uvt == 1) v_par ++;
 if (uvt == 2) t_par ++;
 if (u_par == 4) u_par = 0;
 if (v_par == 4) v_par = 0;
 if (t_par == 4) t_par = 0;
 printf("current_configuration of tricube: (%d,%d,%d) moving: %d, uvt=%d\n", u_par, v_par, t_par, moving_tri, uvt);
}

void tricubic_noise()
{
 GLint i,j,k,l;
 if (!deform) return;
 for (i=0;i<4;i++)
 for (j=0;j<4;j++)
 for (k=0;k<4;k++)
 for (l=0;l<3;l++) tricubic[i][j][k][l] += (GLfloat)((rand() % 2001) - 1000) / 33333.;
 re_interpolate();
}

void change_mapping()
{
 mapping = !mapping;
}


void keyboard(unsigned char key, int x, int y)
{
 switch (key)
   {
        case 27: case 'q':  if (deform) { free_matrix4(tricubic, 4, 4, 4); free_matrix4(tri_mesh, meshD+1, meshD+1, meshD+1); } exit(0); break;
	case 'z': tricubic_noise(); break;
	case '`': write_raylab(); break;
	case 'h': help(); break;
        case 'r': write_bspline(); break;
	case 't': write_triangulation(); break;
	case 'd': move_point(step, 0., 0.); break;
	case 'a': move_point(-step, 0., 0.); break;
	case 's': move_point(0., -step, 0.); break;
	case 'w': move_point(0., step, 0.); break;
	case 'e': move_point(0., 0., step); break;
	case 'x': move_point(0., 0., -step); break;
        case 'i': angX -= 10.; break;
        case 'k': angX += 10.; break;
        case 'j': angY += 10.; break;
        case 'l': angY -= 10.; break;
        case 'o': angZ += 10.; break;
        case 'm': angZ -= 10.; break;
	case '\\': mapping_uvt(); break;
	case '/': change_mapping(); break;
	case '.': compute_trimesh(); break;
	case ',': want_trictl = ! want_trictl; break;		  
	case 'y': moving_tricube(); break;
	case 'u': move_tripoint(); break;
	case 'p': moving_point(); break;
	case 'n': next_point(); break;
	case 'b': prev_point(); break;
	case 'v': right_point(); break;
	case 'c': left_point(); break;
	case '1': want_d = ! want_d; break;
	case '2': want_q = ! want_q; break;
	case '3': want_p = ! want_p; break;
        case '4': auto_move = ! auto_move; break;
        case '-': scalef /= 1.1; break;
	case '=': case '+': scalef *= 1.1; break;
	case '9': lower_deg(); break;
	case '0': higher_deg(); break;
	case '5': if (modify(&epsU, &epsUprev, 2.))  regenerate_pts(); break;
	case '6': if (modify(&epsU, &epsUprev,.5))   regenerate_pts();  break;
	case '7': if (modify(&epsV, &epsVprev, 2.))  regenerate_pts();  break;
	case '8': if (modify(&epsV, &epsVprev, .5))  regenerate_pts();  break;
	case ' ': scalef = 20.; angY = angZ = 0.; angX = -180.; auto_move = 0; break;
   }
}


void menu(int sel)
{
 switch(sel)
 {
     case MENU_AUTO_MOVE: auto_move = ! auto_move; break;
     case MENU_D: want_d = ! want_d; break;
     case MENU_Q: want_q = ! want_q; break;
     case MENU_P: want_p = ! want_p; break;
     case MENU_SPLUS:  scalef *= 1.3; break;
     case MENU_SMINUS: scalef /= 1.3; break;
     case MENU_UPLUS:   if (modify(&epsU, &epsUprev, 2.)) regenerate_pts();  break;
     case MENU_UMINUS:  if (modify(&epsU, &epsUprev, .5)) regenerate_pts();  break;
     case MENU_VPLUS:   if (modify(&epsV, &epsVprev, 2.)) regenerate_pts();  break;
     case MENU_VMINUS:  if (modify(&epsV, &epsVprev, .5)) regenerate_pts();  break;
     case MENU_NORMAL: scalef = 1.; angX = angY = angZ = 0.; auto_move = 0; break;
     case MENU_WRITE_BSP: write_bspline(); break;
     case MENU_LOWER_DEGREE: lower_deg(); break;
     case MENU_HIGHER_DEGREE: higher_deg(); break;
     case MENU_QUIT: exit(0); break;
 }
}


void resize_scene(int w, int h)
{
 wx = w;
 wy = h;
 glViewport(0, 0, (GLsizei)wx, (GLsizei)wy);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 glFrustum(-1., 1., -1., 1., 1.5, 2000.);
 glMatrixMode(GL_MODELVIEW);
}


GLfloat nurbs(GLfloat u, GLfloat v, int xyz)
{
 register int i,j;
 GLfloat up,down;
 GLfloat factor;
 GLfloat b1,b2;
 down = 0.0;
 up = 0.0;
 for (i=0;i<=n;i++)
 {
 b1 = bsp(u, p, i);
 for (j=0;j<=n;j++)
   {
    b2 = bsp(v, p, j);
    factor = P[i][j][xyz]*b1*b2;
    up += factor;
    down += b1*b2;
   }
 }
 if (down==0.0) return 0.0;
 return up/down;
}

GLfloat tricube(GLfloat u, GLfloat v, GLfloat t, int xyz)
{
 register int i,j,k;
 GLfloat up,down;
 GLfloat factor;
 GLfloat b1,b2,b3;
 down = 0.0;
 up = 0.0;
 for (i=0;i<4;i++)
 {
  b1 = tri_bsp(u, 3, i);
  for (j=0;j<4;j++)
   {
    b2 = tri_bsp(v, 3, j);
    for (k=0;k<4;k++)
     {
      b3 = tri_bsp(t, 3, k);
      factor = tricubic[i][j][k][xyz]*b1*b2*b3;
      up += factor;
      down += b1*b2*b3;
     }
   }
 }
 if (down==0.0) return 0.0;
 return up/down;
}


void free_matrix3(GLfloat*** m, int siz1, int siz2)
{
 int i,j;
 if (siz1 <= 0 || siz2 <= 0) return;
 for (i=0;i<siz1;i++)
 for (j=0;j<siz2;j++) free(m[i][j]);
 for (i=0;i<siz1;i++) free(m[i]);
 free(m);
}


void free_matrix4(GLfloat**** m, int siz1, int siz2, int siz3)
{
 int i,j,k;
 if (siz1 <= 0 || siz2 <= 0 || siz3 <= 0) return;
 for (i=0;i<siz1;i++)
 for (j=0;j<siz2;j++)
 for (k=0;k<siz3;k++) free(m[i][j][k]);
 for (i=0;i<siz1;i++)
 for (j=0;j<siz2;j++) free(m[i][j]);
 for (i=0;i<siz1;i++) free(m[i]);
 free(m);
}


void init_ctlpts()
{
 int i,j;
  ctlpts = (GLfloat*)malloc(((n+1)*(n+1)*3+1)*sizeof(GLfloat));
  for (i=0;i<=n;i++)
  for (j=0;j<=n;j++)
    {
     ctlpts[3*i*(n+1)+3*j] = P[i][j][0];
     ctlpts[3*i*(n+1)+3*j+1] = P[i][j][1];
     ctlpts[3*i*(n+1)+3*j+2] = P[i][j][2];
    }
}

GLfloat deform_space(GLfloat a, GLfloat b, GLfloat c, GLint k)
{
 GLfloat u;
 GLfloat v;
 GLfloat t;
 GLfloat Nv;
 u = a - (GLfloat)((int)a);
 v = b - (GLfloat)((int)b);
 t = c - (GLfloat)((int)c);
 if (u < 0.) u += 1.;
 if (v < 0.) v += 1.;
 if (t < 0.) t += 1.;
 Nv = tricube(u, v, t, k);
/* printf("%f %f %f tricube(%f,%f,%f,%d) is %f!\n",a,b,c,u,v,t,k, Nv);*/
 switch (k)
   {
    case 0: Nv += a - u; break;
    case 1: Nv += b - v; break;
    case 2: Nv += c - t; break;
    default: return 0;
   }
 return Nv;
}


void regenerate_pts()
{
 GLint i,j,k;
 GLfloat u,v;
 GLfloat a,b,c;
 GLint maxE;
 maxE = epsU>epsV?epsU:epsV;
 pts = NULL;
 pts = matrix3(maxE+1, maxE+1, s);
 if (!solid)
 {
 for (i=0;i<=epsU;i++)
    {
     u = (GLfloat)i/(GLfloat)epsU;
     for (j=0;j<=epsV;j++)
       {
	v = (GLfloat)j/(GLfloat)epsV;
        for (k=0;k<s;k++)
	{
	  pts[i][j][k] = nurbs(u,v,k);
	}
	if (deform) 
	  {
	   a = pts[i][j][0];
	   b = pts[i][j][1];
	   c = pts[i][j][2];
	   pts[i][j][0] = deform_space(a, b, c, 0);
	   pts[i][j][1] = deform_space(a, b, c, 1);
	   pts[i][j][2] = deform_space(a, b, c, 2);
	  }
       }
      }
#ifdef DEBUG
  printf("bsp#: %lld\n", bcnt);
#endif
 }
 else init_ctlpts();
}


void time_counter()
{
 char tstr[32];
 if (t1 == (time_t)0)
   {
    time(&t1);
    time(&t2);
    return;
   }
 fps++;
 time(&t2);
 if (t2 > t1)
   {
    sprintf(tstr, "N=(%dx%d), Dim=%d, GLScene FPS: %d", n+1, n+1, p, fps/(int)(t2-t1));
    t1 = t2;
    glutSetWindowTitle(tstr);
    fps = 0;
   }
}


void render_scene_solid(void)
{
 int i,j;
 GLfloat v;
 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 /*glPolygonMode(GL_FRONT, GL_LINE);
 glPolygonMode(GL_BACK, GL_LINE);*/
 glLoadIdentity();
 gluLookAt(0., 0., 100., 0. , 0., 0., 0., 100., 0.);
 glScaled(scalef, scalef, scalef);
 glRotatef(angX, 1., 0. , 0.);
 glRotatef(angY, 0., 1. , 0.);
 glRotatef(angZ, 0., 0. , 1.);
 for (i=0;i<=n;i++)
 for (j=0;j<=n;j++)
 {
 if (i == si && j == sj) { glColor3d(1., 1., 0.); v = grid/5.; }
 else { glColor3d(0., .4, .9); v = grid/15.; }
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]+v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]-v);
 glEnd();
 }
 if (want_d)
 {
 glColor3d(1., 0., 0.);
 glBegin(GL_QUADS);
 for (i=0;i<n;i++)
 for (j=0;j<n;j++)
   {
    glVertex3d(D[i][j][0], D[i][j][1], D[i][j][2]);
    glVertex3d(D[i+1][j][0], D[i+1][j][1], D[i+1][j][2]);
    glVertex3d(D[i+1][j+1][0], D[i+1][j+1][1], D[i+1][j+1][2]);
    glVertex3d(D[i][j+1][0], D[i][j+1][1], D[i][j+1][2]);
   }
 glEnd();
 }
 if (want_q)
 {
 glColor3d(0., 1., 0.);
 glBegin(GL_QUADS);
 for (i=0;i<n;i++)
 for (j=0;j<n;j++)
   {
    glVertex3d(Q[i][j][0], Q[i][j][1], Q[i][j][2]);
    glVertex3d(Q[i+1][j][0], Q[i+1][j][1], Q[i+1][j][2]);
    glVertex3d(Q[i+1][j+1][0], Q[i+1][j+1][1], Q[i+1][j+1][2]);
    glVertex3d(Q[i][j+1][0], Q[i][j+1][1], Q[i][j+1][2]);
   }
 glEnd();
 }
 if (want_p)
 {
 glColor3d(0., 0., 1.);
 glBegin(GL_QUADS);
 for (i=0;i<n;i++)
 for (j=0;j<n;j++)
   {
    glVertex3d(P[i][j][0], P[i][j][1], P[i][j][2]);
    glVertex3d(P[i+1][j][0], P[i+1][j][1], P[i+1][j][2]);
    glVertex3d(P[i+1][j+1][0], P[i+1][j+1][1], P[i+1][j+1][2]);
    glVertex3d(P[i][j+1][0], P[i][j+1][1], P[i][j+1][2]);
   }
 glEnd();
 }
 if (want_nurb)
 {
   gluBeginSurface(theNurb);
   gluNurbsSurface(theNurb, m+1, knot, m+1, knot, (n+1)*3, 3, ctlpts, p+1, p+1, GL_MAP2_VERTEX_3);
   gluEndSurface(theNurb);
 }
 if (auto_move)
  {
   angX += 2.3;
   angY += 2.5;
   angZ += 2.7;
  }
 time_counter();
 glFlush();
 glutSwapBuffers();
}

void glCube(GLfloat x, GLfloat y, GLfloat z, GLfloat v)
{
 glBegin(GL_LINE_STRIP);
   glVertex3d(x-v, y-v, z-v);
   glVertex3d(x-v, y+v, z-v);
   glVertex3d(x+v, y+v, z-v);
   glVertex3d(x+v, y-v, z-v);
   glVertex3d(x-v, y-v, z-v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(x-v, y-v, z+v);
   glVertex3d(x-v, y+v, z+v);
   glVertex3d(x+v, y+v, z+v);
   glVertex3d(x+v, y-v, z+v);
   glVertex3d(x-v, y-v, z+v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(x-v, y-v, z-v);
   glVertex3d(x-v, y+v, z-v);
   glVertex3d(x-v, y+v, z+v);
   glVertex3d(x-v, y-v, z+v);
   glVertex3d(x-v, y-v, z-v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(x+v, y-v, z-v);
   glVertex3d(x+v, y+v, z-v);
   glVertex3d(x+v, y+v, z+v);
   glVertex3d(x+v, y-v, z+v);
   glVertex3d(x+v, y-v, z-v);
 glEnd();
}

void render_scene_tricube(void)
{
 int i,j,k;
 GLfloat v;
 glClear(GL_COLOR_BUFFER_BIT);
 glPolygonMode(GL_FRONT, GL_LINE);
 glPolygonMode(GL_BACK, GL_LINE);
 glLoadIdentity();
 gluLookAt(-0., 0., 5., 0. , 0., 0., 0., 5., 0.);
 glScaled(scalef/6., scalef/6., scalef/6.);
 glRotatef(angX, 1., 0. , 0.);
 glRotatef(angY, 0., 1. , 0.);
 glRotatef(angZ, 0., 0. , 1.);
 glTranslated(-.5, -.5, -.5);
 if (want_trictl)
 {
  glColor3d(1., 0., 0.);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  for (k=0;k<4;k++)
   {
    if (i == u_par && j == v_par && k == t_par) v = .03; else v = .01; 
    glCube(tricubic[i][j][k][0], tricubic[i][j][k][1], tricubic[i][j][k][2], v);
   }
 }
 glColor3d(0., 0., 1.);
 glBegin(GL_LINES);
 for (i=0;i<=meshD;i+=16)
 for (j=0;j<=meshD;j+=16)
 for (k=0;k<meshD;k++)
   {
    glVertex3d(tri_mesh[i][j][k][0], tri_mesh[i][j][k][1], tri_mesh[i][j][k][2]);
    glVertex3d(tri_mesh[i][j][k+1][0], tri_mesh[i][j][k+1][1], tri_mesh[i][j][k+1][2]);
   }
 for (i=0;i<=meshD;i+=16)
 for (j=0;j<meshD;j++)
 for (k=0;k<=meshD;k+=16)
   {
    glVertex3d(tri_mesh[i][j][k][0], tri_mesh[i][j][k][1], tri_mesh[i][j][k][2]);
    glVertex3d(tri_mesh[i][j+1][k][0], tri_mesh[i][j+1][k][1], tri_mesh[i][j+1][k][2]);
   }
 for (i=0;i<meshD;i++)
 for (j=0;j<=meshD;j+=16)
 for (k=0;k<=meshD;k+=16)
   {
    glVertex3d(tri_mesh[i][j][k][0], tri_mesh[i][j][k][1], tri_mesh[i][j][k][2]);
    glVertex3d(tri_mesh[i+1][j][k][0], tri_mesh[i+1][j][k][1], tri_mesh[i+1][j][k][2]);
   }
 glEnd();
 if (auto_move)
  {
   angX += 2.3;
   angY += 2.5;
   angZ += 2.7;
  }
 time_counter();
 glFlush();
 glutSwapBuffers();

}


void render_scene_nurbs(void)
{
 int i,j;
 GLfloat u,v;
 glClear(GL_COLOR_BUFFER_BIT);
 /*glEnable(GL_BLEND);
 glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_COLOR);*/
 glPolygonMode(GL_FRONT, GL_LINE);
 glPolygonMode(GL_BACK, GL_LINE);
 glLoadIdentity();
 gluLookAt(0., 0., 100., 0. , 0., 0., 0., 100., 0.);
 glScaled(scalef, scalef, scalef);
 glRotatef(angX, 1., 0. , 0.);
 glRotatef(angY, 0., 1. , 0.);
 glRotatef(angZ, 0., 0. , 1.);
 for (i=0;i<=n;i++)
 for (j=0;j<=n;j++)
 {
 if (i == si && j == sj) { glColor3d(1., 1., 0.); v = grid/5.; }
 else { glColor3d(0., .4, .9); v = grid/15.; }
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]+v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]-v, D[i][j][1]-v, D[i][j][2]-v);
 glEnd();
 glBegin(GL_LINE_STRIP);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]-v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]+v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]+v);
   glVertex3d(D[i][j][0]+v, D[i][j][1]-v, D[i][j][2]-v);
 glEnd();
 }
 if (want_d)
 {
 glColor3d(1., 0., 0.);
 glBegin(GL_QUADS);
 for (i=0;i<n;i++)
 for (j=0;j<n;j++)
   {
    glVertex3d(D[i][j][0], D[i][j][1], D[i][j][2]);
    glVertex3d(D[i+1][j][0], D[i+1][j][1], D[i+1][j][2]);
    glVertex3d(D[i+1][j+1][0], D[i+1][j+1][1], D[i+1][j+1][2]);
    glVertex3d(D[i][j+1][0], D[i][j+1][1], D[i][j+1][2]);
   }
 glEnd();
 }
 if (want_q)
 {
 glColor3d(0., 1., 0.);
 glBegin(GL_QUADS);
 for (i=0;i<n;i++)
 for (j=0;j<n;j++)
   {
    glVertex3d(Q[i][j][0], Q[i][j][1], Q[i][j][2]);
    glVertex3d(Q[i+1][j][0], Q[i+1][j][1], Q[i+1][j][2]);
    glVertex3d(Q[i+1][j+1][0], Q[i+1][j+1][1], Q[i+1][j+1][2]);
    glVertex3d(Q[i][j+1][0], Q[i][j+1][1], Q[i][j+1][2]);
   }
 glEnd();
 }
 if (want_p)
 {
 glColor3d(0., 0., 1.);
 glBegin(GL_QUADS);
 for (i=0;i<n;i++)
 for (j=0;j<n;j++)
   {
    glVertex3d(P[i][j][0], P[i][j][1], P[i][j][2]);
    glVertex3d(P[i+1][j][0], P[i+1][j][1], P[i+1][j][2]);
    glVertex3d(P[i+1][j+1][0], P[i+1][j+1][1], P[i+1][j+1][2]);
    glVertex3d(P[i][j+1][0], P[i][j+1][1], P[i][j+1][2]);
   }
 glEnd();
 }
 if (want_nurb)
 {
  glColor3d(.9, .6, 1);
  for (i=0;i<epsU;i++)
  {
   u = (GLfloat)i/(GLfloat)epsU;
   for (j=0;j<epsV;j++)
    {
     glBegin(GL_LINE_STRIP);
     v = (GLfloat)j/(GLfloat)epsV;
     glVertex3d(pts[i][j][0], pts[i][j][1], pts[i][j][2]);
     glVertex3d(pts[i+1][j][0], pts[i+1][j][1], pts[i+1][j][2]);
     glVertex3d(pts[i+1][j+1][0], pts[i+1][j+1][1], pts[i+1][j+1][2]);
     glVertex3d(pts[i][j+1][0], pts[i][j+1][1], pts[i][j+1][2]);
     glVertex3d(pts[i][j][0], pts[i][j][1], pts[i][j][2]);
     glVertex3d(pts[i+1][j+1][0], pts[i+1][j+1][1], pts[i+1][j+1][2]);
     glEnd();
    }
  }
 }
 if (auto_move)
  {
   angX += 2.3;
   angY += 2.5;
   angZ += 2.7;
  }
 time_counter();
 glFlush();
 glutSwapBuffers();
}

void render_scene(void)
{
 if (!mapping || !deform) render_scene_nurbs();
 else render_scene_tricube();
}


void anim(void)
{
    glutPostRedisplay();
}


void visible(int vis)
{
    if (vis == GLUT_VISIBLE) glutIdleFunc(anim);
    else                     glutIdleFunc(NULL);
}


GLfloat* vector(int siz)
{
 if (siz <= 0) return NULL;
 return (GLfloat*)malloc(siz*sizeof(GLfloat));
}


GLfloat** matrix(int siz)
{
 GLfloat** mem;
 int i;
 if (siz <= 0) return NULL;
 mem = (GLfloat**)malloc(siz*sizeof(GLfloat*));
 for (i=0;i<siz;i++) mem[i] = vector(siz);
 return mem;
}


GLfloat** matrix2(int siz1, int siz2)
{
 GLfloat** mem;
 int i;
 if (siz1 <= 0 || siz2 <= 0) return NULL;
 mem = (GLfloat**)malloc(siz1*sizeof(GLfloat*));
 for (i=0;i<siz1;i++) mem[i] = vector(siz2);
 return mem;
}


GLfloat*** matrix3(int siz1, int siz2, int siz3)
{
 GLfloat*** mem;
 int i;
 if (siz1 <= 0 || siz2 <= 0 || siz3 <= 0) return NULL;
 mem = (GLfloat***)malloc(siz1*sizeof(GLfloat**));
 for (i=0;i<siz1;i++) mem[i] = matrix2(siz2, siz3);
 return mem;
}

GLfloat**** matrix4(int siz1, int siz2, int siz3, int siz4)
{
 GLfloat**** mem;
 int i;
 if (siz1 <= 0 || siz2 <= 0 || siz3 <= 0 || siz4 <= 0) return NULL;
 mem = (GLfloat****)malloc(siz1*sizeof(GLfloat***));
 for (i=0;i<siz1;i++) mem[i] = matrix3(siz2, siz3, siz4);
 return mem;
}

#ifdef DEBUG

void print_vector(GLfloat* v, int siz, char* name)
{
 int j;
 if (name) printf("%s(%p) {", name, (void*)v);
 else printf("NONAME(%p) {\n", (void*)v);
 for (j=0;j<siz;j++) printf("%f\t", v[j]);
 printf("}\n");
}


void print_matrix(GLfloat** m, int siz, char* name)
{
 int i,j;
 if (name) printf(" ==== %s(%p) ====\n", name, (void*)m);
 else printf(" ==== NONAME(%p) ====\n", (void*)m);
 for (i=0;i<siz;i++)
   {
    for (j=0;j<siz;j++) printf("%f\t", m[i][j]);
    printf("\n");
   }
}


void print_matrix2(GLfloat** m, int siz1, int siz2, char* name)
{
 int i,j;
 if (name) printf(" ==== %s(%p) ====\n", name, (void*)m);
 else printf(" ==== NONAME(%p) ====\n", (void*)m);
 for (i=0;i<siz1;i++)
   {
    for (j=0;j<siz2;j++) printf("%f\t", m[i][j]);
    printf("\n");
   }
}


void print_matrix3(GLfloat*** m, int siz1, int siz2, int siz3, char* name)
{
 int i,j,k;
 if (name) printf(" ==== %s(%p) ====\n", name, (void*)m);
 else printf(" ==== NONAME(%p) ====\n", (void*)m);
 for (i=0;i<siz1;i++)
   {
    for (j=0;j<siz2;j++)
      {
       printf("(");
       for (k=0;k<siz3-1;k++) printf("%2.2f,", m[i][j][k]);
       printf("%2.2f)\t", m[i][j][siz3-1]);
      }
    printf("\n");
   }
}

#endif

void random_curve_data(int dim, int npts)
{
    int i,j;
    float dx, dy, sx, sy;
    time_t tm;
    int K;
    if (dim < 1) exit(26);
    time(&tm);
    srand((int)tm);
    if (npts < dim) { printf("Not enough points for dimnesion.\n"); exit(39); }
    sx = -((float)npts)*(grid/2.);
    sy = -((float)npts)*(grid/2.);
    dx = grid;
    dy = grid;
    s = 3;
    p = dim;
    n = npts;
    epsU = n;
    epsV = n;
    D = matrix3(n+1, n+1, s);
    t = (GLfloat*)malloc((n+1)*sizeof(GLfloat));
    K = grid * 300.;
    for (i=0;i<=n;i++)
    {
	for (j=0;j<=n;j++)
	{
	 D[i][j][0] = sx + dx*i;
	 D[i][j][1] = sy + dy*j;
	 D[i][j][2] = (((float)((rand()%(2*K+1))-K))*sqrt((float)p))/800.;
	}
    }
    for (i=0;i<=n;i++) t[i] = (GLfloat)i/(GLfloat)n;
    m = n+p+1;
    knot = (GLfloat*)malloc((m+1)*sizeof(GLfloat));
    for (i=0;i<=p;i++) knot[i] = 0.;
    for (j=1;j<=n-p;j++)
    {
	knot[j+p] = 0.;
	for (i=j;i<=j+p-1;i++) knot[j+p] += t[i];
	knot[j+p] /= (GLfloat)p;
    }
    for (i=m-p;i<=m;i++) knot[i] = 1.;
}


void function3_curve_data(char* fdef, char* fdef2, char* fdef3, int dim, int npts)
{
    int i,j;
    float dx, dy, sx, sy;
    if (dim < 1 || !fdef) exit(22);
    fpar_function(fdef);
    if (!fpar2_ok()) { printf("Parser error!\n"); exit(23); }
    fpar_free();
    fpar_function(fdef2);
    if (!fpar2_ok()) { printf("Parser error!\n"); exit(23); }
    fpar_free();
    fpar_function(fdef3);
    if (!fpar2_ok()) { printf("Parser error!\n"); exit(23); }
    fpar_free();
    fpar_function(fdef);
    if (!fpar2_ok()) { printf("Parser error!\n"); exit(23); }
    if (npts < dim) { printf("Not enough points for dimnesion.\n"); exit(39); }
    sx = -((float)npts)*(grid/2.);
    sy = -((float)npts)*(grid/2.);
    dx = grid;
    dy = grid;
    s = 3;
    p = dim;
    n = npts;
    epsU = n;
    epsV = n;
    D = matrix3(n+1, n+1, s);
    t = (GLfloat*)malloc((n+1)*sizeof(GLfloat));
    for (i=0;i<=n;i++) for (j=0;j<=n;j++) D[i][j][0] = 2.*grid*fpar2_f((1./(2.*grid))*(sx+dx*i), (1./(2.*grid))*(sy+dy*j));
    fpar_free();
    fpar_function(fdef2);
    if (!fpar2_ok()) { printf("Parser error!\n"); exit(23); }
    for (i=0;i<=n;i++) for (j=0;j<=n;j++) D[i][j][1] = 2.*grid*fpar2_f((1./(2.*grid))*(sx+dx*i), (1./(2.*grid))*(sy+dy*j));
    fpar_free();
    fpar_function(fdef3);
    if (!fpar2_ok()) { printf("Parser error!\n"); exit(23); }
    for (i=0;i<=n;i++) for (j=0;j<=n;j++) D[i][j][2] = 2.*grid*fpar2_f((1./(2.*grid))*(sx+dx*i), (1./(2.*grid))*(sy+dy*j));
    for (i=0;i<=n;i++) t[i] = (GLfloat)i/(GLfloat)n;
    fpar_free();
    m = n+p+1;
    knot = (GLfloat*)malloc((m+1)*sizeof(GLfloat));
    for (i=0;i<=p;i++) knot[i] = 0.;
    for (j=1;j<=n-p;j++)
    {
	knot[j+p] = 0.;
	for (i=j;i<=j+p-1;i++) knot[j+p] += t[i];
	knot[j+p] /= (GLfloat)p;
    }
    for (i=m-p;i<=m;i++) knot[i] = 1.;
}


void function_curve_data(char* fdef, int dim, int npts)
{
    int i,j;
    float dx, dy, sx, sy;
    if (dim < 1 || !fdef) exit(22);
    fpar_function(fdef);
    if (!fpar2_ok()) { printf("Parser error!\n"); exit(23); }
    if (npts < dim) { printf("Not enough points for dimnesion.\n"); exit(39); }
    sx = -((float)npts)*(grid/2.);
    sy = -((float)npts)*(grid/2.);
    dx = grid;
    dy = grid;
    s = 3;
    p = dim;
    n = npts;
    epsU = n;
    epsV = n;
    D = matrix3(n+1, n+1, s);
    t = (GLfloat*)malloc((n+1)*sizeof(GLfloat));
    for (i=0;i<=n;i++)
    {
	for (j=0;j<=n;j++)
	{
	 D[i][j][0] = sx + dx*i;
	 D[i][j][1] = sy + dy*j;
	 D[i][j][2] = 2.*grid*fpar2_f((1./(2.*grid))*(sx+dx*i), (1./(2.*grid))*(sy+dy*j));
	}
    }
    fpar_free();
    for (i=0;i<=n;i++) t[i] = (GLfloat)i/(GLfloat)n;
    m = n+p+1;
    knot = (GLfloat*)malloc((m+1)*sizeof(GLfloat));
    for (i=0;i<=p;i++) knot[i] = 0.;
    for (j=1;j<=n-p;j++)
    {
	knot[j+p] = 0.;
	for (i=j;i<=j+p-1;i++) knot[j+p] += t[i];
	knot[j+p] /= (GLfloat)p;
    }
    for (i=m-p;i<=m;i++) knot[i] = 1.;
}


void read_curve_data(char* fn)
{
    int i,j,k;
    FILE* plik;
    plik = fopen(fn, "r");
    if (!plik) { printf("cannot read from: %s\n", fn); exit(9); }
    fscanf(plik,"%d,%d,%d\n", &n, &p, &s);
    if (n < 2 || n > 128) exit(2);
    if (p < 1 || p > 12)  exit(3);
    if ( p < n && (p % 2)) p--;
    epsU = n;
    epsV = n;
    D = matrix3(n+1, n+1, s);
    t = (GLfloat*)malloc((n+1)*sizeof(GLfloat));
    for (i=0;i<=n;i++)
    {
	for (j=0;j<=n;j++)
	{
	fscanf(plik,"(");
	for (k=0;k<s-1;k++) fscanf(plik,"%f,", &D[i][j][k]);
	fscanf(plik,"%f)\n", &D[i][j][s-1]);
	}
    }
    for (i=0;i<=n;i++) t[i] = (GLfloat)i/(GLfloat)n;
    m = n+p+1;
    knot = (GLfloat*)malloc((m+1)*sizeof(GLfloat));
    for (i=0;i<=p;i++) knot[i] = 0.;
    for (j=1;j<=n-p;j++)
    {
	knot[j+p] = 0.;
	for (i=j;i<=j+p-1;i++) knot[j+p] += t[i];
	knot[j+p] /= (GLfloat)p;
    }
    for (i=m-p;i<=m;i++) knot[i] = 1.;
}


void copy_matrix(GLfloat** dst, GLfloat** src, int siz)
{
 int i,j;
 if (!src || !dst || siz<= 0) return;
 for (i=0;i<siz;i++) for (j=0;j<siz;j++) dst[i][j] = src[i][j];
}


void I_matrix(GLfloat** dst, int siz)
{
 int i,j;
 if (!dst || siz<= 0) return;
 for (i=0;i<siz;i++) for (j=0;j<siz;j++)
   {
    if (i == j) dst[i][j] = 1.;
    else dst[i][j] = 0.;
   }
}


int try_swap(GLfloat** m, int idx, int dim)
{
 register int x;
 for (x=idx;x<dim;x++) if (m[idx][x]) return x;
 return -1;
}


GLfloat** invert_matrix(GLfloat** srcC, int dim)
{
 GLfloat** src, **dst;
 GLfloat div, pom;
 register int x,k;
 int i,swit;
 GLfloat* vectmp;
 src = matrix(dim);
 dst = matrix(dim);
 vectmp = vector(dim);
 copy_matrix(src, srcC, dim);
 I_matrix(dst, dim);
 for (i=0;i<dim;i++)
   {
    div = src[i][i];
    if (div == 0.0)
      {
       swit = try_swap(src, i, dim);
       if (swit < 0)
         {
	   printf("UNINVERTABLE_MATRIX\n");
#ifdef DEBUG
	   print_matrix(srcC, dim, "src");
#endif
	   exit(1);
	 }
       for (x=0;x<dim;x++) vectmp[x]    = src[x][i];
       for (x=0;x<dim;x++) src[x][i]    = src[x][swit];
       for (x=0;x<dim;x++) src[x][swit] = vectmp[x];
       for (x=0;x<dim;x++) vectmp[x]    = dst[x][i];
       for (x=0;x<dim;x++) dst[x][i]    = dst[x][swit];
       for (x=0;x<dim;x++) dst[x][swit] = vectmp[x];
       div = src[i][i];
      }
    for (x=0;x<dim;x++)
      {
       src[x][i] /= div;
       dst[x][i] /= div;
      }
    for (k=0;k<dim;k++)
      {
       pom = src[i][k];
       if (k-i)
         {
          for (x=0;x<dim;x++) src[x][k] -= pom* src[x][i];
          for (x=0;x<dim;x++) dst[x][k] -= pom* dst[x][i];
         }
      }
   }
 return dst;
}


GLfloat** matrix_mul(GLfloat** m, GLfloat** n, int ma, int mb, int na, int nb)
{
 GLfloat** dst;
 register int k,j,i;
 if (ma <= 0 || mb  <= 0 || na <= 0 || nb <=0 || mb != na) return NULL;
 if (!m || !n) return NULL;
 dst = (GLfloat**)malloc(ma*sizeof(GLfloat*));
 if (!dst) return NULL;
 for (i=0;i<ma;i++) dst[i] = (GLfloat*)malloc(nb*sizeof(GLfloat));
 for (i=0;i<ma;i++)
 for (j=0;j<nb;j++)
    {
     dst[i][j] = 0.0;
     for (k=0;k<mb;k++) dst[i][j] += m[i][k] * n[k][j];
    }
 return dst;
}


void free_matrix(GLfloat** m, int siz)
{
 int i;
 if (siz <= 0) return;
 for (i=0;i<siz;i++) free(m[i]);
 free(m);
}


void compute_nu()
{
 int i,j;
 GLfloat sum;
 for (i=0;i<=n;i++)
 for (j=0;j<=n;j++) N[i][j] = bsp(t[i], p, j);
 for (i=0;i<=n;i++)
  {
   sum = 0.;
   for (j=0;j<=n;j++) sum += N[i][j];
   for (j=0;j<=n;j++) N[i][j] /= sum;
  }
}


void compute_nv()
{
 int i,j;
 GLfloat sum;
 free_matrix(N, n+1);
 N = matrix(n+1);
 for (i=0;i<=n;i++) for (j=0;j<=n;j++) N[i][j] = bsp(t[i], p, j);
 for (i=0;i<=n;i++)
  {
   sum = 0.;
   for (j=0;j<=n;j++) sum += N[i][j];
   for (j=0;j<=n;j++) N[i][j] /= sum;
  }
}


void construct_Q()
{
 int i,j,k;
 GLfloat*** tmp;
 tmp = (GLfloat***)malloc((n+1)*sizeof(GLfloat**));
 Q = matrix3(n+1, n+1, s);
 for (i=0;i<=n;i++)
   {
    tmp[i] = matrix_mul(NI, D[i], n+1, n+1, n+1, s);
#ifdef DEBUG
    print_matrix2(tmp[i], n+1, s, "Qi");
#endif
   }
 for (i=0;i<=n;i++)
 for (j=0;j<=n;j++)
 for (k=0;k<s;k++) Q[j][i][k] = tmp[i][j][k];
 free_matrix3(tmp, n+1, n+1);
}


void construct_P()
{
 int i;
 P =  (GLfloat***)malloc((n+1)*sizeof(GLfloat**));
 for (i=0;i<=n;i++)
   {
    P[i] = matrix_mul(NI, Q[i], n+1, n+1, n+1, s);
#ifdef DEBUG
    print_matrix2(P[i], n+1, s, "Pi");
#endif
   }
}


void Init(char* fn, char* fdef, char* fdef2, char* fdef3, int dim, int npts)
{
 if (fn) read_curve_data(fn);
 else if (fdef && !fdef2 && !fdef3) function_curve_data(fdef, dim, npts);
 else if (fdef && fdef2 && fdef3)   function3_curve_data(fdef, fdef2, fdef3, dim, npts);
 else random_curve_data(dim, npts);
 if (p <= 0) { printf("Dimnesion must be more/equal: 1\n"); exit(36); }
#ifdef DEBUG
 bcnt = 0;
#endif
 ctlpts = NULL;
 N = matrix(n+1);
 compute_nu();
 NI = invert_matrix(N, n+1);
#ifdef DEBUG
 print_vector(t, n+1, "t");
 print_vector(knot, n+p+1, "knot");
 print_matrix3(D, n+1, n+1, s, "D");
 print_matrix(N, n+1, "N");
 print_matrix(NI, n+1, "NI");
#endif
 construct_Q();
#ifdef DEBUG
 print_matrix3(Q, n+1, n+1, s, "Q");
#endif
 construct_P();
#ifdef DEBUG
 print_matrix3(P, n+1, n+1, s, "P");
#endif
 pts = NULL;
 fflush(stdout);
 if (deform) init_tricubic_engine();
 regenerate_pts();
 printf("done!\n");
 help();
 angY = angZ = 0.;
 angX = -180.;
 si = sj = 0;
}


void glInit()
{
 glClearColor(1.,1.,1.,0.);
 glShadeModel(GL_SMOOTH);
}


void nurbInit()
{
   GLfloat mat_diffuse[] = { 0.6, 0.6, .6, 1.0 };
   GLfloat mat_specular[] = { 1., 1.0, 1., 1.0 };
   GLfloat mat_shininess[] = { 127.0 };
   GLfloat mat_diffuse2[] = { 0.6, 0.6, 0.6, 1.0 };
   GLfloat mat_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat mat_shininess2[] = { 127.0 };
   GLfloat lightPos0[] = { 0., 0., -1., 0.};
   GLfloat lightPos1[] = { 0., 0., 1., 0.};
   GLfloat lightColor[] = { 1.0f, 1.0f, 1.0f, 1.0 };
   GLfloat lightColor1[] = { .7f, 1.0f, 0.0f, 1.0 };
   glClearColor(1., 1., 1., 0.);
   glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glMaterialfv(GL_BACK, GL_DIFFUSE, mat_diffuse2);
   glMaterialfv(GL_BACK, GL_SPECULAR, mat_specular2);
   glMaterialfv(GL_BACK, GL_SHININESS, mat_shininess2);
   glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
   glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor1);
   glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor1);
   glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor);
   glLightfv(GL_LIGHT1, GL_SPECULAR, lightColor);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_LIGHT1);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_AUTO_NORMAL);
   glEnable(GL_NORMALIZE);
   init_ctlpts();
   precision = 15.;
   theNurb = gluNewNurbsRenderer();
   gluNurbsProperty(theNurb, GLU_SAMPLING_TOLERANCE, precision);
   gluNurbsProperty(theNurb, GLU_DISPLAY_MODE, GLU_FILL);
}


int main(int lb, char** par)
{
 wx = 800;
 wy = 600;
 if (lb < 3)
   {
    help();
    exit(20);
   }
 if (par[1][0] == 'd' || par[1][0] == 'D') deform = 1;
 lb--;
 par++;
 if (par[1][1] == 'w') solid = 0;
 else if (par[1][1] == 'o') solid = 1;
 else exit(33);
 if (solid) deform = 0;
 if (par[1][0] == 'l') Init(par[2], NULL, NULL, NULL, -1, -1);
 else if (par[1][0] == 'f' && lb >= 5) Init(NULL, par[4], NULL, NULL, atoi(par[2]), atoi(par[3])-1);
 else if (par[1][0] == 'r' && lb >= 3) Init(NULL, NULL, NULL, NULL, atoi(par[2]), atoi(par[3])-1);
 else if (par[1][0] == 'u' && lb >= 7) Init(NULL, par[4], par[5], par[6], atoi(par[2]), atoi(par[3])-1);
 else exit(21);
 glutInit(&lb, par);
 glutInitDisplayMode(GLUT_DOUBLE);
 glutInitWindowSize(wx, wy);
 glutInitWindowPosition(100, 100);
 glutCreateWindow(par[0]);
 glutCreateMenu(menu);
 glutAddMenuEntry("Write Bspline surface", MENU_WRITE_BSP);
 glutAddMenuEntry("Toggle auto move", MENU_AUTO_MOVE);
 glutAddMenuEntry("Toggle display points D", MENU_D);
 glutAddMenuEntry("Toggle display points P", MENU_P);
 glutAddMenuEntry("Toggle display points Q", MENU_Q);
 glutAddMenuEntry("Zoom", MENU_SPLUS);
 glutAddMenuEntry("Unzoom", MENU_SMINUS);
 glutAddMenuEntry("Start postions", MENU_NORMAL);
 glutAddMenuEntry("More U-Lines", MENU_UPLUS);
 glutAddMenuEntry("Less U-Lines", MENU_UMINUS);
 glutAddMenuEntry("More V-Lines", MENU_VPLUS);
 glutAddMenuEntry("Less V-Lines", MENU_VMINUS);
 glutAddMenuEntry("Lower degree", MENU_LOWER_DEGREE);
 glutAddMenuEntry("Higher degree", MENU_HIGHER_DEGREE);
 glutAddMenuEntry("Quit", MENU_QUIT);
 glutAttachMenu(GLUT_RIGHT_BUTTON);
 glutAttachMenu(GLUT_LEFT_BUTTON);
 if (!solid)
   {
    glInit();
    glutDisplayFunc(render_scene);
   }
 else
   {
    nurbInit();
    glutDisplayFunc(render_scene_solid);
   }
 glutReshapeFunc(resize_scene);
 glutKeyboardFunc(keyboard);
 glutVisibilityFunc(visible);
 glutMainLoop();
 return 0;
}

/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
