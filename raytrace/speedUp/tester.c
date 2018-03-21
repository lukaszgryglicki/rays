#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>

#define REAL long double
int __times;
time_t __t_start, __t_cur;

enum { RED = 0x0606, GREEN, BLUE};

typedef struct _Vertex
{
 REAL x,y,z;
} Vertex;
typedef struct _Vertex Vector;

typedef struct _Material
{
 REAL c,s,t;
} Material;

typedef struct _TexCoord
{
 REAL x,y;
} TexCoord;

typedef struct _Surface
{
 REAL A,B,C,D;
} Surface;

typedef struct _Ray
{
 Vertex P;
 Vector d;
 int r,c;
} Ray;
typedef struct _Screen Texture;

typedef struct _NURBS
{
 int idx;
 int ntri;
 int n1,n2;
 int p1,p2;
 int m1,m2;
 int d1,d2;
 int itex;
 int faces, tid;
 REAL*** D;
 REAL*** values;	/* temporary */
 REAL*** normals;	/* temporary */
 REAL** w;
 REAL* knot1;
 REAL* knot2;
 REAL* t1;		/* temporary */
 REAL* t2;		/* temporary */
} NURBS;

typedef struct _Triangle
{
 Vertex a,b,c;
 Vector na,nb,nc;	/* trojkat moze byc krzywoliniowy, ale */
 Material mra, mrb, mrc;
 Material mga, mgb, mgc;
 Material mba, mbb, mbc;
 TexCoord ta, tb, tc;
 Texture* tex;
 Surface s;		/* jakby trojkat byl plaski */
 REAL mUR,mDR;
 REAL mUG,mDG;
 REAL mUB,mDB;
 REAL ca;
 REAL n_dist;
 int faces;
 int idx;
 int tid;
 int nidx;
 NURBS* nurbs;
 TexCoord nca;
 TexCoord ncb;
 TexCoord ncc;
} Triangle;

typedef struct _Screen
{
 unsigned char* pixels;
 int x,y;
} Screen;

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

typedef struct _Box
{
 REAL minx,miny,minz;
 REAL maxx,maxy,maxz;
 Triangle *t1, *t2;
} Box;

typedef struct _BTree
{
 struct _BTree *l, *r;
 Box *b;
} BTree;

typedef struct _VList
{
 struct _VList *next, *prev;
 Vertex P;
 int idx;
} VList;

typedef struct _PtrList
{
 struct _PtrList *next, *prev;
 void* ptr;
} PtrList;

typedef struct _ListTransform
{
 REAL ** M;
 REAL ** MN;
 int i1, i2;
 int ni1, ni2;
 REAL n_dist;
 int t_id;
} ListTransform;

typedef struct _BList
{
 struct _BList *next, *prev;
 Box b;
} BList;

int nearly_equal(REAL x, REAL y, REAL eps)
{
 return (fabs(x-y) < eps);
}

REAL compute_angle(Vertex* a, Vertex* b, Vertex* c)
{
 register REAL x1,x2,y1,y2,z1,z2,len;
 x1 = b->x - a->x;
 x2 = c->x - a->x;
 y1 = b->y - a->y;
 y2 = c->y - a->y;
 z1 = b->z - a->z;
 z2 = c->z - a->z;
 len = sqrt(x1*x1+y1*y1+z1*z1);
 if (len < 1e-7) return 4.;
/* if (nearly_equal(len, 0., 1e-7)) return 4.;*/
 x1 /= len;
 y1 /= len;
 z1 /=len;
 len = sqrt(x2*x2+y2*y2+z2*z2);
/* if (nearly_equal(len, 0., 1e-7)) return 4.;*/
 if (len < 1e-7) return 4.;
 x2 /= len;
 y2 /= len;
 z2 /= len;
 return (x1*x2+y1*y2+z1*z2);
}

int vertex_in_triangle(Triangle* t, Vertex* w)
{
 REAL kat;
 kat  = compute_angle(w, &t->a, &t->b);
 kat += compute_angle(w, &t->b, &t->c);
 kat += compute_angle(w, &t->c, &t->a);
/* printf("kat = %Lf\n", kat);*/
 if (kat < -0.99999)
   {
    return 1;
   }
 else return 0;
}

int intersection_triangle(Triangle* tr, Ray* r, Vertex* ret)
{
 REAL t, tmp;
 Vertex rett;
 Surface* s;
 s = &tr->s;
 tmp = ( s->A*r->d.x +  s->B*r->d.y + s->C*r->d.z);
 if (!nearly_equal(tmp, 0., 1e-8))
   {
    t = - ( s->A*r->P.x + s->B*r->P.y + s->C*r->P.z + s->D) / tmp;
    rett.x = r->d.x * t + r->P.x;
    rett.y = r->d.y * t + r->P.y;
    rett.z = r->d.z * t + r->P.z;
    if (vertex_in_triangle(tr, &rett))
      {
        if (((rett.x-r->P.x)*r->d.x <= 1e-9) && ((rett.y-r->P.y)*r->d.y <= 1e-9) && ((rett.z-r->P.z)*r->d.z <= 1e-9)) return 0;
	ret->x = rett.x;
	ret->y = rett.y;
	ret->z = rett.z;
	return 1;
      }
    else return 0;
   }
 else return 0;
}


void catch_signal(int signo)
{
 time(&__t_cur);
 __t_cur -= __t_start;
 printf("Executed %d times in %d seconds, %f/s\n", 
	 __times, __t_cur, (float)__times/(float)__t_cur);
 exit(0);
}

void init()
{
 signal(SIGINT, catch_signal);
 signal(SIGALRM, catch_signal);
}

REAL pseudo_length(Vector* v)
{
 return v->x*v->x+v->y*v->y+v->z*v->z;
}

REAL scalar_prod(Vector* v, Vector* w)
{
 return v->x*w->x+v->y*w->y+v->z*w->z;
}

void normalize(Vector* v)
{
 register REAL len;
 len = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
 v->x /= len;
 v->y /= len;
 v->z /= len;
}

REAL calc_surface(Triangle* t)
{
 REAL alfa;
 REAL la,lb;
 Vector a,b;
 a.x = t->c.x - t->a.x;
 a.y = t->c.y - t->a.y;
 a.z = t->c.z - t->a.z;
 b.x = t->b.x - t->a.x;
 b.y = t->b.y - t->a.y;
 b.z = t->b.z - t->a.z;
 la = pseudo_length(&a);
 lb = pseudo_length(&b);
 if (la < 3e-4 || lb < 3e-4)
   {
    printf("Warning: 0 length edge\n");
    return 0.;
   }
 normalize(&a);
 normalize(&b);
 alfa = acos(scalar_prod(&a, &b));
 return sin(alfa) * la * lb;
}

void  compute_surface(Triangle* t)
{
 REAL x1,x2,y1,y2,z1,z2,len;
 REAL nx,ny,nz;
 x1 = t->b.x - t->a.x;
 x2 = t->c.x - t->a.x;
 y1 = t->b.y - t->a.y;
 y2 = t->c.y - t->a.y;
 z1 = t->b.z - t->a.z;
 z2 = t->c.z - t->a.z;
 nx = y1*z2 - z1*y2;
 ny = z1*x2 - x1*z2;
 nz = x1*y2 - y1*x2;
 len = sqrt(nx*nx + ny*ny + nz*nz);
 nx /= len;
 ny /= len;
 nz /= len;
 t->s.A = nx;
 t->s.B = ny;
 t->s.C = nz;
 t->s.D = -t->s.A*t->a.x - t->s.B*t->a.y - t->s.C*t->a.z;
/* printf("surface: %f,%f,%f,%f\n", t->s.A, t->s.B, t->s.C, t->s.D);*/
}

void test_engine()
{
 Triangle t;
 Vertex w;
 Ray ray;
 t.a.x = 0.;
 t.a.y = 0.;
 t.a.z = 0.;
 t.b.x = 100.;
 t.b.y = 0.;
 t.b.z = -50.;
 t.c.x = 50.;
 t.c.y = 80.;
 t.c.z = 50.;
 compute_surface(&t);
 ray.P.x = 50.;
 ray.P.y = 30.;
 ray.P.z = - 70.;
 ray.d.x = 1.;
 ray.d.y = 2.;
 ray.d.z = 40.;
 ray.r = 1;
 ray.c = RED;
 time(&__t_start);
 __times = 0;
 while (1)
   {
    intersection_triangle(&t, &ray, &w);
    __times ++;
    if (!(__times % 10000))
       { 
	 time(&__t_cur);
	 if (__t_cur >= __t_start + 10) goto fine;
       }
   }
fine:
 __t_cur -= __t_start;
 printf("Executed %d times in %d seconds, %f/s\n", 
	 __times, __t_cur, (float)__times/(float)__t_cur);
}

int main()
{
 init();
 test_engine();
 return 0;
}

