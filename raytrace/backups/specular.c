#*/include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define HERE __FILE__,__LINE__
#define BACK_R 0x4F
#define BACK_G 0x1F
#define BACK_B 0x7F
#define MAX_REC 8
#define distorber (double)(0.5)
enum { RED = 0x0606, GREEN, BLUE};

typedef struct _Vertex
{
 double x,y,z;
} Vertex;
typedef struct _Vertex Vector;

typedef struct _Surface
{
 double A,B,C,D;
} Surface;

typedef struct _Ray
{
 Vertex P;
 Vector d;
 int r,c;
} Ray;

typedef struct _Triangle
{
 Vertex a,b,c;
 Vector na,nb,nc;	/* trojkat moze byc krzywoliniowy, ale */
 double tr,tg,tb;	/* zalecana mala krzywizna bo surface plaska */
 double sr,sg,sb;
 double cr,cg,cb;
 Surface s;		/* jakby trojkat byl plaski */
 double mU,mD;
 double ca;
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

typedef struct _VList
{
 struct _VList *next, *prev;
 Vertex P;
 int idx;
} VList;
Vertex light, observer;	/* FIXME: more lights ? */
Screen screen;
int nTriangles;
FILE* dbg;

void panic(char* f, int l)
{
 printf("RayTracer Engine Panic: file: %s, line: %d\n", f,l);
 exit(1);
}


void init_bmp(BMPTag* b)
{
 int i;
 if (!b) panic(HERE);
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


void init_screen(Screen* scr, int x, int y)
{
 if (!scr) panic(HERE);
 scr->pixels = (unsigned char*)malloc((x*y*3+1)*sizeof(unsigned char));
 if (!scr->pixels) panic(HERE);
 scr->x = x;
 scr->y = y;
}


void free_screen(Screen* scr)
{
 if (!scr) panic(HERE);
 free(scr->pixels);
 scr->x = scr->y = 0;
}


void print_ray(Ray* r)
{
 printf("FROM: (%f,%f,%f), DIR: (%f,%f,%f), REC(%d), COL(%d)\n",
	 r->P.x, r->P.y, r->P.z, r->d.x, r->d.y, r->d.z, r->r, r->c);
}


int read_vector(FILE* f, Vector* v)
{
 int nr;
 if (!f || !v) panic(HERE);
 nr = fscanf(f, "Vector: (%lf,%lf,%lf)", &v->x, &v->y, &v->z);
 if (nr == 3) return 1;
 else return 0;
}


int read_vertex(FILE* f, Vertex* v)
{
 int nr;
 if (!f || !v) panic(HERE);
 nr = fscanf(f, "Vertex: (%lf,%lf,%lf)", &v->x, &v->y, &v->z);
 if (nr == 3) return 1;
 else return 0;
}

void normalize(Vector*);

void normalize_color_factors(Triangle* t)
{
 double f;
 f = t->cr + t->sr + t->tr;
 t->cr /= f;
 t->sr /= f;
 t->tr /= f;
 f = t->cg + t->sg + t->tg;
 t->cg /= f;
 t->sg /= f;
 t->tg /= f;
 f = t->cb + t->sb + t->tb;
 t->cb /= f;
 t->sb /= f;
 t->tb /= f;
}


void  compute_normals(Triangle* t)
{
 /* when triangles get clockwise */
 double x1,x2,y1,y2,z1,z2,len;
 double nx,ny,nz;
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
 t->na.x = nx;
 t->nb.x = nx;
 t->nc.x = nx;
 t->na.y = ny;
 t->nb.y = ny;
 t->nc.y = ny;
 t->na.z = nz;
 t->nb.z = nz;
 t->nc.z = nz;
/* printf("normals: %f,%f,%f\n", t->na.x, t->na.y, t->na.z);*/
}


void  compute_surface(Triangle* t)
{
 double x1,x2,y1,y2,z1,z2,len;
 double nx,ny,nz;
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

double length(Vector* v);

int read_triangle(FILE* f, Triangle* t)
{
 int nr,n;
 if (!f || !t) panic(HERE);
 nr = fscanf(f, "Triangle: %d\n{\n", &n);
 if (nr != 1) return 0;
 if (n < 0) return 0;

 fscanf(f, " a: ");
 if (!read_vertex(f, &t->a)) panic(HERE);
 fscanf(f,"\n");
 t->a.x += distorber;
 
 fscanf(f, " b: ");
 if (!read_vertex(f, &t->b)) panic(HERE);
 fscanf(f,"\n");
 t->b.x += distorber;
 
 fscanf(f, " c: ");
 if (!read_vertex(f, &t->c)) panic(HERE);
 fscanf(f,"\n");
 t->c.x += distorber;
 
 fscanf(f, " na: ");
 if (!read_vector(f, &t->na)) panic(HERE);
 fscanf(f,"\n");
 fscanf(f, " nb: ");
 if (!read_vector(f, &t->nb)) panic(HERE);
 fscanf(f,"\n");
 fscanf(f, " nc: ");
 if (!read_vector(f, &t->nc)) panic(HERE);
 if (length(&t->na) <= 1e-9 || length(&t->nb) <= 1e-9 || length(&t->nc) <= 1e-9)
   {
    compute_normals(t);
   }
 else
   {
    normalize(&t->na);
    normalize(&t->nb);
    normalize(&t->nc);
   }
 fscanf(f,"\n");
 nr = fscanf(f," transparency: RGB: (%lf,%lf,%lf)\n", &t->tr, &t->tg, &t->tb);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," specular: RGB: (%lf,%lf,%lf)\n", &t->sr, &t->sg, &t->sb);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," diffuse: RGB: (%lf,%lf,%lf)\n", &t->cr, &t->cg, &t->cb);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," surface: ABCD: (%lf,%lf,%lf,%lf)\n", &t->s.A, &t->s.B, &t->s.C, &t->s.D);
 if (nr != 4) panic(HERE);
 normalize_color_factors(t);
 if (t->s.A == 0. && t->s.B == 0. && t->s.C == 0. && t->s.D == 0.) compute_surface(t);
 nr = fscanf(f, " transparencyFact: (%lf,%lf)\n", &t->mU, &t->mD);
 if (nr != 2) panic(HERE);
 nr = fscanf(f, " specularFact: %lf\n", &t->ca);
 if (nr != 1) panic(HERE);
 fscanf(f,"}\n");
 return 1;
}


void load_scene(Triangle** ts, char* scenefile)
{
 FILE* f;
 int nr,sx,sy,n,i;
 f = fopen(scenefile,"r");
 if (!f) panic(HERE);
 nr = fscanf(f, "Screen: (%d,%d)\n", &sx, &sy);
 if (nr != 2) panic(HERE);
 if (sx <= 0 || sy <= 0);
 init_screen(&screen, sx, sy);
 fscanf(f, "Observer: ");
 if (!read_vertex(f, &observer)) panic(HERE);
 fscanf(f,"\n");
 fscanf(f, "Light: ");
 if (!read_vertex(f, &light)) panic(HERE);
 fscanf(f,"\n");
 nr = fscanf(f, "nTriangles: %d\n", &n);
 if (nr != 1) panic(HERE);
 if (n <= 0) panic(HERE);
 *ts = (Triangle*)malloc(n*sizeof(Triangle));
 for (i=0;i<n;i++)
   {
    if (!read_triangle(f, &((*ts)[i]))) panic(HERE);
/*    printf("read %d triangles\n", i);*/
   }
 nTriangles = n;
 fclose(f);
 /*t = *ts;
 t->a.x = -200.;
 t->a.y = -200.;
 t->a.z =  100;
 t->b.x = -200.;
 t->b.y =  100.;
 t->b.z =  100.;
 t->c.x =  300.;
 t->c.y =  200.;
 t->c.z =  100.;
 t->na.x =  0.;
 t->na.y =  0.;
 t->na.z =  -1.;
 t->nb.x =  0.;
 t->nb.y =  0.;
 t->nb.z =  -1.;
 t->nc.x =  0.;
 t->nc.y =  0.;
 t->nc.z =  -1.;
 t->tr = .33;
 t->tg = .33;
 t->tb = .33;
 t->sr = .33;
 t->sg = .33;
 t->sb = .33;
 t->cr = .33;
 t->cg = .33;
 t->cb = .33;
 t->s.A = 0.;
 t->s.B = 0.;
 t->s.C = 1.;
 t->s.D = -100.;
 t->mU = 1.;
 t->mD = 1.8;
 t->ca = 100.;
 *ts = t;*/
}


void free_scene(Triangle** ts)
{
 free((void*)(*ts));
}


void set_color(Screen* s, int x, int y, int r, int g, int b)
{
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    s->pixels[3*(s->y * x + y)   ] = r;
    s->pixels[3*(s->y * x + y) +1] = g;
    s->pixels[3*(s->y * x + y) +2] = b;
   }
 else panic(HERE);
}


void color_background(Screen* s, int x, int y)
{
 set_color(s, x, y, BACK_R, BACK_G, BACK_B);
}


int nearly_equal(double x, double y, double eps)
{
 if (fabs(x-y) < eps) return 1;
 else return 0;
}


double compute_angle(Vertex* a, Vertex* b, Vertex* c)
{
 double x1,x2,y1,y2,z1,z2,len;
 x1 = b->x - a->x;
 x2 = c->x - a->x;
 y1 = b->y - a->y;
 y2 = c->y - a->y;
 z1 = b->z - a->z;
 z2 = c->z - a->z;
 len = sqrt(x1*x1+y1*y1+z1*z1);
 if (nearly_equal(len, 0., 1e-9)) return 4.;
 x1 /= len;
 y1 /= len;
 z1 /=len;
 len = sqrt(x2*x2+y2*y2+z2*z2);
 if (nearly_equal(len, 0., 1e-9)) return 4.;
 x2 /= len;
 y2 /= len;
 z2 /= len;
 return (x1*x2+y1*y2+z1*z2);
}


double compute_angle2(Vertex* a, Vertex* b, Vertex* c)
{
 double x1,x2,y1,y2,z1,z2,len;
 x1 = b->x - a->x;
 x2 = c->x - a->x;
 y1 = b->y - a->y;
 y2 = c->y - a->y;
 z1 = b->z - a->z;
 z2 = c->z - a->z;
 len = sqrt(x1*x1+y1*y1+z1*z1);
 x1 /= len;
 y1 /= len;
 z1 /= len;
 len = sqrt(x2*x2+y2*y2+z2*z2);
 x2 /= len;
 y2 /= len;
 z2 /= len;
 return acos(x1*x2+y1*y2+z1*z2);
}


int vertex_in_triangle(Triangle* t, Vertex* w)
{
 double kat;
 kat  = compute_angle(w, &t->a, &t->b);
 kat += compute_angle(w, &t->b, &t->c);
 kat += compute_angle(w, &t->c, &t->a);
/* printf("kat = %f\n", kat);*/
 if (kat < -0.999999)		/* FIXME: is this test OK? */
   {
    return 1;
   }
 else return 0;
}


void vlist_add(VList** head, Vertex* v, int idx)
{
 VList* temp;
 if (*head == NULL)
   {
    *head = (VList*)(malloc(sizeof(VList)));
    if (!(*head)) panic(HERE);
    (*head)->P.x = v->x;
    (*head)->P.y = v->y;
    (*head)->P.z = v->z;
    (*head)->idx = idx;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }
  temp = (VList*)(malloc(sizeof(VList)));
  if (!temp) panic(HERE);
  temp->P.x = v->x;
  temp->P.y = v->y;
  temp->P.z = v->z;
  temp->idx = idx;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
}


int intersection_triangle(Triangle* tr, Ray* r, Vertex* ret)
{
 double t, tmp;
 Vertex rett;
 Surface* s;
 if (!ret) panic(HERE);
 s = &tr->s;
 tmp = ( s->A*r->d.x +  s->B*r->d.y + s->C*r->d.z);
/* printf("ABCD=(%f,%f,%f,%f) tmp = %f\n", s->A, s->B, s->C, s->D, tmp);*/
 if (!nearly_equal(tmp, 0., 1e-9))
   {
    t = - ( s->A*r->P.x + s->B*r->P.y + s->C*r->P.z + s->D) / tmp;
/*    printf("t = %f\n", t);*/
    rett.x = r->d.x * t + r->P.x;
    rett.y = r->d.y * t + r->P.y;
    rett.z = r->d.z * t + r->P.z;
/*    printf("(%f,%f,%f)\n", rett.x, rett.y, rett.z);*/
    if (vertex_in_triangle(tr, &rett))
      {
/*	  printf("almost_inside\n");*/
        if (((rett.x-r->P.x)*r->d.x <= 1e-9) && ((rett.y-r->P.y)*r->d.y <= 1e-9) && ((rett.z-r->P.z)*r->d.z <= 1e-9)) return 0;
	ret->x = rett.x;
	ret->y = rett.y;
	ret->z = rett.z;
/*	  printf("check the same: (%f,%f,%f)\n", (rett.x-r->P.x)*r->d.x,(rett.y-r->P.y)*r->d.y,(rett.z-r->P.z)*r->d.z);*/
/*	printf("passed!\n");*/
	return 1;
      }
    else return 0;
   }
 return 0;
}

double distance(Vertex*, Vertex*);

void vlist_free(VList** head)
{
 VList* temp;
 if (!head) return;
/* printf("freeying\n");*/
 while (*head)
   {
    temp = *head;
/*    printf("%p -> %p\n", temp, temp->next);*/
    *head = (*head)->next;
    free(temp);
/*    printf("freed.\n");*/
   }
 *head = NULL;
}


int intersection(Triangle* tr, Ray* r, Vertex* ret, int* idx)
{
 int i;
 int sol;
 double dist,ndist;
 int ts;
 VList* head;
 VList* hd;
 VList* amb;
 VList* tamb;
 head = NULL;
 amb = NULL;
 sol = 0;
/* printf("TRACING:\n");*/
/* print_ray(r);*/
/* printf("ray> (%f,%f,%f) --> (%f,%f,%f)\n", r->P.x, r->P.y, r->P.z, r->d.x, r->d.y, r->d.z);*/
 for (i=0;i<nTriangles;i++)
   {
    if (intersection_triangle(&tr[i], r, ret))
      {
       vlist_add(&head, ret, i);
/*       printf("intersection with triangle: %d --> (%f,%f,%f)\n", i,ret->x, ret->y, ret->z);*/
       sol ++;
       *idx = i;
      }
   }
/*printf("intersections: %d\n", sol);*/
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", i,ret->x, ret->y, ret->z);*/
 if (sol <=1) { vlist_free(&head); return sol; }
/* printf("multiple soluions: searching nearest.\n");*/
 dist = 1e12;
 hd = head;
 while (head)
   {
    if ((ndist = distance(&r->P, &head->P)) < dist)
      {
/*       printf("new distance is: %f\n", ndist);*/
       dist = ndist;
       *idx = head->idx;
       ret->x = head->P.x;
       ret->y = head->P.y;
       ret->z = head->P.z;
      }
    head = head->next;
/*    printf("element.\n");*/
/*    printf("processed solution.\n");*/
   }
 head = hd;
 ts = 0;
 while (head)
   {
    if (nearly_equal(distance(&r->P, &head->P), dist, 1e-9)) 
      {
       ts++;
       vlist_add(&amb, &head->P, head->idx);
      }
    head = head->next;
   }
 if (ts > 1)
 {
/*  printf("Ambiguous solutions.\n");*/
     printf("*");
  printf("\nAbmigous solution for ray:\n");
  print_ray(r);
  tamb = amb;
  while (amb)
    {
     printf("With triangle: %d, solution: (%f,%f,%f,), distance: %f\n",
	     amb->idx, amb->P.x, amb->P.y, amb->P.z, distance(&r->P, &amb->P));
     amb = amb->next;
    }
  amb = tamb;
  vlist_free(&amb);
  vlist_free(&head);
  /* FIXME: must be handled somehow */ 
/*  panic(HERE);*/
  return 0;
 }
 /* FIXME: choose nearest solution */
 /* and free list and point to proper index */
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", *idx,ret->x, ret->y, ret->z);*/
/* printf("big_free\n");*/
 head = hd;
 vlist_free(&head);
 return sol;
}


int line_intersect(Vertex* p, Ray* l1, Ray* l2)
{
 double x1,y1,z1,x2,y2,z2;
 double dx1,dy1,dz1,dx2,dy2,dz2;
 double K,L;
 double xr,yr,zr;
/* printf("LINE_INTERSECT:\n");*/
/* print_ray(l1);*/
/* print_ray(l2);*/
/* printf("LINE_INTERSECT_ENDS\n");*/
 if (!p || !l1 || !l2) panic(HERE);
 x1 = l1->P.x;
 y1 = l1->P.y;
 z1 = l1->P.z;
 x2 = l2->P.x;
 y2 = l2->P.y;
 z2 = l2->P.z;
 dx1 = l1->d.x;
 dy1 = l1->d.y;
 dz1 = l1->d.z;
 dx2 = l2->d.x;
 dy2 = l2->d.y;
 dz2 = l2->d.z;
 if (nearly_equal(dx1, 0., 1e-9) &&
	 nearly_equal(dy1, 0., 1e-9) &&
	 nearly_equal(dz1, 0., 1e-9)) return -1;	/* degraded line */
 if (nearly_equal(dx2, 0., 1e-9) &&
	 nearly_equal(dy2, 0., 1e-9) &&
	 nearly_equal(dz2, 0., 1e-9)) return -1;	/* degraded line */
 if ( (nearly_equal(dx1, 0., 1e-9) &&
	 nearly_equal(dy1, 0., 1e-9)) ||
	 (nearly_equal(dx2, 0., 1e-9) &&
	 nearly_equal(dy2, 0., 1e-9)))
    {
     /* jedna z linii to  z = const */
     /* inny rzut*/
     if ( (nearly_equal(dx1, 0., 1e-9) &&
	 nearly_equal(dz1, 0., 1e-9)) ||
	 (nearly_equal(dx2, 0., 1e-9) &&
	 nearly_equal(dz2, 0., 1e-9)))
        {
         /* jedna z linii to y = const */
	 /* oblicz z rzutu na x = const bo napewno obie ok */
	 x_case:
         if (nearly_equal(dy1*dz2-dy2*dz1, 0., 1e-9))
          {
           /* rownoleglosc w x = const */
           if (!nearly_equal(dx1*dz2-dx2*dz1, 0., 1e-9)) goto y_case;
           if (!nearly_equal(dx1*dy2-dx2*dy1, 0., 1e-9)) goto z_case;
	   /* real parallel in 3d */
           goto p_case;
          }
         L = (dz1*(y2-y1)-dy1*(z2-z1))/(dy1*dz2-dy2*dz1);
         if (!nearly_equal(dy1, 0., 1e-9)) K = (y2-y1+L*dy2)/dy1;
         else if (!nearly_equal(dz1, 0., 1e-9)) K = (z2-z1+L*dz2)/dz1;
         if (!nearly_equal(x1+K*dx1, x2+L*dx2, 1e-9)) return 0; /* not 3d intersect */
         p->x = x1 + K*dx1;
         p->y = y1 + K*dy1;
         p->z = z1 + K*dz1;
         return 1;	/* intersection OK */
        }
     /* oblicz z rzutu na y = const */
      y_case:
      if (nearly_equal(dx1*dz2-dx2*dz1, 0., 1e-9))
       {
        /* rownoleglosc w y = const */
           if (!nearly_equal(dy1*dz2-dy2*dz1, 0., 1e-9)) goto x_case;
           if (!nearly_equal(dx1*dy2-dx2*dy1, 0., 1e-9)) goto z_case;
	   /* real parallel in 3d */
           goto p_case;
       }
     L = (dz1*(x2-x1)-dx1*(z2-z1))/(dx1*dz2-dx2*dz1);
     if (!nearly_equal(dx1, 0., 1e-9)) K = (x2-x1+L*dx2)/dx1;
     else if (!nearly_equal(dz1, 0., 1e-9)) K = (z2-z1+L*dz2)/dz1;
     if (!nearly_equal(y1+K*dy1, y2+L*dy2, 1e-9)) return 0; /* not 3d intersect */
     p->x = x1 + K*dx1;
     p->y = y1 + K*dy1;
     p->z = z1 + K*dz1;
     return 1;	/* intersection OK */
    }
 /* oblicz z rzutu na z = const */
 z_case:
 if (nearly_equal(dx1*dy2-dx2*dy1, 0., 1e-9))
   {
   /* rownoleglosc w z = const */
   if (!nearly_equal(dy1*dz2-dy2*dz1, 0., 1e-9)) goto x_case;
   if (!nearly_equal(dx1*dz2-dx2*dz1, 0., 1e-9)) goto y_case;
   /* real parallel in 3d */
   goto p_case;
   }
 L = (dy1*(x2-x1)-dx1*(y2-y1))/(dx1*dy2-dx2*dy1);
 if (!nearly_equal(dx1, 0., 1e-9)) K = (x2-x1+L*dx2)/dx1;
 else if (!nearly_equal(dy1, 0., 1e-9)) K = (y2-y1+L*dy2)/dy1;
 if (!nearly_equal(z1+K*dz1, z2+L*dz2, 1e-9)) return 0; /* not 3d intersect */
 p->x = x1 + K*dx1;
 p->y = y1 + K*dy1;
 p->z = z1 + K*dz1;
 return 1;	/* intersection OK */
 p_case:
 xr = (x2 - x1) / dx1;
 yr = (y2 - y1) / dy1;
 zr = (z2 - z1) / dz1;
 if (!nearly_equal(xr, yr, 1e-9)) return 0; /* parallel disjoint */
 if (!nearly_equal(xr, zr, 1e-9)) return 0; /* parallel disjoint */
 if (!nearly_equal(yr, zr, 1e-9)) return 0; /* parallel disjoint */
 p->x = x1;
 p->y = y1;
 p->z = z1;
 return 2;	/* parallel the same*/
}


void make_ray(Ray* l, Vertex* a, Vertex* b)
{
 if (!a || !b || !l) panic(HERE);
 l->P.x = a->x;
 l->P.y = a->y;
 l->P.z = a->z;
 l->d.x = b->x - a->x;
 l->d.y = b->y - a->y;
 l->d.z = b->z - a->z;
 l->r = 0;
 l->c = RED;
}


double distance(Vertex* a, Vertex* b)
{
 return sqrt((b->x-a->x)*(b->x-a->x)+(b->y-a->y)*(b->y-a->y)+(b->z-a->z)*(b->z-a->z));
}

void normalize(Vector* v)
{
 double len;
 len = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
 if (nearly_equal(len, 0., 1e-9)) panic(HERE);
 v->x /= len;
 v->y /= len;
 v->z /= len;
}


void get_normal(Vector* n, Triangle* t, Vertex* p)
{
 Ray bc;
 Ray ap;
 Vertex d;
 Vector tn;
 double lenBC, lenBD;
 double fact;
 make_ray(&bc, &t->b, &t->c);
 make_ray(&ap, &t->a, p);
 /*printf("P = (%f,%f,%f)\n", p->x, p->y, p->z);
 print_ray(&bc);
 print_ray(&ap);*/;
 if (line_intersect(&d, &bc, &ap) != 1) panic(HERE);
 if (nearly_equal(t->a.x, d.x, 1e-9) &&
	 nearly_equal(t->a.y, d.y, 1e-9) &&
	 nearly_equal(t->a.y, d.y, 1e-9))
   {
    n->x = t->na.x;
    n->y = t->na.y;
    n->z = t->na.z;
    return ;
   }
 lenBC = distance(&t->c, &t->b);
 if (nearly_equal(lenBC, 0., 1e-9)) panic(HERE);
 lenBD = distance(&d, &t->b);
 fact = lenBD / lenBC;
 tn.x = (1.-fact)* t->nb.x + fact*t->nc.x;
 tn.y = (1.-fact)* t->nb.y + fact*t->nc.y;
 tn.z = (1.-fact)* t->nb.z + fact*t->nc.z;
 lenBC = distance(&d, &t->a);
 if (nearly_equal(lenBC, 0., 1e-9)) panic(HERE);
 lenBD = distance(p, &t->a);
 fact = lenBD / lenBC;
 n->x = (1.-fact)* t->na.x + fact*tn.x;
 n->y = (1.-fact)* t->na.y + fact*tn.y;
 n->z = (1.-fact)* t->na.z + fact*tn.z;
 normalize(n);
/* printf("gn: (%f,%f,%f)\n",p->x,p->y,p->z);*/
/* printf("intersect: (%3.0f,%3.0f,%3.0f),( %3.0f,%3.0f,%3.0f), %3.0f,%3.0f; (%1.2f,%1.2f,%1.2f)\n", p->x, p->y, p->z, d.x, d.y, d.z, lenBC, lenBD, n->x, n->y, n->z);*/
}


double scalar_prod(Vector* v, Vector* w)
{
 return v->x*w->x+v->y*w->y+v->z*w->z;
}


double length(Vector* v)
{
 return sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
}


double recurse_color(Ray* r, Triangle* t, Vertex* v, Vector* n, double f)
{
 double rv;
 double ca;
 double fact, fact2;
 double a2;
 double nf;
 double add;
 int idx;
 Vector obsV, nn, tr;
 Vertex nv;
 Ray vi;
/* print_ray(r);*/
 if (r->r >= MAX_REC) return 0.;
/* printf("intersect: (%f,%f,%f)\n", v->x, v->y, v->z);*/
/* printf("N: (%f,%f,%f)\n", n->x, n->y, n->z);*/
 rv = 0.;
 tr.x = light.x - v->x;
 tr.y = light.y - v->y;
 tr.z = light.z - v->z;
/* printf("L: (%f,%f,%f)\n", tr.x, tr.y, tr.z);*/
 obsV.x = r->P.x - v->x;
 obsV.y = r->P.y - v->y;
 obsV.z = r->P.z - v->z;
/* printf("V: (%f,%f,%f)\n", obsV.x, obsV.y, obsV.z);*/
 if (length(&obsV) < 1e-9) return 0.;
 if (length(&tr) < 1e-9) return 0.;
 normalize(&tr);
 normalize(&obsV);
/* printf("nL: (%f,%f,%f)\n", tr.x, tr.y, tr.z);*/
/* printf("nV: (%f,%f,%f)\n", obsV.x, obsV.y, obsV.z);*/
 ca = scalar_prod(&tr, n);
 if (ca < 0.) ca = 0.;
/* printf("ca = %f\n", ca);*/
 if (r->c == RED)   fact = t->cr;
 if (r->c == GREEN) fact = t->cg;
 if (r->c == BLUE)  fact = t->cb;
/* printf("fact = %f\n", fact);*/
 a2 = 2.*scalar_prod(&obsV, n);
/* printf("a2 = %f\n", a2);*/
 vi.P.x = v->x;
 vi.P.y = v->y;
 vi.P.z = v->z;
 vi.d.x = a2*n->x - obsV.x;
 vi.d.y = a2*n->y - obsV.y;
 vi.d.z = a2*n->z - obsV.z;
 vi.r = r->r + 1;
 vi.c = r->c;
/* print_ray(&vi);*/
 normalize(&vi.d);
 if (ca < 0.) ca = 0.;
 if (ca > 1.) ca = 1.;
/* printf("ca = %f\n", ca);*/
/* printf("diffuse adds: %f\n", fact * ca);*/
 rv += fact * ca;
 /* FIXME: is this specular flash OK? */
 rv += pow(scalar_prod(&vi.d, &tr), t->ca);
 if (rv > 1.) rv = 1.;
/* printf("rv = %f\n", rv);*/
 if (r->c == RED)   fact = t->sr /*+ t->tr*/;
 if (r->c == GREEN) fact = t->sg /*+ t->tg*/;
 if (r->c == BLUE)  fact = t->sb /*+ t->tb*/;
 nf = f * fact;
 /*if (rv < 0.01) 
   {
     printf("rv = %f, fact = %f, ca = %f\n", rv, fact, ca);
     printf("r> (%f,%f,%f) -> (%f,%f,%f)\n", r->P.x, r->P.y, r->P.z, r->d.x, r->d.y, r->d.z);
     printf("v> (%f,%f,%f)\n", v->x, v->y, v->z);
     printf("t> (%f,%f,%f) (%f,%f,%f) (%f,%f,%f)\n", t->a.x, t->a.y, t->a.z, t->b.x, t->b.y, t->b.z, t->c.x, t->c.y, t->c.z);
   }*/
/* printf("nf = %f\n", nf);*/
/* rv = 0.;*/
 if (nf > 1./256. && (rv <1.-1/256.))
   {
/*       return rv;*/
    /* FIXME: continue here! */
    /*oblicz odbicie...*/
    if (!intersection(t, &vi, &nv, &idx))
       {
        if (r->c == RED)   fact2 = BACK_R;
        if (r->c == GREEN) fact2 = BACK_G;
        if (r->c == BLUE)  fact2 = BACK_B;
	rv += nf * (fact2/255.);
/*	rv = 1.;*/
/*	printf("rec_level = %d, f = %f, col_bkgnd = %f\n", r->r, f, fact*(fact2/255.));*/
/*	fprintf(dbg,"rec_level = %d, f = %f, col_bkgnd = %f\n", r->r, f, fact*(fact2/255.));*/
/*	printf("background.\n");*/
/*	printf("no intersection --> background\n");*/
       }
    else
       {
/*	   printf("specular computation.\n");*/
/*       printf("intersect: (%d,%d)\n", x,y);*/
         get_normal(&nn, &t[idx], &nv);
	 add = nf * recurse_color(&vi, t, &nv, &nn, nf);
	 rv += add;
       }
   }
/* else printf("skipped!\n");*/
 if (rv <= 1./256.) return 0.;
 if (rv > 1.) rv = 1.;
/* printf("FINAL rv = %f\n", rv);*/
 return rv;
}


void calculate_color(Screen* s, Ray* r, Triangle* t, int x, int y)
{
 Vertex v;
 Vector n;
 double rr,gg,bb;
 int idx;
 /* moze byc wiele przeciec z obiektami, posortowac w/g kolejnosci */
 if (!intersection(t, r, &v, &idx))
   {
/*     printf("-");*/
/*       printf("error with:\n");*/
/*       print_ray(r);*/
       fprintf(dbg,"no intr at: (%d,%d)\n", x,y);
       fprintf(dbg, "(%f,%f,%f) --> (%f,%f,%f)\n", r->P.x, r->P.y, r->P.z, r->d.x, r->d.y, r->d.z);
       color_background(s, x, y);
/*       set_color(s, x, y, 0,0,0);*/
   }
 else
   {
/*       printf("intersect: (%d,%d)\n", x,y);*/
/*       printf("V = (%f,%f,%f)\n", v.x, v.y, v.z);*/
       get_normal(&n, &t[idx], &v);
       r->c = RED;
       rr = recurse_color(r, t, &v, &n, 1.);
       r->c = GREEN;
       gg = recurse_color(r, t, &v, &n, 1.);
       r->c = BLUE;
       bb = recurse_color(r, t, &v, &n, 1.);
/*       if (rr > .99)*/
/*       fprintf(dbg,"COMPUTED: RGB: (%f,%f,%f) at (%d,%d)\n", rr,gg,bb,x,y);*/
       set_color(s, x, y, (int)(rr*255.), (int)(gg*255.), (int)(bb*255.));
/*       panic(HERE);*/		/* BIG_PANIC */
/*       printf("*");*/
   }
}


void raytrace(Screen* s, char* scenefile)
{
 Ray r;
 Triangle* ts;
 int i,j;
 dbg = fopen("debug.out", "w");
 if (!dbg) panic(HERE);
 load_scene(&ts, scenefile);
 r.P.x = observer.x;
 r.P.y = observer.y;
 r.P.z = observer.z;;
 r.r = 0;
 r.d.z = (s->x+s->y)/3.;	/* FIXME: mozna zmieniac aby dokladniej zobaczyc */
 printf("\n");
 /*calculate_color(s, &r, ts, 400, 300);*/
 /*printf("TEST1\n");
 r.d.x = 40 - s->x/2.;
 r.d.y = 33 - s->y/2.;
 calculate_color(s, &r, ts, 40, 33);
 printf("TEST2\n");
 r.d.x = 40 - s->x/2.;
 r.d.y = 34 - s->y/2.;
 calculate_color(s, &r, ts, 40, 34);
 panic(HERE);*/
 for (i=0;i<s->x;i++)
   {
    r.d.x = i - s->x/2.;
    printf("%04d/%04d ", i+1,s->x);
    for (j=0;j<s->y;j++)
      {
       if (!(j % 32)) { printf("."); fflush(stdout); }
       r.d.y = j - s->y/2.;
       calculate_color(s, &r, ts, i, j);  
      }
    printf(" (%3.3f%%)\n", ((double)i*100.)/(double)s->x);
/*    printf("\n");*/
   }
 printf("\n");
 free_scene(&ts);
 fclose(dbg);
}


int get_red(Screen* s, int x, int y)
{
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)];
   }
 else panic(HERE);
 return -1;
}


int get_green(Screen* s, int x, int y)
{
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)+1];
   }
 else panic(HERE);
 return -1;
}


int get_blue(Screen* s, int x, int y)
{
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)+2];
   }
 else panic(HERE);
 return -1;
}


void wrt_bmp(Screen* s, char* out_f)
{
 BMPTag bm_handle;
 FILE* plik;
 int i,j;
 init_bmp(&bm_handle);
 plik = fopen(out_f,"w");
 if (!plik) { printf("Error writing BMP: %s\n", out_f); panic(HERE); }
 fprintf(plik,"%c%c",'B', 'M');
 bm_handle.bm_y = s->y;
 bm_handle.bm_x = s->x;
 bm_handle.fsize = sizeof(BMPTag)+(bm_handle.bm_y*bm_handle.bm_x*3);
 fwrite(&bm_handle.fsize,4,1,plik);
 fwrite(&bm_handle.dummy,4,1,plik);
 bm_handle.offset=sizeof(BMPTag);
 bm_handle.planes=1;
 bm_handle.bpp=24;
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
 for (i=0;i<bm_handle.bm_y;i++)  for (j=0;j<bm_handle.bm_x;j++)
   fprintf(plik,"%c%c%c", get_blue(s, j, i), get_green(s, j, i), get_red(s, j, i));
 fclose(plik);
}


int main(int lb, char** par)
{
 if (lb < 2) { printf("usage: %s scenefile\n",par[0]); return 1; }
 raytrace(&screen, par[1]);
 wrt_bmp(&screen, "screen.bmp");
 free_screen(&screen);
 return 0;
}

