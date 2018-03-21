#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265
#define REAL double

enum { RED = 0x0606, GREEN, BLUE};

typedef struct _Screen
{
 unsigned char* pixels;
 int x,y,z;
} Screen;

typedef struct _Vertex
{
 REAL x,y,z,v;
} Vertex;

typedef struct _Vertex Vector;

typedef struct _Material
{
 REAL c,s,t;
} Material;

typedef struct _TexCoord
{
 REAL x,y,z;
} TexCoord;

typedef struct _Surface
{
 REAL A,B,C,D,E;
} Surface;

typedef struct _Ray
{
 Vertex P;
 Vector d;
 int r,c;
} Ray;

typedef struct _Screen Texture;

typedef struct _Triangle
{
 Vertex a,b,c,d;
 Vector na,nb,nc,nd;	
 Material mra, mrb, mrc, mrd;
 Material mga, mgb, mgc, mgd;
 Material mba, mbb, mbc, mbd;
 TexCoord ta, tb, tc, td;
 Texture* tex;
 Surface s;		
 REAL mUR,mDR;
 REAL mUG,mDG;
 REAL mUB,mDB;
 REAL caR, caG, caB;
 int idx;
 int tid;
} Triangle;

typedef struct _Box
{
 REAL minx,miny,minz,minv;
 REAL maxx,maxy,maxz,maxv;
 Triangle *t1, *t2;
} Box;

typedef struct _BTree
{
 struct _BTree *l, *r;
 Box *b;
} BTree;

typedef struct _BList
{
 struct _BList *next, *prev;
 Box b;
} BList;

typedef struct _IdxTree
{
 struct _IdxTree *s, *h, *t;
 int idx;
} IdxTree;

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

typedef struct _Voxel
{
 PtrList* tris;
 Box box;
 int ntbox;
 struct _Voxel* varr;
} Voxel;

#define HERE __FILE__,__LINE__



///////////////////////////////

Triangle* g_ts = NULL;
PtrList* pbv_head = NULL;
Screen screen;
Vertex observer, light;
REAL precision = 1e-6;
BTree* btree = NULL;
BList* boxes = NULL;
int nTriangles, bkgnd, max_rec;
REAL slimit_min = 10.;
REAL slimit_max = 1000.;
REAL slimit_step = 1.;
int g_idx = 0;
Texture* texture = NULL;
int nTex = 1;
IdxTree* itree = NULL;
REAL lookz = 0.25;

//////////////////////////////

void panic(char* fname, int fline)
{
	printf("panic at: %s, line %d\n", fname, fline);
	exit(0xFFFF);
}

void ptrlist_add(PtrList** head, void* t)
{
 PtrList* temp;
 if (*head == NULL)
   {
    *head = (PtrList*)(malloc(sizeof(PtrList)));

    (*head)->ptr = (void*)t;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }
  temp = (PtrList*)(malloc(sizeof(PtrList)));
  temp->ptr = (void*)t;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
}

void ptrlist_delete(PtrList** head, PtrList* ptr)
{
 PtrList* temp;
 temp = *head;
 while (temp)
   {
    if (temp == ptr)
      {
       if (temp == *head) *head = temp->next;
       if (temp->prev) temp->prev->next = temp->next;
       if (temp->next) temp->next->prev = temp->prev;
       free(temp);
       return;
      }
    temp = temp->next;
   }
}

void ptrlist_free(PtrList** head)
{
 PtrList* temp;
 if (!head) return;
 while (*head)
   {
    temp = *head;
    *head = (*head)->next;
    free(temp);
   }
}

PtrList* ptrlist_get_tail(PtrList* head)
{
 if (!head) return NULL;
 while (head->next)
   {
    head = head->next;
   }
 return head;
}

int count_items(PtrList* head)
{
 int item;
 item = 0;
 while (head) { item++; head = head->next; }
 return item;
}

IdxTree* alloc_itree()
{
 IdxTree* it;
 it = (IdxTree*)malloc(sizeof(IdxTree));
 it->idx = -1;
 it->s = NULL;
 it->h = NULL;
 it->t = NULL;
 return it;
}

REAL scalar_prod(Vector* v, Vector* w)
{
 return v->x*w->x+v->y*w->y+v->z*w->z+v->v*w->v;
}


REAL length(Vector* v)
{
 return sqrt(v->x*v->x+v->y*v->y+v->z*v->z+v->v*v->v);
}


REAL pseudo_length(Vector* v)
{
 return v->x*v->x+v->y*v->y+v->z*v->z+v->v*v->v;
}

void normalize(Vector* v)
{
 register REAL len;
 len = sqrt(v->x*v->x+v->y*v->y+v->z*v->z+v->v*v->v);
 v->x /= len;
 v->y /= len;
 v->z /= len;
 v->v /= len;
}

void normalize_color_factors(Triangle* t)
{
 REAL f;
 f = t->mra.c + t->mra.s + t->mra.t;
 if (f > 1.)
 {
 t->mra.c /= f;
 t->mra.s /= f;
 t->mra.t /= f;
 }
 f = t->mga.c + t->mga.s + t->mga.t;
 if (f > 1.)
 {
 t->mga.c /= f;
 t->mga.s /= f;
 t->mga.t /= f;
 }
 f = t->mba.c + t->mba.s + t->mba.t;
 if (f > 1.)
 {
 t->mba.c /= f;
 t->mba.s /= f;
 t->mba.t /= f;
 }
 f = t->mrb.c + t->mrb.s + t->mrb.t;
 if (f > 1.)
 {
 t->mrb.c /= f;
 t->mrb.s /= f;
 t->mrb.t /= f;
 }
 f = t->mgb.c + t->mgb.s + t->mgb.t;
 if (f > 1.)
 {
 t->mgb.c /= f;
 t->mgb.s /= f;
 t->mgb.t /= f;
 }
 f = t->mbb.c + t->mbb.s + t->mbb.t;
 if (f > 1.)
 {
 t->mbb.c /= f;
 t->mbb.s /= f;
 t->mbb.t /= f;
 }
 f = t->mrc.c + t->mrc.s + t->mrc.t;
 if (f > 1.)
 {
 t->mrc.c /= f;
 t->mrc.s /= f;
 t->mrc.t /= f;
 }
 f = t->mgc.c + t->mgc.s + t->mgc.t;
 if (f > 1.)
 {
 t->mgc.c /= f;
 t->mgc.s /= f;
 t->mgc.t /= f;
 }
 f = t->mbc.c + t->mbc.s + t->mbc.t;
 if (f > 1.)
 {
 t->mbc.c /= f;
 t->mbc.s /= f;
 t->mbc.t /= f;
 }
 f = t->mrd.c + t->mrd.s + t->mrd.t;
 if (f > 1.)
 {
 t->mrd.c /= f;
 t->mrd.s /= f;
 t->mrd.t /= f;
 }
 f = t->mgd.c + t->mgd.s + t->mgd.t;
 if (f > 1.)
 {
 t->mgd.c /= f;
 t->mgd.s /= f;
 t->mgd.t /= f;
 }
 f = t->mbd.c + t->mbd.s + t->mbd.t;
 if (f > 1.)
 {
 t->mbd.c /= f;
 t->mbd.s /= f;
 t->mbd.t /= f;
 }
}

void  compute_normals(Triangle* t)
{
 REAL ax,ay,az,av;
 REAL bx,by,bz,bv;
 REAL cx,cy,cz,cv;
 REAL nx,ny,nz,nv;
 REAL A, B, C, D, E, F, len;     

 // Vectors: a: AB, b: AC, c: AD
 ax = t->b.x - t->a.x;
 bx = t->c.x - t->a.x;
 cx = t->d.x - t->a.x;

 ay = t->b.y - t->a.y;
 by = t->c.y - t->a.y;
 cy = t->d.y - t->a.y;

 az = t->b.z - t->a.z;
 bz = t->c.z - t->a.z;
 cz = t->d.z - t->a.z;

 av = t->b.v - t->a.v;
 bv = t->c.v - t->a.v;
 cv = t->d.v - t->a.v;

 // Normal is a 4D orthogonal vector to 3 4D vectors in 4D space
 // Intermediate values A..F
 A = (bx * cy) - (by * cx);
 B = (bx * cz) - (bz * cx);
 C = (bx * cv) - (bv * cx);
 D = (by * cz) - (bz * cy);
 E = (by * cv) - (bv * cy);
 F = (bz * cv) - (bv * cz);

 // Final normal coords
 nx =   (ay * F) - (az * E) + (av * D);
 ny = - (ax * F) + (az * C) - (av * B);
 nz =   (ax * E) - (ay * C) + (av * A);
 nv = - (ax * D) + (ay * B) - (az * A);

 // Normalize part
 len = sqrt(nx*nx + ny*ny + nz*nz + nv*nv);

 nx /= len;
 ny /= len;
 nz /= len;
 nv /= len;

 // All 4 vertices has the same normal 
 // (hypertriangle is hyperflat thus lying within 3D hypersurface in 4D)
 t->na.x = nx;
 t->nb.x = nx;
 t->nc.x = nx;
 t->nd.x = nx;
 t->na.y = ny;
 t->nb.y = ny;
 t->nc.y = ny;
 t->nd.y = ny;
 t->na.z = nz;
 t->nb.z = nz;
 t->nc.z = nz;
 t->nd.z = nz;
 t->na.v = nv;
 t->nb.v = nv;
 t->nc.v = nv;
 t->nd.v = nv;
}

void  compute_surface(Triangle* t)
{
 // Computing 4D hypertriange's 3d hypersurface which its lying within
 REAL ax,ay,az,av;
 REAL bx,by,bz,bv;
 REAL cx,cy,cz,cv;
 REAL nx,ny,nz,nv;
 REAL A, B, C, D, E, F, len;     

 ax = t->b.x - t->a.x;
 bx = t->c.x - t->a.x;
 cx = t->d.x - t->a.x;

 ay = t->b.y - t->a.y;
 by = t->c.y - t->a.y;
 cy = t->d.y - t->a.y;

 az = t->b.z - t->a.z;
 bz = t->c.z - t->a.z;
 cz = t->d.z - t->a.z;

 av = t->b.v - t->a.v;
 bv = t->c.v - t->a.v;
 cv = t->d.v - t->a.v;

 A = (bx * cy) - (by * cx);
 B = (bx * cz) - (bz * cx);
 C = (bx * cv) - (bv * cx);
 D = (by * cz) - (bz * cy);
 E = (by * cv) - (bv * cy);
 F = (bz * cv) - (bv * cz);

 nx =   (ay * F) - (az * E) + (av * D);
 ny = - (ax * F) + (az * C) - (av * B);
 nz =   (ax * E) - (ay * C) + (av * A);
 nv = - (ax * D) + (ay * B) - (az * A);

 len = sqrt(nx*nx + ny*ny + nz*nz + nv*nv);

 nx /= len;
 ny /= len;
 nz /= len;
 nv /= len;

 t->s.A = nx;
 t->s.B = ny;
 t->s.C = nz;
 t->s.D = nv;

 t->s.E = -t->s.A*t->a.x - t->s.B*t->a.y - t->s.C*t->a.z - t->s.D*t->a.v;
}


void init_screen(Screen* scr, int x, int y, int z)
{
 int i,n;
 if (!scr) panic(HERE);
 scr->pixels = (unsigned char*)malloc((x*y*z*3+1)*sizeof(unsigned char));
 if (!scr->pixels) panic(HERE);
 scr->x = x;
 scr->y = y;
 scr->z = z;
 n = x * y * z * 3;
 for (i=0;i<=n;i++) scr->pixels[i] = 0;
}

int read_vertex(FILE* f, Vertex* v)
{
 int nr;
 if (!f || !v) panic(HERE);
 nr = fscanf(f, "(%lf,%lf,%lf,%lf)", &v->x, &v->y, &v->z, &v->v);
 if (nr == 4) return 1;
 else return 0;
}

#define read_vector read_vertex

int read_texcoord(FILE* f, TexCoord* t)
{
 int nr;
 if (!f || !t) panic(HERE);
 nr = fscanf(f, "(%lf,%lf,%lf)", &t->x, &t->y, &t->z);
 if (nr != 3) return 0;
 if (t->x < 0. || t->x > 1.) return 0;
 if (t->y < 0. || t->y > 1.) return 0;
 if (t->z < 0. || t->z > 1.) return 0;
 return 1;
}

int read_triangle(FILE* f, Triangle* t, int idx)
{
 unsigned long pos;
 int nr;
 char str[1024];

 if (!f || !t) panic(HERE);
 
 pos = ftell(f);
 
 fscanf(f, "{\n");
 t->idx = idx;
 
 fscanf(f, " a: ");
 if (!read_vertex(f, &t->a)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " b: ");
 if (!read_vertex(f, &t->b)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " c: ");
 if (!read_vertex(f, &t->c)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " d: ");
 if (!read_vertex(f, &t->d)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " texA: ");
 if (!read_texcoord(f, &t->ta)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " texB: ");
 if (!read_texcoord(f, &t->tb)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " texC: ");
 if (!read_texcoord(f, &t->tc)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " texD: ");
 if (!read_texcoord(f, &t->td)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " na: ");
 if (!read_vector(f, &t->na)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " nb: ");
 if (!read_vector(f, &t->nb)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " nc: ");
 if (!read_vector(f, &t->nc)) panic(HERE);
 fscanf(f,"\n");

 fscanf(f, " nd: ");
 if (!read_vector(f, &t->nd)) panic(HERE);
 fscanf(f,"\n");

 if (length(&t->na) <= precision || length(&t->nb) <= precision || length(&t->nc) <= precision)
   {
    compute_normals(t);
   }
 else
   {
    normalize(&t->na);
    normalize(&t->nb);
    normalize(&t->nc);
	normalize(&t->nd);
   }

 fscanf(f,"\n");

 nr = fscanf(f," ca: (%lf,%lf,%lf)\n", &t->mra.c, &t->mga.c, &t->mba.c);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," cb: (%lf,%lf,%lf)\n", &t->mrb.c, &t->mgb.c, &t->mbb.c);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," cc: (%lf,%lf,%lf)\n", &t->mrc.c, &t->mgc.c, &t->mbc.c);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," cd: (%lf,%lf,%lf)\n", &t->mrd.c, &t->mgd.c, &t->mbd.c);
 if (nr != 3) panic(HERE);

 nr = fscanf(f," sa: (%lf,%lf,%lf)\n", &t->mra.s, &t->mga.s, &t->mba.s);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," sb: (%lf,%lf,%lf)\n", &t->mrb.s, &t->mgb.s, &t->mbb.s);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," sc: (%lf,%lf,%lf)\n", &t->mrc.s, &t->mgc.s, &t->mbc.s);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," sd: (%lf,%lf,%lf)\n", &t->mrd.s, &t->mgd.s, &t->mbd.s);
 if (nr != 3) panic(HERE);

 nr = fscanf(f," ta: (%lf,%lf,%lf)\n", &t->mra.t, &t->mga.t, &t->mba.t);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," tb: (%lf,%lf,%lf)\n", &t->mrb.t, &t->mgb.t, &t->mbb.t);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," tc: (%lf,%lf,%lf)\n", &t->mrc.t, &t->mgc.t, &t->mbc.t);
 if (nr != 3) panic(HERE);
 nr = fscanf(f," td: (%lf,%lf,%lf)\n", &t->mrd.t, &t->mgd.t, &t->mbd.t);
 if (nr != 3) panic(HERE);

 pos = ftell(f);
 nr = fscanf(f," surface: (%lf,%lf,%lf,%lf,%lf)\n", &t->s.A, &t->s.B, &t->s.C, &t->s.D, &t->s.E);
 if (nr != 5)
   {
    t->s.A = t->s.B = t->s.C = t->s.D = t->s.E = 0.;
    fseek(f, pos, SEEK_SET);
   }

 normalize_color_factors(t);

 if (t->s.A == 0. && t->s.B == 0. && t->s.C == 0. && t->s.D == 0. && t->s.E == 0.) compute_surface(t);
 
 nr = fscanf(f, " %s (%lf,%lf)\n", str, &t->mUR, &t->mDR);
 if (nr != 3) panic(HERE);
 if (!strcmp(str, "tf:"))
   {
    t->mUG = t->mUR;
    t->mUB = t->mUR;
    t->mDG = t->mDR;
    t->mDB = t->mDR;
   }
 else if (!strcmp(str, "tfr:"))
   {
    nr = fscanf(f, " tfg: (%lf,%lf)\n", &t->mUG, &t->mDG);
    if (nr != 2) panic(HERE);
    nr = fscanf(f, " tfb: (%lf,%lf)\n", &t->mUB, &t->mDB);
    if (nr != 2) panic(HERE);
   }
 else panic(HERE);

 if (t->mUR < .05 || t->mUR > 20.) panic(HERE);
 if (t->mDR < .05 || t->mDR > 20.) panic(HERE);
 if (t->mUG < .05 || t->mUG > 20.) panic(HERE);
 if (t->mDG < .05 || t->mDG > 20.) panic(HERE);
 if (t->mUB < .05 || t->mUB > 20.) panic(HERE);
 if (t->mDB < .05 || t->mDB > 20.) panic(HERE);

 pos = ftell(f);
 nr = fscanf(f, " sf: (%lf,%lf,%lf)\n", &t->caR, &t->caG, &t->caB);
 if (nr != 3) panic(HERE);
 if (t->caR < -1.) panic(HERE);
 if (t->caG < -1.) panic(HERE);
 if (t->caB < -1.) panic(HERE);

 nr = fscanf(f, " texture: %d\n", &t->tid);
 if (nr != 1) panic(HERE);
 t->tex = &texture[0];
 fscanf(f, "}\n");
 return 1;
}

void compute_slimits(int n)
{
 if (n < 1000)
   {
    slimit_max = n/2;
    slimit_min = n/2;
    slimit_step = 0.0;
   }
 else if (n >= 1000 && n < 2000)
   {
    slimit_max = n/4;
    slimit_min = n/8;
    slimit_step = (4.*(slimit_max - slimit_min))/n;
   }
 else if (n >= 2000 && n < 4000)
   {
    slimit_max = n/10;
    slimit_min = n/40;
    slimit_step = (3.5*(slimit_max - slimit_min))/n;
   }
 else if (n >= 4000 && n < 8000)
   {
    slimit_max = n/25;
    slimit_min = n/125;
    slimit_step = (3.*(slimit_max - slimit_min))/n;
   }
 else if (n >= 8000 && n < 16000)
   {
    slimit_max = n/100;
    slimit_min = n/500;
    slimit_step = (2.5*(slimit_max - slimit_min))/n;
   }
 else if (n >= 16000 && n < 32000)
   {
    slimit_max = sqrt((REAL)n);
    slimit_min = n/1600;
    slimit_step = (2.0*(slimit_max - slimit_min))/n;
   }
 else if (n >= 32000 && n < 64000)
   {
    slimit_max = sqrt((REAL)n)/2.;
    slimit_min = 10;
    slimit_step = (1.75*(slimit_max - slimit_min))/n;
   }
 else
   {
    slimit_max = sqrt((REAL)n)/2.;
    slimit_min = 10;
    slimit_step = (1.5*(slimit_max - slimit_min))/n;
   }
}

// FIXME: is it OK, hypersurface involved
REAL calc_surface(Triangle* t)
{
 REAL alfa,beta,gamma;
 REAL la,lb,lc,res;
 Vector a,b,c;
 a.x = t->c.x - t->a.x;
 a.y = t->c.y - t->a.y;
 a.z = t->c.z - t->a.z;
 a.v = t->c.v - t->a.v;
 b.x = t->b.x - t->a.x;
 b.y = t->b.y - t->a.y;
 b.z = t->b.z - t->a.z;
 b.v = t->b.v - t->a.v;
 c.x = t->d.x - t->a.x;
 c.y = t->d.y - t->a.y;
 c.z = t->d.z - t->a.z;
 c.v = t->d.v - t->a.v;
 la = pseudo_length(&a);
 lb = pseudo_length(&b);
 lc = pseudo_length(&c);
 if (la < 3e-6 || lb < 3e-6 || lc < 3e-6)	
   {
    return 0.;
   }
 normalize(&a);
 normalize(&b);
 normalize(&c);
 alfa  = acos(scalar_prod(&a, &b));
 beta  = acos(scalar_prod(&a, &c));
 gamma = acos(scalar_prod(&b, &c));

 res = sin(alfa) * sin(beta) * sin(gamma) * la * lb * lc / 1000000.;
 //printf("hypersurface = %lf\n", res);
 return res;
}

REAL calc_box(Triangle* t1, Triangle* t2)
{
 REAL minx,miny,minz,minv;
 REAL maxx,maxy,maxz,maxv;
 REAL dx,dy,dz,dv,res;

 minx = t1->a.x;
 if (t1->b.x < minx) minx = t1->b.x;
 if (t1->c.x < minx) minx = t1->c.x;
 if (t1->d.x < minx) minx = t1->d.x;
 if (t2->a.x < minx) minx = t2->a.x;
 if (t2->b.x < minx) minx = t2->b.x;
 if (t2->c.x < minx) minx = t2->c.x;
 if (t2->d.x < minx) minx = t2->d.x;
 miny = t1->a.y;
 if (t1->b.y < miny) miny = t1->b.y;
 if (t1->c.y < miny) miny = t1->c.y;
 if (t1->d.y < miny) miny = t1->d.y;
 if (t2->a.y < miny) miny = t2->a.y;
 if (t2->b.y < miny) miny = t2->b.y;
 if (t2->c.y < miny) miny = t2->c.y;
 if (t2->d.y < miny) miny = t2->d.y;
 minz = t1->a.z;
 if (t1->b.z < minz) minz = t1->b.z;
 if (t1->c.z < minz) minz = t1->c.z;
 if (t1->d.z < minz) minz = t1->d.z;
 if (t2->a.z < minz) minz = t2->a.z;
 if (t2->b.z < minz) minz = t2->b.z;
 if (t2->c.z < minz) minz = t2->c.z;
 if (t2->d.z < minz) minz = t2->d.z;
 minv = t1->a.v;
 if (t1->b.v < minv) minv = t1->b.v;
 if (t1->c.v < minv) minv = t1->c.v;
 if (t1->d.v < minv) minv = t1->d.v;
 if (t2->a.v < minv) minv = t2->a.v;
 if (t2->b.v < minv) minv = t2->b.v;
 if (t2->c.v < minv) minv = t2->c.v;
 if (t2->d.v < minv) minv = t2->d.v;
 maxx = t1->a.x;
 if (t1->b.x > maxx) maxx = t1->b.x;
 if (t1->c.x > maxx) maxx = t1->c.x;
 if (t1->d.x > maxx) maxx = t1->d.x;
 if (t2->a.x > maxx) maxx = t2->a.x;
 if (t2->b.x > maxx) maxx = t2->b.x;
 if (t2->c.x > maxx) maxx = t2->c.x;
 if (t2->d.x > maxx) maxx = t2->d.x;
 maxy = t1->a.y;
 if (t1->b.y > maxy) maxy = t1->b.y;
 if (t1->c.y > maxy) maxy = t1->c.y;
 if (t1->d.y > maxy) maxy = t1->d.y;
 if (t2->a.y > maxy) maxy = t2->a.y;
 if (t2->b.y > maxy) maxy = t2->b.y;
 if (t2->c.y > maxy) maxy = t2->c.y;
 if (t2->d.y > maxy) maxy = t2->d.y;
 maxz = t1->a.z;
 if (t1->b.z > maxz) maxz = t1->b.z;
 if (t1->c.z > maxz) maxz = t1->c.z;
 if (t1->d.z > maxz) maxz = t1->d.z;
 if (t2->a.z > maxz) maxz = t2->a.z;
 if (t2->b.z > maxz) maxz = t2->b.z;
 if (t2->c.z > maxz) maxz = t2->c.z;
 if (t2->d.z > maxz) maxz = t2->d.z;
 maxv = t1->a.v;
 if (t1->b.v > maxv) maxv = t1->b.v;
 if (t1->c.v > maxv) maxv = t1->c.v;
 if (t1->d.v > maxv) maxv = t1->d.v;
 if (t2->a.v > maxv) maxv = t2->a.v;
 if (t2->b.v > maxv) maxv = t2->b.v;
 if (t2->c.v > maxv) maxv = t2->c.v;
 if (t2->d.v > maxv) maxv = t2->d.v;
 dx = maxx - minx;
 dy = maxy - miny;
 dz = maxz - minz;
 dv = maxv - minv;

 res = dx*dx + dy*dy + dz*dz + dv*dv;
 //printf("calcbox = %lf\n", res);
 return res;
}


void blist_add(BList** head, REAL ix, REAL iy, REAL iz, REAL iv, REAL ax, REAL ay, REAL az, REAL av, Triangle* t1, Triangle* t2)
{
 BList* temp;
 if (*head == NULL)
   {
    *head = (BList*)(malloc(sizeof(BList)));

    (*head)->b.t1 = t1;
    (*head)->b.t2 = t2;
    (*head)->b.minx = ix;
    (*head)->b.miny = iy;
    (*head)->b.minz = iz;
	(*head)->b.minv = iv;
    (*head)->b.maxx = ax;
    (*head)->b.maxy = ay;
    (*head)->b.maxz = az;
	(*head)->b.maxv = av;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }

  temp = (BList*)(malloc(sizeof(BList)));

  temp->b.t1 = t1;
  temp->b.t2 = t2;
  temp->b.minx = ix;
  temp->b.miny = iy;
  temp->b.minz = iz;
  temp->b.minv = iv;
  temp->b.maxx = ax;
  temp->b.maxy = ay;
  temp->b.maxz = az;
  temp->b.maxv = av;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
}

REAL calc_boxes(BTree* b1, BTree* b2)
{
 REAL minx,miny,minz,minv;
 REAL maxx,maxy,maxz,maxv;
 REAL dx,dy,dz,dv,res;
 minx = b1->b->minx;
 if (b2->b->minx < minx) minx = b2->b->minx;
 miny = b1->b->miny;
 if (b2->b->miny < miny) miny = b2->b->miny;
 minz = b1->b->minz;
 if (b2->b->minz < minz) minz = b2->b->minz;
 minv = b1->b->minv;
 if (b2->b->minv < minv) minv = b2->b->minv;
 maxx = b1->b->maxx;
 if (b2->b->maxx > maxx) maxx = b2->b->maxx;
 maxy = b1->b->maxy;
 if (b2->b->maxy > maxy) maxy = b2->b->maxy;
 maxz = b1->b->maxz;
 if (b2->b->maxz > maxz) maxz = b2->b->maxz;
 maxv = b1->b->maxv;
 if (b2->b->maxv > maxv) maxv = b2->b->maxv;
 dx = maxx - minx;
 dy = maxy - miny;
 dz = maxz - minz;
 dv = maxv - minv;
 res = dx*dx+dy*dy+dz*dz+dv*dv;
 return res;
}

void create_btree_box(Box** b, BTree* b1, BTree* b2)
{
 Box* temp;
 REAL minx,miny,minz,minv;
 REAL maxx,maxy,maxz,maxv;
 minx = b1->b->minx;
 if (b2->b->minx < minx) minx = b2->b->minx;
 miny = b1->b->miny;
 if (b2->b->miny < miny) miny = b2->b->miny;
 minz = b1->b->minz;
 if (b2->b->minz < minz) minz = b2->b->minz;
 minv = b1->b->minv;
 if (b2->b->minv < minv) minv = b2->b->minv;
 maxx = b1->b->maxx;
 if (b2->b->maxx > maxx) maxx = b2->b->maxx;
 maxy = b1->b->maxy;
 if (b2->b->maxy > maxy) maxy = b2->b->maxy;
 maxz = b1->b->maxz;
 if (b2->b->maxz > maxz) maxz = b2->b->maxz;
 maxv = b1->b->maxv;
 if (b2->b->maxv > maxv) maxv = b2->b->maxv;
 temp = (Box*)malloc(sizeof(Box));
 temp->t1 = temp->t2 = NULL;	/* abstract box-> node */
 temp->minx = minx;
 temp->miny = miny;
 temp->minz = minz;
 temp->minv = minv;
 temp->maxx = maxx;
 temp->maxy = maxy;
 temp->maxz = maxz;
 temp->maxv = maxv;
 *b = temp;
}


void add_box(BList** head, Triangle* t1, Triangle* t2)
{
 REAL minx,miny,minz,minv;
 REAL maxx,maxy,maxz,maxv;
/* if (!head || ! t1) spanic(HERE, ":", HERE);*/
 minx = t1->a.x;
 if (t1->b.x < minx) minx = t1->b.x;
 if (t1->c.x < minx) minx = t1->c.x;
 if (t1->d.x < minx) minx = t1->d.x;
 if (t2)
   {
    if (t2->a.x < minx) minx = t2->a.x;
    if (t2->b.x < minx) minx = t2->b.x;
    if (t2->c.x < minx) minx = t2->c.x;
	if (t2->d.x < minx) minx = t2->d.x;
   }
 miny = t1->a.y;
 if (t1->b.y < miny) miny = t1->b.y;
 if (t1->c.y < miny) miny = t1->c.y;
 if (t1->d.y < miny) miny = t1->d.y;
 if (t2)
   {
    if (t2->a.y < miny) miny = t2->a.y;
    if (t2->b.y < miny) miny = t2->b.y;
    if (t2->c.y < miny) miny = t2->c.y;
	if (t2->d.y < miny) miny = t2->d.y;
   }
 minz = t1->a.z;
 if (t1->b.z < minz) minz = t1->b.z;
 if (t1->c.z < minz) minz = t1->c.z;
 if (t1->d.z < minz) minz = t1->d.z;
 if (t2)
   {
    if (t2->a.z < minz) minz = t2->a.z;
    if (t2->b.z < minz) minz = t2->b.z;
    if (t2->c.z < minz) minz = t2->c.z;
	if (t2->d.z < minz) minz = t2->d.z;
   }
 minv = t1->a.v;
 if (t1->b.v < minv) minv = t1->b.v;
 if (t1->c.v < minv) minv = t1->c.v;
 if (t1->d.v < minv) minv = t1->d.v;
 if (t2)
   {
    if (t2->a.v < minv) minv = t2->a.v;
    if (t2->b.v < minv) minv = t2->b.v;
    if (t2->c.v < minv) minv = t2->c.v;
	if (t2->d.v < minv) minv = t2->d.v;
   }
 maxx = t1->a.x;
 if (t1->b.x > maxx) maxx = t1->b.x;
 if (t1->c.x > maxx) maxx = t1->c.x;
 if (t1->d.x > maxx) maxx = t1->d.x;
 if (t2)
   {
    if (t2->a.x > maxx) maxx = t2->a.x;
    if (t2->b.x > maxx) maxx = t2->b.x;
    if (t2->c.x > maxx) maxx = t2->c.x;
	if (t2->d.x > maxx) maxx = t2->d.x;
   }
 maxy = t1->a.y;
 if (t1->b.y > maxy) maxy = t1->b.y;
 if (t1->c.y > maxy) maxy = t1->c.y;
 if (t1->d.y > maxy) maxy = t1->d.y;
 if (t2)
   {
    if (t2->a.y > maxy) maxy = t2->a.y;
    if (t2->b.y > maxy) maxy = t2->b.y;
    if (t2->c.y > maxy) maxy = t2->c.y;
	if (t2->d.y > maxy) maxy = t2->d.y;
   }
 maxz = t1->a.z;
 if (t1->b.z > maxz) maxz = t1->b.z;
 if (t1->c.z > maxz) maxz = t1->c.z;
 if (t1->d.z > maxz) maxz = t1->d.z;
 if (t2)
   {
    if (t2->a.z > maxz) maxz = t2->a.z;
    if (t2->b.z > maxz) maxz = t2->b.z;
    if (t2->c.z > maxz) maxz = t2->c.z;
	if (t2->d.z > maxz) maxz = t2->d.z;
   }
 maxv = t1->a.v;
 if (t1->b.v > maxv) maxv = t1->b.v;
 if (t1->c.v > maxv) maxv = t1->c.v;
 if (t1->d.v > maxv) maxv = t1->d.v;
 if (t2)
   {
    if (t2->a.v > maxv) maxv = t2->a.v;
    if (t2->b.v > maxv) maxv = t2->b.v;
    if (t2->c.v > maxv) maxv = t2->c.v;
	if (t2->d.v > maxv) maxv = t2->d.v;
   }
 blist_add(head, minx, miny, minz, minv, maxx, maxy, maxz, maxv, t1, t2);
}

BTree* find_nearest_boxes_smart(PtrList* head, PtrList* tail, PtrList** b1, PtrList** b2, int swap, int slimit)
{
 BTree* nbox;
 PtrList *p1, *p2;
 PtrList *bl, *br;
 REAL minD,d;
 int i,j;
 p1 = head;
 minD = 1e15;
 *b1 = *b2 = NULL;
 bl = br = NULL;
  i = j = 0;
  while (p1)
   {
    p2 = p1->next;
    while (p2)
      {
       if ( (d = calc_boxes((BTree*)p1->ptr, (BTree*)p2->ptr)) < minD)
         {
	  minD = d;
          bl = p1;
	  br = p2;
         }
       j++;
       p2 = p2->next;
       if (swap == 0 && j >= slimit) p2 = NULL;	/* skipper for speed */
      }
    j = 0;
    i++;
    p1 = p1->next;
    if (swap == 1 && i >= slimit) p1 = NULL;
   }
 if (!bl || !br) panic(HERE);
 nbox = (BTree*)malloc(sizeof(BTree));
 if (swap)
   {
    nbox->l = (BTree*)bl->ptr;
    nbox->r = (BTree*)br->ptr;
    create_btree_box(&nbox->b, nbox->l, nbox->r);
   }
 else
   {
    nbox->r = (BTree*)bl->ptr;
    nbox->l = (BTree*)br->ptr;
    create_btree_box(&nbox->b, nbox->r, nbox->l);
   }
 *b1 = bl;
 *b2 = br;
 return nbox;
}



void create_btree(BList* head, BTree** a_btree)
{
 BList *hd;
 PtrList* phead;
 PtrList* ptail;
 BTree* temp;
 PtrList *bt1, *bt2;
 int no;
 REAL slimit;
 hd = head;
 phead = NULL;
 bt1 = bt2 = NULL;

 no = 0;
 while (head)
   {
    temp = (BTree*)malloc(sizeof(BTree));
    temp->l = temp->r = NULL;	/* real box: leaf */
    temp->b = &head->b;		/* uses a real box containing triangle(s) */
    ptrlist_add(&phead, (void*)temp);
    head = head->next;
   }
 no = count_items(phead);
 ptail = ptrlist_get_tail(phead);
 if (no != count_items(phead)) panic(HERE);
 slimit = slimit_min;

 while (no > 1)
   {
    temp = find_nearest_boxes_smart(phead, ptail, &bt1, &bt2, (no % 2), (int)slimit);
    ptrlist_delete(&phead, bt1);
    ptrlist_delete(&phead, bt2);
    ptrlist_add(&phead, (void*)temp);
    slimit += slimit_step;
    if (slimit > slimit_max) slimit = slimit_max;
    if (slimit > no) slimit = no;
    no --;
   }
 if (!phead) phead = ptail;
 if (!a_btree) btree = (BTree*)phead->ptr;
 else *a_btree = (BTree*)phead->ptr;
 ptrlist_free(&phead);
}

void process_voxel(Voxel* vx, int whendiv, int ndiv, short* sarray)
{
 PtrList* temp;
 Triangle *t, *t1, *t2;
 int n, idx;
 int ix,iy,iz,iv,i,ndiv2,ndiv3,ii;
 REAL minS,s;
 REAL x,y,z,v;
 REAL dx,dy,dz,dv;
 PtrList *head, *hd, *minT;
 BList * bhead;
 BTree* b_head;
 temp = vx->tris;
 
 vx->varr = NULL;

 if (vx->ntbox > whendiv)
   {
    n = ndiv * ndiv * ndiv * ndiv;
    vx->varr = (Voxel*)malloc(n*sizeof(Voxel));
    ndiv2 = ndiv * ndiv;
	ndiv3 = ndiv2 * ndiv;
    for (i=0;i<n;i++)
      {
       vx->varr[i].tris = NULL;
       dx = (vx->box.maxx - vx->box.minx);
       dy = (vx->box.maxy - vx->box.miny);
       dz = (vx->box.maxz - vx->box.minz);
	   dv = (vx->box.maxv - vx->box.minv);

	   //FIXME: unchecked!
       iv = i % ndiv;
       iz = (i / ndiv) % ndiv;
	   iy = (i / ndiv2) % ndiv;
       ix = (i / (ndiv3));
       
       vx->varr[i].box.minx = vx->box.minx + ((REAL)ix / (REAL)ndiv) * dx;
       vx->varr[i].box.miny = vx->box.miny + ((REAL)iy / (REAL)ndiv) * dy;
       vx->varr[i].box.minz = vx->box.minz + ((REAL)iz / (REAL)ndiv) * dz;
	   vx->varr[i].box.minv = vx->box.minv + ((REAL)iv / (REAL)ndiv) * dv;
       
       vx->varr[i].box.maxx = vx->box.minx + ((REAL)(ix+1) / (REAL)ndiv) * dx;
       vx->varr[i].box.maxy = vx->box.miny + ((REAL)(iy+1) / (REAL)ndiv) * dy;
       vx->varr[i].box.maxz = vx->box.minz + ((REAL)(iz+1) / (REAL)ndiv) * dz;
	   vx->varr[i].box.maxv = vx->box.minv + ((REAL)(iv+1) / (REAL)ndiv) * dv;

       vx->varr[i].ntbox = 0;
       vx->varr[i].varr = NULL;
      }
    while (temp)
      {
       t = (Triangle*)temp->ptr;
       
	   //FIXME: unchecked, test this carefully!
       x = (t->a.x + t->b.x + t->c.x + t->d.x) / 4.;
       y = (t->a.y + t->b.y + t->c.y + t->d.y) / 4.;
       z = (t->a.z + t->b.z + t->c.z + t->d.z) / 4.;
	   v = (t->a.v + t->b.v + t->c.v + t->d.v) / 4.;
       
       ix = (int)(((x - vx->box.minx) / (vx->box.maxx - vx->box.minx)) * ndiv);
       iy = (int)(((y - vx->box.miny) / (vx->box.maxy - vx->box.miny)) * ndiv);
       iz = (int)(((z - vx->box.minz) / (vx->box.maxz - vx->box.minz)) * ndiv);
	   iv = (int)(((v - vx->box.minv) / (vx->box.maxv - vx->box.minv)) * ndiv);

       if (ix < 0) ix = 0;
       if (iy < 0) iy = 0;
       if (iz < 0) iz = 0;
	   if (iv < 0) iv = 0;
       if (ix >= ndiv) ix = ndiv - 1;
       if (iy >= ndiv) iy = ndiv - 1;
       if (iz >= ndiv) iz = ndiv - 1;
	   if (iv >= ndiv) iv = ndiv - 1;

       idx = ndiv3 * ix + ndiv2 * iy + ndiv * iz + iv;

       ptrlist_add(&(vx->varr[idx].tris), t);
       vx->varr[idx].ntbox ++;
       
       temp = temp->next;
      }
    for (i=0;i<n;i++)
      {
       if (vx->varr[i].ntbox > 0) process_voxel(&vx->varr[i], whendiv, ndiv, sarray);
      }
   }
 else  
   {
	   //FIXME: wych 3 * vx->ntbox ?
       compute_slimits(3*vx->ntbox);
       head = vx->tris;
       bhead = NULL;
       ii = 0;
       
       while (head)
       {
        ii++;

        t1 = t2 = NULL;
        minS = 1e12;
        s = minS;
        hd = head;
        minT = head;
        while (head)
          {
           if ((s = calc_surface((Triangle*)head->ptr)) < minS)
             {
              minS = s;
	      minT = head;
             }
          head = head->next;
         } 
      t1 = (Triangle*)minT->ptr;
      head = hd;
      ptrlist_delete(&head, minT);
      hd = head;
      minT = head;
      minS = 1e12;
      s = minS;
      if (head)
        {
         while (head)
           {
            if ((s = calc_box((Triangle*)head->ptr, t1)) < minS)
              {
               minS = s;
               minT = head;
              }
            head = head->next;
           }
       t2 = (Triangle*)minT->ptr;
       head = hd;
       ptrlist_delete(&head, minT);
       hd = head;
      }
    add_box(&bhead, t1, t2);
   }
   ptrlist_free(&head);

   create_btree(bhead, &b_head);
   ptrlist_add(&pbv_head, (void*)b_head);
  }
}

void merge_voxel_trees(PtrList* phead, int no)
{
 PtrList* ptail;
 BTree* temp;
 PtrList *bt1, *bt2;
 REAL slimit;

 ptail = ptrlist_get_tail(phead);
 if (no != count_items(phead)) panic(HERE);
 slimit = slimit_min;

 while (no > 1)
   {    
    temp = find_nearest_boxes_smart(phead, ptail, &bt1, &bt2, (no % 2), (int)slimit);
    ptrlist_delete(&phead, bt1);
    ptrlist_delete(&phead, bt2);
    ptrlist_add(&phead, (void*)temp);
    slimit += slimit_step;
    if (slimit > slimit_max) slimit = slimit_max;
    if (slimit > no) slimit = no;
    no --;
   }
 if (!phead) phead = ptail;
 if (phead)btree = (BTree*)phead->ptr;
 else panic(HERE);
 ptrlist_free(&phead);
}

void preprocess_scene()
{
 short *tv_idx;
 int i,vx_ndiv,vx_whendiv;
 REAL mi, mI;
 Voxel v_root;
 vx_ndiv = 2;
 vx_whendiv = 5000;

 tv_idx = (short*)malloc(nTriangles*sizeof(short));
 v_root.box.minx = 1e10;
 v_root.box.miny = 1e10;
 v_root.box.minz = 1e10;
 v_root.box.minv = 1e10;
 v_root.box.maxx = -1e10;
 v_root.box.maxy = -1e10;
 v_root.box.maxz = -1e10;
 v_root.box.maxv = -1e10;
 v_root.tris = NULL;
 v_root.ntbox = 0;
 v_root.varr = NULL;
 
 for (i=0;i<nTriangles;i++) tv_idx[i] = 0;
 for (i=0;i<nTriangles;i++) 
   {
    ptrlist_add(&v_root.tris, &g_ts[i]);

    mI = mi = g_ts[i].a.x;
    if (g_ts[i].b.x > mI) mI = g_ts[i].b.x;
    if (g_ts[i].c.x > mI) mI = g_ts[i].c.x;
	if (g_ts[i].d.x > mI) mI = g_ts[i].d.x;
    if (mI > v_root.box.maxx) v_root.box.maxx = mI;
    if (g_ts[i].b.x < mi) mi = g_ts[i].b.x;
    if (g_ts[i].c.x < mi) mi = g_ts[i].c.x;
	if (g_ts[i].d.x < mi) mi = g_ts[i].d.x;
    if (mi < v_root.box.minx) v_root.box.minx = mi;

    mI = mi = g_ts[i].a.y;
    if (g_ts[i].b.y > mI) mI = g_ts[i].b.y;
    if (g_ts[i].c.y > mI) mI = g_ts[i].c.y;
	if (g_ts[i].d.y > mI) mI = g_ts[i].d.y;
    if (mI > v_root.box.maxy) v_root.box.maxy = mI;
    if (g_ts[i].b.y < mi) mi = g_ts[i].b.y;
    if (g_ts[i].c.y < mi) mi = g_ts[i].c.y;
	if (g_ts[i].d.y < mi) mi = g_ts[i].d.y;
    if (mi < v_root.box.miny) v_root.box.miny = mi;

    mI = mi = g_ts[i].a.z;
    if (g_ts[i].b.z > mI) mI = g_ts[i].b.z;
    if (g_ts[i].c.z > mI) mI = g_ts[i].c.z;
	if (g_ts[i].d.z > mI) mI = g_ts[i].d.z;
    if (mI > v_root.box.maxz) v_root.box.maxz = mI;
    if (g_ts[i].b.z < mi) mi = g_ts[i].b.z;
    if (g_ts[i].c.z < mi) mi = g_ts[i].c.z;
	if (g_ts[i].d.z < mi) mi = g_ts[i].d.z;
    if (mi < v_root.box.minz) v_root.box.minz = mi;

    mI = mi = g_ts[i].a.v;
    if (g_ts[i].b.v > mI) mI = g_ts[i].b.v;
    if (g_ts[i].c.v > mI) mI = g_ts[i].c.v;
	if (g_ts[i].d.v > mI) mI = g_ts[i].d.v;
    if (mI > v_root.box.maxv) v_root.box.maxv = mI;
    if (g_ts[i].b.v < mi) mi = g_ts[i].b.v;
    if (g_ts[i].c.v < mi) mi = g_ts[i].c.v;
	if (g_ts[i].d.v < mi) mi = g_ts[i].d.v;
    if (mi < v_root.box.minv) v_root.box.minv = mi;
   }
 v_root.ntbox = nTriangles;

 pbv_head = NULL;
 process_voxel(&v_root, vx_whendiv, vx_ndiv, tv_idx);
 i = count_items(pbv_head);

 compute_slimits(i);
 merge_voxel_trees(pbv_head, i);
}

int get_red(Screen* s, int x, int y, int z)
{
 if (x == s->x) x--;
 if (y == s->y) y--;
 if (z == s->z) z--;
 return s->pixels[3*(s->z * s->y * x + s->y * y + z)];
}


int get_green(Screen* s, int x, int y, int z)
{
 if (x == s->x) x--;
 if (y == s->y) y--;
 if (z == s->z) z--;
 return s->pixels[3*(s->z * s->y * x + s->y * y + z)+1];
}


int get_blue(Screen* s, int x, int y, int z)
{
 if (x == s->x) x--;
 if (y == s->y) y--;
 if (z == s->z) z--;
 return s->pixels[3*(s->z * s->y * x + s->y * y + z)+2];
}

void get_texturei(Texture* t, REAL x, REAL y, REAL z, int *r, int* g, int* b)
{
 int i,j,k;
 i = (int)(x*(REAL)t->x);
 j = (int)(y*(REAL)t->y);
 k = (int)(z*(REAL)t->z);

 *r = get_red(t, i, j, k);
 *g = get_green(t, i, j, k);
 *b = get_blue(t, i, j, k);
}

void set_red(Texture* s, int x, int y, int z, unsigned char c)
{
 s->pixels[3*(s->z * s->y * x + s->y * y + z)] = c;
}

void set_green(Texture* s, int x, int y, int z, unsigned char c)
{
 s->pixels[3*(s->z * s->y * x + s->y * y + z)+1] = c;
}

void set_blue(Texture* s, int x, int y, int z, unsigned char c)
{
 s->pixels[3*(s->z * s->y * x + s->y * y + z)+2] = c;
}

void set_texture(Texture* t, int i, int j, int k, unsigned char r, unsigned char g, unsigned char b)
{
	set_red(t, i, j, k, r);
	set_green(t, i, j, k, g);
	set_blue(t, i, j, k, b);
}

void make_stub_texture(Texture* t)
{
 int n,i,j,k;
 t->x = 0x80;
 t->y = 0x80;
 t->z = 0x80;
 n = t->x * t->y * t->z * 3;
 t->pixels = (unsigned char*)malloc( n );
 if (!t->pixels) panic(HERE);

 for (i=0;i<0x80;i++) for (j=0;j<0x80;j++) for (k=0;k<0x80;k++)
  {
	 set_texture(t, i, j, k, i<<1, j<<1, k<<1);
  }
}

void make_textures()
{
	nTex = 1;
	texture = (Texture*)malloc(nTex*sizeof(Texture));
	make_stub_texture(&texture[0]);
}


void load_scene(char* scenefile)
{
 int sx, sy, sz;
 int i, n, nr;
 unsigned long pos;
 FILE* f;
 char str[1024];

 max_rec = 16;

 f = fopen(scenefile,"rb");
 if (!f) panic(HERE);

 nr = fscanf(f, "Screen: (%d,%d,%d)\n", &sx, &sy, &sz);
 if (nr != 3) panic(HERE);

 if (sx <= 0 || sy <= 0 || sz <= 0) panic(HERE);
 
 pos = ftell(f);
 nr = fscanf(f, "%s", str);
 if (nr != 1) panic(HERE);
 if (!strcmp(str, "Background:"))
   {
    nr = fscanf(f, " %d\n", &bkgnd);
    if (nr != 1)   panic(HERE);
    if (bkgnd < 0) panic(HERE);
   }
 else 
   {
    fseek(f, pos, SEEK_SET);
    bkgnd = 0;
   }

 pos = ftell(f);
 nr = fscanf(f, "%s", str);
 if (nr != 1) panic(HERE);
 if (!strcmp(str, "MaxRecurse:"))
   {
    nr = fscanf(f, " %d\n", &max_rec);
    if (nr != 1) panic(HERE);
    if (max_rec < 0) max_rec = 0;
   }
 else fseek(f, pos, SEEK_SET);

 make_textures();
 init_screen(&screen, sx, sy, sz);

 
 fscanf(f, "Observer: ");
 if (!read_vertex(f, &observer)) panic(HERE);
 fscanf(f,"\n");

 pos = ftell(f);
 nr = fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "LookZ:"))
   {
    nr = fscanf(f, "%lf%%\n", &lookz);
    if (nr != 1) panic(HERE);
    lookz /= 100.;
   }
 else fseek(f, pos, SEEK_SET);

 pos = ftell(f);
 nr = fscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "Light:"))
   {
    if (!read_vertex(f, &light)) panic(HERE);
   }
 else panic(HERE);
 fscanf(f,"\n");

 nr = fscanf(f, "nTriangles: %d\n", &n);
 if (nr != 1) panic(HERE);
 if (n < 0) panic(HERE);
 nTriangles = n;
 if (n > 0) g_ts = (Triangle*)malloc(n*sizeof(Triangle));

 for (i=0;i<n;i++) if (!read_triangle(f, &(g_ts[i]), i)) panic(HERE);

 fclose(f);
 fprintf(stdout, "\nFinally %d triangles.\n", nTriangles);
}

void blist_add_box(BList** head, Box* b)
{
 BList* temp;
 if (*head == NULL)
   {

    *head = (BList*)(malloc(sizeof(BList)));

    (*head)->b.t1 = NULL;	/* this is STUB for -f mode only */
    (*head)->b.t2 = NULL;	/* to be used with '-=' travel mode */
    (*head)->b.minx = b->minx;
    (*head)->b.miny = b->miny;
    (*head)->b.minz = b->minz;
    (*head)->b.minv = b->minv;
    (*head)->b.maxx = b->maxx;
    (*head)->b.maxy = b->maxy;
    (*head)->b.maxz = b->maxz;
    (*head)->b.maxv = b->maxv;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }
  temp = (BList*)(malloc(sizeof(BList)));

  temp->b.t1 = NULL;		/* this is a STUB */
  temp->b.t2 = NULL;		/* read above */
  temp->b.minx = b->minx;
  temp->b.miny = b->miny;
  temp->b.minz = b->minz;
  temp->b.minv = b->minv;
  temp->b.maxx = b->maxx;
  temp->b.maxy = b->maxy;
  temp->b.maxz = b->maxz;
  temp->b.maxv = b->maxv;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
}

int check_btree_fit(BTree* bt)
{
 int err;
 err = 0;
 if (bt->b->t1)
 {
     if (bt->b->t1->a.x + precision < bt->b->minx || bt->b->t1->a.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t1.a.x not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->b.x + precision < bt->b->minx || bt->b->t1->b.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t1.b.x not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->c.x + precision < bt->b->minx || bt->b->t1->c.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t1.c.x not fits within box min/max\n"); err = 1; }
	 if (bt->b->t1->d.x + precision < bt->b->minx || bt->b->t1->d.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t1.d.x not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->a.y + precision < bt->b->miny || bt->b->t1->a.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t1.a.y not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->b.y + precision < bt->b->miny || bt->b->t1->b.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t1.b.y not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->c.y + precision < bt->b->miny || bt->b->t1->c.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t1.c.y not fits within box min/max\n"); err = 1; }
	 if (bt->b->t1->d.y + precision < bt->b->miny || bt->b->t1->d.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t1.d.y not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->a.z + precision < bt->b->minz || bt->b->t1->a.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t1.a.z not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->b.z + precision < bt->b->minz || bt->b->t1->b.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t1.b.z not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->c.z + precision < bt->b->minz || bt->b->t1->c.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t1.c.z not fits within box min/max\n"); err = 1; }
	 if (bt->b->t1->d.z + precision < bt->b->minz || bt->b->t1->d.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t1.d.z not fits within box min/max\n"); err = 1; }
	 if (bt->b->t1->a.v + precision < bt->b->minv || bt->b->t1->a.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t1.a.v not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->b.v + precision < bt->b->minv || bt->b->t1->b.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t1.b.v not fits within box min/max\n"); err = 1; }
     if (bt->b->t1->c.v + precision < bt->b->minv || bt->b->t1->c.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t1.c.v not fits within box min/max\n"); err = 1; }
	 if (bt->b->t1->d.v + precision < bt->b->minv || bt->b->t1->d.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t1.d.v not fits within box min/max\n"); err = 1; }
     if (err)
     {
       fprintf(stdout, "a=(%lf,%lf,%lf), b=(%lf,%lf,%lf), c=(%lf,%lf,%lf), d=(%lf,%lf,%lf,%lf)\n"
	     , bt->b->t1->a.x, bt->b->t1->a.y, bt->b->t1->a.z, bt->b->t1->a.v
		 , bt->b->t1->b.x, bt->b->t1->b.y, bt->b->t1->b.z, bt->b->t1->b.v
		 , bt->b->t1->c.x, bt->b->t1->c.y, bt->b->t1->c.z, bt->b->t1->c.v
		 , bt->b->t1->d.x, bt->b->t1->d.y, bt->b->t1->d.z, bt->b->t1->d.v);
       fprintf(stdout, "X[%lf %lf] Y[%lf %lf] Z[%lf %lf] V[%lf %lf]\n"
		   , bt->b->minx, bt->b->maxx
		   , bt->b->miny, bt->b->maxy
		   , bt->b->minz, bt->b->maxz
		   , bt->b->minv, bt->b->maxv);
       return 0;
     }
 }
 if (bt->b->t2)
 {
     if (bt->b->t2->a.x + precision < bt->b->minx || bt->b->t2->a.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t2.a.x not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->b.x + precision < bt->b->minx || bt->b->t2->b.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t2.b.x not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->c.x + precision < bt->b->minx || bt->b->t2->c.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t2.c.x not fits within box min/max\n"); err = 1; }
	 if (bt->b->t2->d.x + precision < bt->b->minx || bt->b->t2->d.x - precision > bt->b->maxx) { fprintf(stdout, "check_btree_fit: t2.d.x not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->a.y + precision < bt->b->miny || bt->b->t2->a.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t2.a.y not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->b.y + precision < bt->b->miny || bt->b->t2->b.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t2.b.y not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->c.y + precision < bt->b->miny || bt->b->t2->c.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t2.c.y not fits within box min/max\n"); err = 1; }
	 if (bt->b->t2->d.y + precision < bt->b->miny || bt->b->t2->d.y - precision > bt->b->maxy) { fprintf(stdout, "check_btree_fit: t2.d.y not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->a.z + precision < bt->b->minz || bt->b->t2->a.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t2.a.z not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->b.z + precision < bt->b->minz || bt->b->t2->b.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t2.b.z not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->c.z + precision < bt->b->minz || bt->b->t2->c.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t2.c.z not fits within box min/max\n"); err = 1; }
	 if (bt->b->t2->d.z + precision < bt->b->minz || bt->b->t2->d.z - precision > bt->b->maxz) { fprintf(stdout, "check_btree_fit: t2.d.z not fits within box min/max\n"); err = 1; }
	 if (bt->b->t2->a.v + precision < bt->b->minv || bt->b->t2->a.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t2.a.v not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->b.v + precision < bt->b->minv || bt->b->t2->b.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t2.b.v not fits within box min/max\n"); err = 1; }
     if (bt->b->t2->c.v + precision < bt->b->minv || bt->b->t2->c.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t2.c.v not fits within box min/max\n"); err = 1; }
	 if (bt->b->t2->d.v + precision < bt->b->minv || bt->b->t2->d.v - precision > bt->b->maxv) { fprintf(stdout, "check_btree_fit: t2.d.v not fits within box min/max\n"); err = 1; }
     if (err)
     {
       fprintf(stdout, "a=(%lf,%lf,%lf), b=(%lf,%lf,%lf), c=(%lf,%lf,%lf), d=(%lf,%lf,%lf,%lf)\n"
	     , bt->b->t2->a.x, bt->b->t2->a.y, bt->b->t2->a.z, bt->b->t2->a.v
		 , bt->b->t2->b.x, bt->b->t2->b.y, bt->b->t2->b.z, bt->b->t2->b.v
		 , bt->b->t2->c.x, bt->b->t2->c.y, bt->b->t2->c.z, bt->b->t2->c.v
		 , bt->b->t2->d.x, bt->b->t2->d.y, bt->b->t2->d.z, bt->b->t2->d.v);
       fprintf(stdout, "X[%lf %lf] Y[%lf %lf] Z[%lf %lf] V[%lf %lf]\n"
		   , bt->b->minx, bt->b->maxx
		   , bt->b->miny, bt->b->maxy
		   , bt->b->minz, bt->b->maxz
		   , bt->b->minv, bt->b->maxv);
       return 0;
     }
 }
 return 1;
}


int load_preprocessed(char* scenef)
{
 char fn[2048];
 FILE* f;
 int i, ii, nr, n, root, r, l;

 strcpy(fn, scenef);
 if (!strstr(fn, ".")) strcat(fn, ".btree");
 else
   {
    i = 0;
    while (fn[i] != '.') i++;
    fn[i] = 0;
    strcat(fn, ".btree");
   }

 f = fopen(fn, "r");
 if (!f) return 0;
 
 nr = fscanf(f, "nElem: %d\n", &n);
 if (nr != 1) panic(HERE);

 nr = fscanf(f, "R: %d\n", &root);
 if (nr != 1) panic(HERE);
 
 btree = (BTree*)malloc(n*sizeof(BTree));
 boxes = NULL;

 for (i=0;i<n;i++)
   {
	nr = fscanf(f, "BT: %d\n{\n", &ii);
    if (nr != 1 || i != ii) panic(HERE);

    nr = fscanf(f, " R: %d\n L: %d\n", &r, &l);
    if (nr != 2) panic(HERE);
    if (r >= 0) btree[i].r = &btree[r];
    else btree[i].r = NULL;

    if (l >= 0) btree[i].l = &btree[l];
    else btree[i].l = NULL;

    btree[i].b = (Box*)malloc(sizeof(Box));

    nr = fscanf(f," B:\n {\n  T1: %d\n  T2: %d\n", &r, &l);
    if (nr != 2) panic(HERE);

    if (r >= 0 && r < nTriangles) btree[i].b->t1 = &g_ts[r];
    else btree[i].b->t1 = NULL;

    if (l >= 0 && l < nTriangles) btree[i].b->t2 = &g_ts[l];
    else btree[i].b->t2 = NULL;

    nr = fscanf(f,"  (%lf,%lf,%lf,%lf)\n", &btree[i].b->minx, &btree[i].b->miny, &btree[i].b->minz, &btree[i].b->minv);
    if (nr != 4) panic(HERE);

    nr = fscanf(f,"  (%lf,%lf,%lf,%lf)\n", &btree[i].b->maxx, &btree[i].b->maxy, &btree[i].b->maxz, &btree[i].b->maxv);
    if (nr != 4) panic(HERE);

    fscanf(f, " }\n}\n");
	if (!check_btree_fit(&btree[i])) panic(HERE);
    blist_add_box(&boxes, btree[i].b);
   }
 
 btree = &btree[root];
 fclose(f);

 return 1;
}

void assoc_idx_tree(BTree* t, void** pi)
{
 if (t->r) assoc_idx_tree(t->r, pi);
 if (t->l) assoc_idx_tree(t->l, pi);
 pi[g_idx] = t;
 g_idx ++;
}

void create_assoc_table(void** pi)
{
 BTree* tmp;
 g_idx = 0;
 tmp = btree;
 assoc_idx_tree(tmp, pi);
}

void travel_tree(BTree* t)
{
 if (t->r) travel_tree(t->r);
 if (t->l) travel_tree(t->l);
 g_idx ++;
}

void count_indices()
{
 BTree* tmp;
 g_idx = 0;
 tmp = btree;
 travel_tree(tmp);
}

void save_box_to_file(FILE* f, Box* b)
{
 fprintf(f, " B:\n {\n");
 if (b->t1) fprintf(f, "  T1: %d\n", b->t1->idx);
 else fprintf(f, "  T1: -1\n");
 if (b->t2) fprintf(f, "  T2: %d\n", b->t2->idx);
 else fprintf(f, "  T2: -1\n");
 fprintf(f, "  (%1.12lf,%1.12lf,%1.12lf,%1.12lf)\n  (%1.12lf,%1.12lf,%1.12lf,%1.12lf)\n",
	 b->minx,b->miny,b->minz,b->minv,b->maxx,b->maxy,b->maxz,b->maxv);
 fprintf(f, " }\n");
}

void write_btree_to_file(BTree* t, FILE* f, void** p_idx, int n)
{
 int i;

 if (t->r) write_btree_to_file(t->r, f, p_idx, n);
 if (t->l) write_btree_to_file(t->l, f, p_idx, n);

 fprintf(f, "BT: %d\n{\n", g_idx);
 if (t->r)
   {
    for (i=0;i<n;i++) if (p_idx[i] == (void*)t->r)
      {
       fprintf(f, " R: %d\n", i);
       goto r_done;
      }
    panic(HERE);
   }
 else fprintf(f, " R: -1\n");
 r_done:
 if (t->l)
   {
    for (i=0;i<n;i++) if (p_idx[i] == (void*)t->l)
      {
       fprintf(f, " L: %d\n", i);
       goto l_done;
      }
    panic(HERE);
   }
 else fprintf(f, " L: -1\n");
 l_done:
 save_box_to_file(f, t->b);
 fprintf(f, "}\n");
 g_idx ++;
}

void save_preprocessed(char* scenef)
{
 void **p_idx;
 int n,i;
 FILE* f;
 char fn[2048];
 BTree* tmp;
 
 count_indices();
 n = g_idx;
 
 p_idx = (void**)malloc(n*sizeof(void*));
 create_assoc_table(p_idx);
 strcpy(fn, scenef);
 if (!strstr(fn, ".")) strcat(fn, ".btree");
 else
   {
    i = 0;
    while (fn[i] != '.') i++;
    fn[i] = 0;
    strcat(fn, ".btree");
   }
 
 f = fopen(fn, "w");
 if (!f) panic(HERE);

 tmp = btree;
 g_idx = 0;
 fprintf(f, "nElem: %d\n", n);

 //FIXME: terribly slow saving, but currently left over
 for (i=0;i<n;i++) if (p_idx[i] == (void*)btree)
   {
    fprintf(f,"R: %d\n", i);
    break;
   }
 write_btree_to_file(tmp, f, p_idx, n);
 fclose(f);
}

void free_texture(Texture* t)
{
 free((void*)(t->pixels));
}

void free_scene()
{
 int i;
 free((void*)g_ts);
 for (i=0;i<nTex;i++) free_texture(&texture[i]);
 free(texture);
}

void set_color(Screen* s, int x, int y, int z, int r, int g, int b)
{
 if (x >= 0 && x < s->x && y >= 0 && y < s->y && z >= 0 && z < s->z)
   {
    s->pixels[3*(s->z * s->y * x + s->y * y + z)]     = r;
    s->pixels[3*(s->z * s->y * x + s->y * y + z) + 1] = g;
    s->pixels[3*(s->z * s->y * x + s->y * y + z) + 2] = b;
   }
 else panic(HERE);
}

void color_backgroundr(Ray* ray, int* r, int* g, int* b)
{
 Texture* t;
 Vector v;
 REAL x,y,z,w;				/* normalized ray vector */
 REAL tx, ty, tz;			/* final texcoords */
 REAL xyz_len, xyz_x, xyz_y, xyz_z;	/* relative to cast 4D -> 3D(x,y,z),w*/
 REAL xz_len, xz_x, xz_z, xz_a;		/* relative to cast 3D -> 2D(x,z),y */


 t = &texture[0];
 memcpy(&v, &ray->d, sizeof(Vector));
 normalize(&v);

 x = v.x;
 y = v.y;
 z = v.z;
 w = v.v;
 /* (x,y,z,w) original vector */
 /* now cast into 3D and w-part */

 xyz_len = sqrt(x*x + y*y + z*z);
 if (xyz_len <= 0.) xyz_len = 1e-9;
 /* now we have length of 4D vector cast into 3D(x,y,z) */
 tz = 1. - (acos( w ) / PI);
 /* tz is "height" in 4th dimension */
 /* and there is 3rd coord of a texture */
 xyz_x = x / xyz_len;
 xyz_y = y / xyz_len;
 xyz_z = z / xyz_len;
 /* there are coords of casted 3D vector, using it we'll */
 /* get remaining tx and ty coords */

 xz_len = sqrt( xyz_x*xyz_x + xyz_z*xyz_z );
 if (xz_len <= 0.) xz_len = 1e-9;
 /* there we cast our 3D vector (after cast from 4D) into 2D on x,z plane */
 /* leaving y to compute the 2nd text coord */
 ty = 1. - (acos( xyz_y ) / PI);
 /* ty is a second text coord computed from 3D projection */

 xz_x = xyz_x / xz_len;
 xz_z = xyz_z / xz_len;
 /* there is z length cast onto 0Z axis */
 xz_a = acos( xz_z );
 /* two directions of hyperspherical texture */
 if (xyz_x < 0.) xz_a *= -1.;
 tx = (xz_a + PI) / (2. * PI);
 /* tx is first texture coord - such should be repetable in 3D texture */
 /* it's creating hyperspherical texture coords with tx,ty,tz */
 /* only X in tex should be circular, Y and Z are cylindrical */

	
 if (tx <= 0.) 
	{
	    tx = 0.000000001;
	}
	
 if (ty <= 0.) 
	{
	    ty = 0.000000001;
	}
	
 if (tz <= 0.) 
	{
	    tz = 0.000000001;
	}

 if (tx >= 1.) 
	{
	    tx = 0.999999999;
	}

 if (ty >= 1.) 
	{
	    ty = 0.999999999;
	}

 if (tz >= 1.) 
	{
	    tz = 0.999999999;
	}
	
 get_texturei(t, tx, ty, tz, r, g, b);
}

void calculate_color(Screen* s, Ray* r, Triangle* t, int x, int y, int z)
{
 Vertex v;
 Vector n;
 Material Rm, Gm, Bm;
 TexCoord tC;
 REAL rr,gg,bb;
 int idx, ir,ig,ib;

 if (!intersection(t, r, &v, &idx, itree))
   {
       color_backgroundr(r, &ir, &ig, &ib);
       set_color(s, x, y, z, ir, ig, ib);
   }
 else
   {
	   get_normal(&n, &t[idx], &v, &Rm, &Gm, &Bm, &tC);
       r->c = GREEN;    
       gg = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1., itree);
       
       
       r->c = RED;
       rr = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1., itree);

       r->c = BLUE;
	   bb = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1., itree);

	   //FIXME: wtf?
	   /*
       if (rr < 0.0036 && gg < 0.0036 && bb < 0.036 && y > 0) 
	   {
            ir = get_red(s, x, y-1);
            ig = get_green(s, x, y-1);
            ib = get_blue(s, x, y-1);
	    set_color(s, x, y, ir, ig, ib);
	   } else */
	   set_color(s, x, y, z, (int)(rr*255.), (int)(gg*255.), (int)(bb*255.));
   }
}

// get_normal			= get_normal_old
// intersection			= intersection_new
// intersection_itree	= intersection_itree_spec_shadow
void raytrace()
{
 Ray r;
 Screen* s;
 int i,j,k;

 s = &screen;

 r.P.x = observer.x;
 r.P.y = observer.y;
 r.P.z = observer.z;
 r.P.v = observer.v;
 r.r = 0;
 r.d.z = (s->x + s->y + s->z) * lookz;	

 itree = alloc_itree();
 for (i=0;i<s->x;i++)
   {
    r.d.x = i - s->x/2.;
    for (j=0;j<s->y;j++)
      {
       r.d.y = j - s->y/2.;
	   for (k=0;k<s->z;k++)
		 {
		  r.d.z = k - s->z/2.;
		  calculate_color(s, &r, g_ts, i, j, k);
		 }       
      }
   }

 free_scene();
}


int rt4d(char* scenefile)
{
	int prep_loaded;
	load_scene(scenefile);
	if (!(prep_loaded = load_preprocessed(scenefile))) 
	{
		preprocess_scene(scenefile);
		save_preprocessed(scenefile);
	}

	raytrace();
	return 0;
}

int main(int lb, char** par)
{
	if (lb < 2)
	{
		printf("%s scenefile.4d\n", par[0]);
		return 1;
	}
	return rt4d(par[1]);
}
