#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define HERE __FILE__,__LINE__
#define BACK_R 0x4F
#define BACK_G 0x1F
#define BACK_B 0x7F
#define MAX_REC 6
#define distorber1 (double)(1./17.)
#define distorber2 (double)(1./19.)
enum { RED = 0x0606, GREEN, BLUE};

typedef struct _Vertex
{
 double x,y,z;
} Vertex;
typedef struct _Vertex Vector;

typedef struct _Material
{
 double c,s,t;
} Material;

typedef struct _TexCoord
{
 double x,y;
} TexCoord;

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
typedef struct _Screen Texture;

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
 double mUR,mDR;
 double mUG,mDG;
 double mUB,mDB;
 double ca;
 int faces;
 int idx;
 int tid;
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
 double minx,miny,minz;
 double maxx,maxy,maxz;
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

typedef struct _BList
{
 struct _BList *next, *prev;
 Box b;
} BList;
char textureDir[512];
Vertex light, observer;	/* FIXME: more lights ? */
Screen screen;
int nTriangles;
FILE* dbg;
Texture* texture;
int nTex;
double step;
BList* boxes;
BTree* btree;
int glob_sol;
int proc_tr;
int proc_tr2;
int (*intersection)(Triangle*, Ray*, Vertex*, int*);
void (*get_normal)(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);
/*int nmalloc,nfree;*/
void wrt_bmp(Screen* s, char* out_f);
int intersection_new(Triangle*, Ray*, Vertex*, int*);
int intersection_old(Triangle*, Ray*, Vertex*, int*);
void get_normal_new(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);
void get_normal_old(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);

void panic(char* f, int l)
{
 printf("RayTracer Engine Panic: file: %s, line: %d\n", f,l);
 wrt_bmp(&screen, "panic.bmp");
 exit(1);
}


void init_bmp(BMPTag* b)	/*FIXME: add option to save JPEG */
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


int read_texcoord(FILE* f, TexCoord* t)
{
 int nr;
 if (!f || !t) panic(HERE);
 nr = fscanf(f, "TexCoord: (%lf,%lf)", &t->x, &t->y);
/* printf("%d %f,%f\n", nr, t->x, t->y);*/
 if (nr != 2) return 0;
 if (t->x < 0. || t->x > 1.) return 0;
 if (t->y < 0. || t->y > 1.) return 0;
 return 1;
}

void normalize(Vector*);

void normalize_color_factors(Triangle* t)
{
 double f;
 /*f = t->cr + t->sr + t->tr;
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
 t->tb /= f;*/
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
/* printf("normal: (%f,%f,%f)\n", nx, ny, nz);*/
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

int read_triangle(FILE* f, Triangle* t, int idx)
{
 int nr,n,ti;
 char str[128];
 if (!f || !t) panic(HERE);
 nr = fscanf(f, "Triangle: %d\n{\n", &n);
/* printf("nr = %d\n", nr);*/
 if (nr != 1) return 0;
 if (n < 0) return 0;
 t->idx = idx;
 fscanf(f, " a: ");
 if (!read_vertex(f, &t->a)) panic(HERE);
 fscanf(f,"\n");
 t->a.x += distorber1;
 t->a.y += distorber2;
 fscanf(f, " b: ");
 if (!read_vertex(f, &t->b)) panic(HERE);
 fscanf(f,"\n");
 t->b.x += distorber1;
 t->b.y += distorber2;
 fscanf(f, " c: ");
 if (!read_vertex(f, &t->c)) panic(HERE);
 fscanf(f,"\n");
 t->c.x += distorber1;
 t->c.y += distorber2;
 fscanf(f, " texA: ");
 if (!read_texcoord(f, &t->ta)) panic(HERE);
 fscanf(f,"\n");
 fscanf(f, " texB: ");
 if (!read_texcoord(f, &t->tb)) panic(HERE);
 fscanf(f,"\n");
 fscanf(f, " texC: ");
 if (!read_texcoord(f, &t->tc)) panic(HERE);
 fscanf(f,"\n");
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
 nr = fscanf(f," %s RGB: (%lf,%lf,%lf)\n", str, &t->mra.t, &t->mga.t, &t->mba.t);
/* printf("nr = %d\n", nr);*/
 if (nr != 4) panic(HERE);
 if (!strcmp(str, "transparency:"))
   {
    nr = fscanf(f," specular: RGB: (%lf,%lf,%lf)\n", &t->mra.s, &t->mga.s, &t->mba.s);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," diffuse: RGB: (%lf,%lf,%lf)\n", &t->mra.c, &t->mga.c, &t->mba.c);
    if (nr != 3) panic(HERE);
    t->mrb.t = t->mra.t;
    t->mrc.t = t->mra.t;
    t->mgb.t = t->mga.t;
    t->mgc.t = t->mga.t;
    t->mbb.t = t->mba.t;
    t->mbc.t = t->mba.t;
    t->mrb.s = t->mra.s;
    t->mrc.s = t->mra.s;
    t->mgb.s = t->mga.s;
    t->mgc.s = t->mga.s;
    t->mbb.s = t->mba.s;
    t->mbc.s = t->mba.s;
    t->mrb.c = t->mra.c;
    t->mrc.c = t->mra.c;
    t->mgb.c = t->mga.c;
    t->mgc.c = t->mga.c;
    t->mbb.c = t->mba.c;
    t->mbc.c = t->mba.c;
   }
 else if (!strcmp(str, "transparencyA:"))
   {
    nr = fscanf(f," specularA: RGB: (%lf,%lf,%lf)\n", &t->mra.s, &t->mga.s, &t->mba.s);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," diffuseA: RGB: (%lf,%lf,%lf)\n", &t->mra.c, &t->mga.c, &t->mba.c);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," transparencyB: RGB: (%lf,%lf,%lf)\n", &t->mrb.t, &t->mgb.t, &t->mbb.t);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," specularB: RGB: (%lf,%lf,%lf)\n", &t->mrb.s, &t->mgb.s, &t->mbb.s);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," diffuseB: RGB: (%lf,%lf,%lf)\n", &t->mrb.c, &t->mgb.c, &t->mbb.c);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," transparencyC: RGB: (%lf,%lf,%lf)\n", &t->mrc.t, &t->mgc.t, &t->mbc.t);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," specularC: RGB: (%lf,%lf,%lf)\n", &t->mrc.s, &t->mgc.s, &t->mbc.s);
    if (nr != 3) panic(HERE);
    nr = fscanf(f," diffuseC: RGB: (%lf,%lf,%lf)\n", &t->mrc.c, &t->mgc.c, &t->mbc.c);
    if (nr != 3) panic(HERE);
   }
 else panic(HERE);
 nr = fscanf(f," surface: ABCD: (%lf,%lf,%lf,%lf)\n", &t->s.A, &t->s.B, &t->s.C, &t->s.D);
 if (nr != 4) panic(HERE);
 normalize_color_factors(t);
 if (t->s.A == 0. && t->s.B == 0. && t->s.C == 0. && t->s.D == 0.) compute_surface(t);
 nr = fscanf(f, " %s (%lf,%lf)\n", str, &t->mUR, &t->mDR);
 if (nr != 3) panic(HERE);
 if (!strcmp(str, "transparencyFact:"))
   {
    t->mUG = t->mUR;
    t->mUB = t->mUR;
    t->mDG = t->mDR;
    t->mDB = t->mDR;
   }
 else if (!strcmp(str, "transparencyFactR:"))
   {
    nr = fscanf(f, " transparencyFactG: (%lf,%lf)\n", &t->mUG, &t->mDG);
    if (nr != 2) panic(HERE);
    nr = fscanf(f, " transparencyFactB: (%lf,%lf)\n", &t->mUB, &t->mDB);
    if (nr != 2) panic(HERE);
   }
 else panic(HERE);
 if (t->mUR < .05 || t->mUR > 20.) panic(HERE);
 if (t->mDR < .05 || t->mDR > 20.) panic(HERE);
 if (t->mUG < .05 || t->mUG > 20.) panic(HERE);
 if (t->mDG < .05 || t->mDG > 20.) panic(HERE);
 if (t->mUB < .05 || t->mUB > 20.) panic(HERE);
 if (t->mDB < .05 || t->mDB > 20.) panic(HERE);
 nr = fscanf(f, " specularFact: %lf\n", &t->ca);
 if (nr != 1) panic(HERE);
 if (t->ca < -1.) panic(HERE);
 nr = fscanf(f, " faces: %d\n", &t->faces);
 if (nr != 1) panic(HERE);
 if (t->faces != 1 && t->faces != 2) panic(HERE);
 nr = fscanf(f, " texture: %d\n", &ti);
 if (nr != 1) panic(HERE);
 if (n < 0) panic(HERE);
 if (ti > 0 && ti <= nTex) { t->tex = &texture[ti-1]; t->tid = ti; }
 else { t->tex = NULL; t->tid = 0; }
 fscanf(f,"}\n");
/* if (t->tex) printf("texture = %p(%d,%d)\n", t->tex->pixels, t->tex->x, t->tex->y);*/
 return 1;
}

void set_color(Screen*, int, int, int, int, int);

void create_texture(Texture* t, char* dir, int n)
{
 BMPTag bm_handle;
 FILE* plik;
 char fn[512];
 int i,j;
 char r,g,b,m;
 init_bmp(&bm_handle);
 sprintf(fn, "./%s/%d.bmp", dir, n+1);
 printf("Opening texture file: %s\n", fn);
 plik = fopen(fn,"r");
 if (!plik) { printf("Error opening BMP: %s\n", fn); panic(HERE); }
 i = fscanf(plik,"%c%c",&b,&m);
 if (i != 2) panic(HERE);
 if (b != 'B' || m != 'M') panic(HERE);
 fread(&bm_handle.fsize,4,1,plik);
 fread(&bm_handle.dummy,4,1,plik);
 fread(&bm_handle.offset,4,1,plik);
 fread(&bm_handle.dummy2,4,1,plik);
 fread(&bm_handle.bm_x,4,1,plik);
 fread(&bm_handle.bm_y,4,1,plik);
 fread(&bm_handle.planes,2,1,plik);
 fread(&bm_handle.bpp,2,1,plik);
 if (bm_handle.bpp != 24) panic(HERE);
 fseek(plik,bm_handle.offset,SEEK_SET);
 t->pixels = (unsigned char*)malloc(3*bm_handle.bm_y*bm_handle.bm_x*sizeof(unsigned char));
 t->x = bm_handle.bm_x;
 t->y = bm_handle.bm_y;
 for (i=0;i<bm_handle.bm_y;i++)  for (j=0;j<bm_handle.bm_x;j++)
    {
     fscanf(plik,"%c%c%c", &b,&g,&r);
     set_color(t, j, i, r, g, b);
    }
/* printf("Created texture %d from file %s: %p(%d,%d)\n", n, fn, t->pixels, t->x, t->y);*/
 fclose(plik);
}


void free_texture(Texture* t)
{
 free((void*)(t->pixels));
}


void create_textures(char* dir)
{
 int i;
 if (!strcmp(dir, "")) panic(HERE);
 texture = (Texture*)malloc(nTex*sizeof(Texture));
 if (!texture) panic(HERE);
 for (i=0;i<nTex;i++) create_texture(&texture[i], dir, i);
}


int is_binary(FILE* f)
{
 char str[64];
 fread(str,4,1,f);
 fseek(f, 0, SEEK_SET);
 if (!strcmp(str,"BINT")) return 1;
 else return 0;
}


void load_binary_scene(FILE* f, Triangle** ts)
{
 int tlen,ti,i;
 int sx,sy;
 char texDir[256];
 Triangle* t;
/* int nTex;*/
 fread(texDir, 4, 1, f);
 strcpy(texDir, "");
 fread(&sx,sizeof(int),1,f);
 fread(&sy,sizeof(int),1,f);
 if (sx <= 0 || sy <= 0) panic(HERE);
 init_screen(&screen, sx, sy);
 fread(&observer.x,sizeof(double),1,f);
 fread(&observer.y,sizeof(double),1,f);
 fread(&observer.z,sizeof(double),1,f);
 fread(&light.x,sizeof(double),1,f);
 fread(&light.y,sizeof(double),1,f);
 fread(&light.z,sizeof(double),1,f);
 fread(&tlen,sizeof(int),1,f);
 fread(texDir,tlen,1,f);
 fread(&nTex,sizeof(int),1,f);
 if (nTex < 0) panic(HERE);
 create_textures(texDir);
 fread(&nTriangles,sizeof(int),1,f);
 if (nTriangles <= 0) panic(HERE);
 *ts = (Triangle*)malloc(nTriangles*sizeof(Triangle));
 t = *ts;
 for (i=0;i<nTriangles;i++)
   {
    fread(&ti,sizeof(int),1,f);
    if (ti < 0) panic(HERE);
    t[i].idx = i;
    fread(&t[i].a.x,sizeof(double),1,f);
    fread(&t[i].a.y,sizeof(double),1,f);
    fread(&t[i].a.z,sizeof(double),1,f);
    t[i].a.x += distorber1;
    t[i].a.y += distorber2;
    fread(&t[i].b.x,sizeof(double),1,f);
    fread(&t[i].b.y,sizeof(double),1,f);
    fread(&t[i].b.z,sizeof(double),1,f);
    t[i].b.x += distorber1;
    t[i].b.y += distorber2;
    fread(&t[i].c.x,sizeof(double),1,f);
    fread(&t[i].c.y,sizeof(double),1,f);
    fread(&t[i].c.z,sizeof(double),1,f);
    t[i].c.x += distorber1;
    t[i].c.y += distorber2;
    fread(&t[i].ta.x,sizeof(double),1,f);
    fread(&t[i].ta.y,sizeof(double),1,f);
    fread(&t[i].tb.x,sizeof(double),1,f);
    fread(&t[i].tb.y,sizeof(double),1,f);
    fread(&t[i].tc.x,sizeof(double),1,f);
    fread(&t[i].tc.y,sizeof(double),1,f);
    fread(&t[i].na.x,sizeof(double),1,f);
    fread(&t[i].na.y,sizeof(double),1,f);
    fread(&t[i].na.z,sizeof(double),1,f);
    fread(&t[i].nb.x,sizeof(double),1,f);
    fread(&t[i].nb.y,sizeof(double),1,f);
    fread(&t[i].nb.z,sizeof(double),1,f);
    fread(&t[i].nc.x,sizeof(double),1,f);
    fread(&t[i].nc.y,sizeof(double),1,f);
    fread(&t[i].nc.z,sizeof(double),1,f);
    if (length(&t[i].na) <= 1e-9 || length(&t[i].nb) <= 1e-9 || length(&t[i].nc) <= 1e-9)
      {
       compute_normals(&t[i]);
      }
    else
      {
       normalize(&t[i].na);
       normalize(&t[i].nb);
       normalize(&t[i].nc);
      }
    fread(&t[i].mra.t,sizeof(double),1,f);
    fread(&t[i].mga.t,sizeof(double),1,f);
    fread(&t[i].mba.t,sizeof(double),1,f);
    fread(&t[i].mra.s,sizeof(double),1,f);
    fread(&t[i].mga.s,sizeof(double),1,f);
    fread(&t[i].mba.s,sizeof(double),1,f);
    fread(&t[i].mra.c,sizeof(double),1,f);
    fread(&t[i].mga.c,sizeof(double),1,f);
    fread(&t[i].mba.c,sizeof(double),1,f);
    fread(&t[i].mrb.t,sizeof(double),1,f);
    fread(&t[i].mgb.t,sizeof(double),1,f);
    fread(&t[i].mbb.t,sizeof(double),1,f);
    fread(&t[i].mrb.s,sizeof(double),1,f);
    fread(&t[i].mgb.s,sizeof(double),1,f);
    fread(&t[i].mbb.s,sizeof(double),1,f);
    fread(&t[i].mrb.c,sizeof(double),1,f);
    fread(&t[i].mgb.c,sizeof(double),1,f);
    fread(&t[i].mbb.c,sizeof(double),1,f);
    fread(&t[i].mrc.t,sizeof(double),1,f);
    fread(&t[i].mgc.t,sizeof(double),1,f);
    fread(&t[i].mbc.t,sizeof(double),1,f);
    fread(&t[i].mrc.s,sizeof(double),1,f);
    fread(&t[i].mgc.s,sizeof(double),1,f);
    fread(&t[i].mbc.s,sizeof(double),1,f);
    fread(&t[i].mrc.c,sizeof(double),1,f);
    fread(&t[i].mgc.c,sizeof(double),1,f);
    fread(&t[i].mbc.c,sizeof(double),1,f);
    fread(&t[i].s.A,sizeof(double),1,f);
    fread(&t[i].s.B,sizeof(double),1,f);
    fread(&t[i].s.C,sizeof(double),1,f);
    fread(&t[i].s.D,sizeof(double),1,f);
    normalize_color_factors(t);
    if (t[i].s.A == 0. && t[i].s.B == 0. && t[i].s.C == 0. && t[i].s.D == 0.) compute_surface(&t[i]);
    fread(&t[i].mUR,sizeof(double),1,f);
    fread(&t[i].mDR,sizeof(double),1,f);
    fread(&t[i].mUG,sizeof(double),1,f);
    fread(&t[i].mDG,sizeof(double),1,f);
    fread(&t[i].mUB,sizeof(double),1,f);
    fread(&t[i].mDB,sizeof(double),1,f);
    if (t[i].mUR < .05 || t[i].mUR > 20.) panic(HERE);
    if (t[i].mDR < .05 || t[i].mDR > 20.) panic(HERE);
    if (t[i].mUG < .05 || t[i].mUG > 20.) panic(HERE);
    if (t[i].mDG < .05 || t[i].mDG > 20.) panic(HERE);
    if (t[i].mUB < .05 || t[i].mUB > 20.) panic(HERE);
    if (t[i].mDB < .05 || t[i].mDB > 20.) panic(HERE);
    fread(&t[i].ca,sizeof(double),1,f);
    if (t[i].ca < -1.) panic(HERE);
    fread(&t[i].faces,sizeof(double),1,f);
    if (t[i].faces != 1 && t[i].faces != 2) panic(HERE);
    fread(&ti,sizeof(double),1,f);
    if (ti > 0 && ti <= nTex) t[i].tex = &texture[ti-1];
    else t[i].tex = NULL;
    t[i].tid = ti;
   }
 strcpy(textureDir, texDir);
 fclose(f);
}


void save_binary_scene(char* name, Triangle* t)
{
 int tlen,i;
 char out[512];
 char sig[5];
 FILE* f;
 strcpy(out,name);
 strcat(out,".bin");
 f = fopen(out, "w");
 if (!f) panic(HERE);
/* int nTex;*/
 strcpy(sig, "BINT");
 fwrite(sig,4,1,f);
 fwrite(&screen.x,sizeof(int),1,f);
 fwrite(&screen.y,sizeof(int),1,f);
 fwrite(&observer.x,sizeof(double),1,f);
 fwrite(&observer.y,sizeof(double),1,f);
 fwrite(&observer.z,sizeof(double),1,f);
 fwrite(&light.x,sizeof(double),1,f);
 fwrite(&light.y,sizeof(double),1,f);
 fwrite(&light.z,sizeof(double),1,f);
 tlen = strlen(textureDir);
/* printf("textureDir: %s, len = %d\n", textureDir, tlen);*/
 fwrite(&tlen,sizeof(int),1,f);
 fwrite(textureDir,tlen,1,f);
 fwrite(&nTex,sizeof(int),1,f);
 fwrite(&nTriangles,sizeof(int),1,f);
 for (i=0;i<nTriangles;i++)
   {
    fwrite(&t[i].idx,sizeof(int),1,f);
    t[i].a.x -= distorber1;
    t[i].a.y -= distorber2;
    fwrite(&t[i].a.x,sizeof(double),1,f);
    fwrite(&t[i].a.y,sizeof(double),1,f);
    fwrite(&t[i].a.z,sizeof(double),1,f);
    t[i].a.x += distorber1;
    t[i].a.y += distorber2;
    t[i].b.x -= distorber1;
    t[i].b.y -= distorber2;
    fwrite(&t[i].b.x,sizeof(double),1,f);
    fwrite(&t[i].b.y,sizeof(double),1,f);
    fwrite(&t[i].b.z,sizeof(double),1,f);
    t[i].b.x += distorber1;
    t[i].b.y += distorber2;
    t[i].c.x -= distorber1;
    t[i].c.y -= distorber2;
    fwrite(&t[i].c.x,sizeof(double),1,f);
    fwrite(&t[i].c.y,sizeof(double),1,f);
    fwrite(&t[i].c.z,sizeof(double),1,f);
    t[i].c.x += distorber1;
    t[i].c.y += distorber2;
    fwrite(&t[i].ta.x,sizeof(double),1,f);
    fwrite(&t[i].ta.y,sizeof(double),1,f);
    fwrite(&t[i].tb.x,sizeof(double),1,f);
    fwrite(&t[i].tb.y,sizeof(double),1,f);
    fwrite(&t[i].tc.x,sizeof(double),1,f);
    fwrite(&t[i].tc.y,sizeof(double),1,f);
    fwrite(&t[i].na.x,sizeof(double),1,f);
    fwrite(&t[i].na.y,sizeof(double),1,f);
    fwrite(&t[i].na.z,sizeof(double),1,f);
    fwrite(&t[i].nb.x,sizeof(double),1,f);
    fwrite(&t[i].nb.y,sizeof(double),1,f);
    fwrite(&t[i].nb.z,sizeof(double),1,f);
    fwrite(&t[i].nc.x,sizeof(double),1,f);
    fwrite(&t[i].nc.y,sizeof(double),1,f);
    fwrite(&t[i].nc.z,sizeof(double),1,f);
    fwrite(&t[i].mra.t,sizeof(double),1,f);
    fwrite(&t[i].mga.t,sizeof(double),1,f);
    fwrite(&t[i].mba.t,sizeof(double),1,f);
    fwrite(&t[i].mra.s,sizeof(double),1,f);
    fwrite(&t[i].mga.s,sizeof(double),1,f);
    fwrite(&t[i].mba.s,sizeof(double),1,f);
    fwrite(&t[i].mra.c,sizeof(double),1,f);
    fwrite(&t[i].mga.c,sizeof(double),1,f);
    fwrite(&t[i].mba.c,sizeof(double),1,f);
    fwrite(&t[i].mrb.t,sizeof(double),1,f);
    fwrite(&t[i].mgb.t,sizeof(double),1,f);
    fwrite(&t[i].mbb.t,sizeof(double),1,f);
    fwrite(&t[i].mrb.s,sizeof(double),1,f);
    fwrite(&t[i].mgb.s,sizeof(double),1,f);
    fwrite(&t[i].mbb.s,sizeof(double),1,f);
    fwrite(&t[i].mrb.c,sizeof(double),1,f);
    fwrite(&t[i].mgb.c,sizeof(double),1,f);
    fwrite(&t[i].mbb.c,sizeof(double),1,f);
    fwrite(&t[i].mrc.t,sizeof(double),1,f);
    fwrite(&t[i].mgc.t,sizeof(double),1,f);
    fwrite(&t[i].mbc.t,sizeof(double),1,f);
    fwrite(&t[i].mrc.s,sizeof(double),1,f);
    fwrite(&t[i].mgc.s,sizeof(double),1,f);
    fwrite(&t[i].mbc.s,sizeof(double),1,f);
    fwrite(&t[i].mrc.c,sizeof(double),1,f);
    fwrite(&t[i].mgc.c,sizeof(double),1,f);
    fwrite(&t[i].mbc.c,sizeof(double),1,f);
    fwrite(&t[i].s.A,sizeof(double),1,f);
    fwrite(&t[i].s.B,sizeof(double),1,f);
    fwrite(&t[i].s.C,sizeof(double),1,f);
    fwrite(&t[i].s.D,sizeof(double),1,f);
    fwrite(&t[i].mUR,sizeof(double),1,f);
    fwrite(&t[i].mDR,sizeof(double),1,f);
    fwrite(&t[i].mUG,sizeof(double),1,f);
    fwrite(&t[i].mDG,sizeof(double),1,f);
    fwrite(&t[i].mUB,sizeof(double),1,f);
    fwrite(&t[i].mDB,sizeof(double),1,f);
    fwrite(&t[i].ca,sizeof(double),1,f);
    fwrite(&t[i].faces,sizeof(double),1,f);
    fwrite(&t[i].tid,sizeof(double),1,f);
   }
 fclose(f);
}


void load_scene(Triangle** ts, char* scenefile)
{
 FILE* f;
 char texDir[512];
 int nr,sx,sy,n,i;
 f = fopen(scenefile,"r");
 if (!f) panic(HERE);
 if (is_binary(f)) { load_binary_scene(f, ts); return; }
 nr = fscanf(f, "Screen: (%d,%d)\n", &sx, &sy);
 if (nr != 2) panic(HERE);
 if (sx <= 0 || sy <= 0) panic(HERE);
 init_screen(&screen, sx, sy);
 fscanf(f, "Observer: ");
 if (!read_vertex(f, &observer)) panic(HERE);
 fscanf(f,"\n");
 fscanf(f, "Light: ");
 if (!read_vertex(f, &light)) panic(HERE);
 fscanf(f,"\n");
 nr = fscanf(f, "TexDirectory: %s\n", texDir);
 if (nr != 1) panic(HERE);
 nr = fscanf(f, "NumTextures: %d\n", &nTex);
 if (nr != 1) panic(HERE);
 if (nTex < 0) panic(HERE);
 create_textures(texDir);
 nr = fscanf(f, "nTriangles: %d\n", &n);
 if (nr != 1) panic(HERE);
 if (n <= 0) panic(HERE);
 *ts = (Triangle*)malloc(n*sizeof(Triangle));
 for (i=0;i<n;i++)
   {
    if (!read_triangle(f, &((*ts)[i]), i)) panic(HERE);
/*    printf("read %d triangles\n", i);*/
   }
 nTriangles = n;
 fclose(f);
 strcpy(textureDir, texDir);
 save_binary_scene(scenefile, *ts);
}


void free_scene(Triangle** ts)
{
 int i;
 free((void*)(*ts));
 for (i=0;i<nTex;i++) free_texture(&texture[i]);
 free(texture);
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
 return (fabs(x-y) < eps);
}


double compute_angle(Vertex* a, Vertex* b, Vertex* c)
{
 register double x1,x2,y1,y2,z1,z2,len;
 x1 = b->x - a->x;
 x2 = c->x - a->x;
 y1 = b->y - a->y;
 y2 = c->y - a->y;
 z1 = b->z - a->z;
 z2 = c->z - a->z;
 len = sqrt(x1*x1+y1*y1+z1*z1);
 if (nearly_equal(len, 0., 1e-7)) return 4.;
 x1 /= len;
 y1 /= len;
 z1 /=len;
 len = sqrt(x2*x2+y2*y2+z2*z2);
 if (nearly_equal(len, 0., 1e-7)) return 4.;
 x2 /= len;
 y2 /= len;
 z2 /= len;
 return (x1*x2+y1*y2+z1*z2);
}


/*double compute_angle2(Vertex* a, Vertex* b, Vertex* c)
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
}*/

int vertex_in_triangle(Triangle* t, Vertex* w)
{
 double kat;
 kat  = compute_angle(w, &t->a, &t->b);
 kat += compute_angle(w, &t->b, &t->c);
 kat += compute_angle(w, &t->c, &t->a);
/* printf("kat = %f\n", kat);*/
 if (kat < -0.99999)		/* FIXME: is this test OK? */
   {
    return 1;
   }
 else return 0;
}


void blist_add(BList** head, double ix, double iy, double iz, double ax, double ay, double az, Triangle* t1, Triangle* t2)
{
 BList* temp;
 if (*head == NULL)
   {
/*     nmalloc++;*/
    *head = (BList*)(malloc(sizeof(BList)));
    /* FIXME: test disabled for speed\n"); */
/*    if (!(*head)) panic(HERE);*/
    (*head)->b.t1 = t1;
    (*head)->b.t2 = t2;
    (*head)->b.minx = ix;
    (*head)->b.miny = iy;
    (*head)->b.minz = iz;
    (*head)->b.maxx = ax;
    (*head)->b.maxy = ay;
    (*head)->b.maxz = az;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }
/*  nmalloc++;*/
  temp = (BList*)(malloc(sizeof(BList)));
/*  if (!temp) panic(HERE);*/
  temp->b.t1 = t1;
  temp->b.t2 = t2;
  temp->b.minx = ix;
  temp->b.miny = iy;
  temp->b.minz = iz;
  temp->b.maxx = ax;
  temp->b.maxy = ay;
  temp->b.maxz = az;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
}


void ptrlist_add(PtrList** head, Triangle* t)
{
 PtrList* temp;
 if (*head == NULL)
   {
/*     nmalloc++;*/
    *head = (PtrList*)(malloc(sizeof(PtrList)));
    /* FIXME: test disabled for speed\n"); */
/*    if (!(*head)) panic(HERE);*/
    (*head)->ptr = (void*)t;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }
/*  nmalloc++;*/
  temp = (PtrList*)(malloc(sizeof(PtrList)));
/*  if (!temp) panic(HERE);*/
  temp->ptr = (void*)t;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
}


void vlist_add(VList** head, Vertex* v, int idx)
{
 VList* temp;
 if (*head == NULL)
   {
/*     nmalloc++;*/
    *head = (VList*)(malloc(sizeof(VList)));
    /* FIXME: test disabled for speed\n"); */
/*    if (!(*head)) panic(HERE);*/
    (*head)->P.x = v->x;
    (*head)->P.y = v->y;
    (*head)->P.z = v->z;
    (*head)->idx = idx;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }
/*  nmalloc++;*/
  temp = (VList*)(malloc(sizeof(VList)));
/*  if (!temp) panic(HERE);*/
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
 /* FIXME: some tests skipped for speed.\n"); */
/* if (!ret) panic(HERE);*/
 s = &tr->s;
 tmp = ( s->A*r->d.x +  s->B*r->d.y + s->C*r->d.z);
/* printf("ABCD=(%f,%f,%f,%f) tmp = %f\n", s->A, s->B, s->C, s->D, tmp);*/
 if (!nearly_equal(tmp, 0., 1e-8)) /* FIXME: was 1e-9 */
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
 else return 0;
}

double distance(Vertex*, Vertex*);
double pseudo_distance(Vertex*, Vertex*);

void ptrlist_free(PtrList** head)
{
 PtrList* temp;
 if (!head) return;
/* printf("freeying\n");*/
 while (*head)
   {
    temp = *head;
/*    printf("%p -> %p\n", temp, temp->next);*/
    *head = (*head)->next;
/*    nfree++;*/
    free(temp);
/*    if (*head) (*head)->prev = NULL;*/
/*    printf("freed.\n");*/
   }
/* *head = NULL;*/
}


void blist_free(BList** head)
{
 BList* temp;
 if (!head) return;
/* printf("freeying\n");*/
 while (*head)
   {
    temp = *head;
/*    printf("%p -> %p\n", temp, temp->next);*/
    *head = (*head)->next;
/*    nfree++;*/
    free(temp);
/*    if (*head) (*head)->prev = NULL;*/
/*    printf("freed.\n");*/
   }
/* *head = NULL;*/
}


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
/*    nfree++;*/
    free(temp);
/*    if (*head) (*head)->prev = NULL;*/
/*    printf("freed.\n");*/
   }
/* *head = NULL;*/
}


void intersection_btree(VList** head, Ray* r, BTree* tree)
{
 Box* b;
 register double K,x,y,z;
 Vertex ret;
/* if ((tree->l && !tree->r) || (!tree->l && tree->r)) panic(HERE);*/
 if (!tree->l && !tree->r)
 {
  proc_tr += 2;
  if (intersection_triangle(tree->b->t1, r, &ret))
      {
       vlist_add(head, &ret, tree->b->t1->idx);
       glob_sol ++;
      }
  if (tree->b->t2 && intersection_triangle(tree->b->t2, r, &ret))
      {
       vlist_add(head, &ret, tree->b->t2->idx);
       glob_sol ++;
      }
   return;
  }
  b = tree->r->b;
  K = (b->minz - r->P.z) / r->d.z;
  if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       y = r->P.y + K * r->d.y;
       if (x >= b->minx && x <= b->maxx && y >= b->miny && y <= b->maxy)
         {
          intersection_btree(head, r, tree->r);
          goto lcase;
         }
      }
    K = (b->miny - r->P.y) / r->d.y;
    if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       z = r->P.z + K * r->d.z;
       if (x >= b->minx && x <= b->maxx && z >= b->minz && z <= b->maxz)
         {
          intersection_btree(head, r, tree->r);
          goto lcase;
	 }
      }
    K = (b->minx - r->P.x) / r->d.x;
    if (K >= 1e-9)
      {
       y = r->P.y + K * r->d.y;
       z = r->P.z + K * r->d.z;
       if (y >= b->miny && y <= b->maxy && z >= b->minz && z <= b->maxz)
         {
          intersection_btree(head, r, tree->r);
          goto lcase;
	 }
      }
    K = (b->maxz - r->P.z) / r->d.z;
    if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       y = r->P.y + K * r->d.y;
       if (x >= b->minx && x <= b->maxx && y >= b->miny && y <= b->maxy)
         {
          intersection_btree(head, r, tree->r);
          goto lcase;
	 }
      }
    K = (b->maxy - r->P.y) / r->d.y;
    if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       z = r->P.z + K * r->d.z;
       if (x >= b->minx && x <= b->maxx && z >= b->minz && z <= b->maxz)
         {
          intersection_btree(head, r, tree->r);
          goto lcase;
	 }
      }
    K = (b->maxx - r->P.x) / r->d.x;
    if (K >= 1e-9)
      {
       y = r->P.y + K * r->d.y;
       z = r->P.z + K * r->d.z;
       if (y >= b->miny && y <= b->maxy && z >= b->minz && z <= b->maxz)
          intersection_btree(head, r, tree->r);
          goto lcase;
      }
 lcase:
    b = tree->l->b;
    K = (b->minz - r->P.z) / r->d.z;
    if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       y = r->P.y + K * r->d.y;
       if (x >= b->minx && x <= b->maxx && y >= b->miny && y <= b->maxy)
         {
          intersection_btree(head, r, tree->l);
          return;
	 }
      }
    K = (b->miny - r->P.y) / r->d.y;
    if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       z = r->P.z + K * r->d.z;
       if (x >= b->minx && x <= b->maxx && z >= b->minz && z <= b->maxz)
         {
          intersection_btree(head, r, tree->l);
          return;
	 }
      }
    K = (b->minx - r->P.x) / r->d.x;
    if (K >= 1e-9)
      {
       y = r->P.y + K * r->d.y;
       z = r->P.z + K * r->d.z;
       if (y >= b->miny && y <= b->maxy && z >= b->minz && z <= b->maxz)
         {
          intersection_btree(head, r, tree->l);
          return;
	 }
      }
    K = (b->maxz - r->P.z) / r->d.z;
    if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       y = r->P.y + K * r->d.y;
       if (x >= b->minx && x <= b->maxx && y >= b->miny && y <= b->maxy)
         {
          intersection_btree(head, r, tree->l);
          return;
	 }
      }
    K = (b->maxy - r->P.y) / r->d.y;
    if (K >= 1e-9)
      {
       x = r->P.x + K * r->d.x;
       z = r->P.z + K * r->d.z;
       if (x >= b->minx && x <= b->maxx && z >= b->minz && z <= b->maxz)
         {
          intersection_btree(head, r, tree->l);
          return;
	 }
      }
    K = (b->maxx - r->P.x) / r->d.x;
    if (K >= 1e-9)
      {
       y = r->P.y + K * r->d.y;
       z = r->P.z + K * r->d.z;
       if (y >= b->miny && y <= b->maxy && z >= b->minz && z <= b->maxz)
         {
          intersection_btree(head, r, tree->l);
          return;
	 }
      }
}


int intersection_new(Triangle* tr, Ray* r, Vertex* ret, int* idx)
{
 double dist,ndist;
 int ts;
 int sol;
 VList* head;
 VList* hd;
 VList* amb;
/* VList* tamb;*/
 head = NULL;
 glob_sol = 0;
/* proc_tr = 0;*/
 intersection_btree(&head, r, btree);
/* printf("%f %% (%d/%d)\n", (double)(proc_tr * 100)/(double)nTriangles, proc_tr, nTriangles);*/
 sol = glob_sol;
 proc_tr2 += nTriangles;
 if (sol == 1)
 {
    ret->x = head->P.x;
    ret->y = head->P.y;
    ret->z = head->P.z;
    *idx = head->idx;
   }
/*printf("intersections: %d\n", sol);*/
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", i,ret->x, ret->y, ret->z);*/
 if (sol <=1) { vlist_free(&head); return sol; }
/* printf("multiple soluions: searching nearest.\n");*/
 dist = 1e12;
 amb = NULL;
 hd = head;
 while (head)
   {
    if ((ndist = pseudo_distance(&r->P, &head->P)) < dist)
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
 /* FIXME: start form here, 3lines tottaly ignore ambiguous slutions for speed. */
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
 /* comment 3 lines above to have ambiguous solution checking */
 ts = 0;
 while (head)
   {
    if (nearly_equal(distance(&r->P, &head->P), dist, 1e-10))	/* FIXME: was 1e-9, the smallest ossible */
      {								/* to avoid ambiguous solutions */
       ts++;
       vlist_add(&amb, &head->P, head->idx);
       /* FIXME: break here to avoid ambiguous solutions */
/*       break;*/
      }
    head = head->next;
   }
 if (ts > 1)
 {
/*  tamb = amb;*/
/*  printf("Ambiguous solutions.\n");*/
     printf("*");
     fflush(stdout);
  /*printf("\nAbmigous solution for ray:\n");
  print_ray(r);
  while (amb)
    {
     printf("With triangle: %d, solution: (%f,%f,%f,), distance: %f\n",
	     amb->idx, amb->P.x, amb->P.y, amb->P.z, distance(&r->P, &amb->P));
     amb = amb->next;
    }
  amb = tamb;*/
  head = hd;
  vlist_free(&amb);
  vlist_free(&head);
  /* FIXME: must be handled somehow */
/*  panic(HERE);*/
  return sol;
 }
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", *idx,ret->x, ret->y, ret->z);*/
/* printf("big_free\n");*/
 head = hd;
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
}


int intersection_old(Triangle* tr, Ray* r, Vertex* ret, int* idx)
{
 register int i;
 int sol;
 double dist,ndist;
 int ts;
 VList* head;
 VList* hd;
 VList* amb;
/* VList* tamb;*/
 head = NULL;
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
 proc_tr += nTriangles;
 proc_tr2 += nTriangles;
/*printf("intersections: %d\n", sol);*/
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", i,ret->x, ret->y, ret->z);*/
 if (sol <=1) { vlist_free(&head); return sol; }
/* printf("multiple soluions: searching nearest.\n");*/
 dist = 1e12;
 amb = NULL;
 hd = head;
 while (head)
   {
    if ((ndist = pseudo_distance(&r->P, &head->P)) < dist)
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
 /* FIXME: start form here, 3lines tottaly ignore ambiguous slutions for speed. */
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
 /* comment 3 lines above to have ambiguous solution checking */
 ts = 0;
 while (head)
   {
    if (nearly_equal(pseudo_distance(&r->P, &head->P), dist, 1e-8))	/* FIXME: was 1e-9, the smallest ossible */
      {								/* to avoid ambiguous solutions */
       ts++;
       vlist_add(&amb, &head->P, head->idx);
       /* FIXME: break here to avoid ambiguous solutions */
/*       break;*/
      }
    head = head->next;
   }
 if (ts > 1)
 {
/*  tamb = amb;*/
/*  printf("Ambiguous solutions.\n");*/
     printf("*");
     fflush(stdout);
  /*printf("\nAbmigous solution for ray:\n");
  print_ray(r);
  while (amb)
    {
     printf("With triangle: %d, solution: (%f,%f,%f,), distance: %f\n",
	     amb->idx, amb->P.x, amb->P.y, amb->P.z, distance(&r->P, &amb->P));
     amb = amb->next;
    }
  amb = tamb;*/
  head = hd;
  vlist_free(&amb);
  vlist_free(&head);
  /* FIXME: must be handled somehow */
/*  panic(HERE);*/
  return sol;
 }
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", *idx,ret->x, ret->y, ret->z);*/
/* printf("big_free\n");*/
 head = hd;
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
}


int line_intersect(Vertex* p, Ray* l1, Ray* l2)
{
 double x1,y1,z1,x2,y2,z2;
 double dx1,dy1,dz1,dx2,dy2,dz2;
 register double K,L;
 double xr,yr,zr;
/* printf("LINE_INTERSECT:\n");*/
/* print_ray(l1);*/
/* print_ray(l2);*/
/* printf("LINE_INTERSECT_ENDS\n");*/
 /* FIXME: skip some check for speed */
/* if (!p || !l1 || !l2) panic(HERE);*/
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
 /* for speed suppose lines are not degraded */
 /*if (nearly_equal(dx1, 0., 1e-9) &&
	 nearly_equal(dy1, 0., 1e-9) &&
	 nearly_equal(dz1, 0., 1e-9)) return -1;
 if (nearly_equal(dx2, 0., 1e-9) &&
	 nearly_equal(dy2, 0., 1e-9) &&
	 nearly_equal(dz2, 0., 1e-9)) return -1; */
 if ( (nearly_equal(dx1, 0., 1e-7) &&		/* FIXME: was 1-e9 */
	 nearly_equal(dy1, 0., 1e-7)) ||
	 (nearly_equal(dx2, 0., 1e-7) &&
	 nearly_equal(dy2, 0., 1e-7)))
    {
     /* jedna z linii to  z = const */
     /* inny rzut*/
     if ( (nearly_equal(dx1, 0., 1e-7) &&
	 nearly_equal(dz1, 0., 1e-7)) ||
	 (nearly_equal(dx2, 0., 1e-7) &&
	 nearly_equal(dz2, 0., 1e-7)))
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
         if (!nearly_equal(dy1, 0., 1e-7)) K = (y2-y1+L*dy2)/dy1;
         else if (!nearly_equal(dz1, 0., 1e-9)) K = (z2-z1+L*dz2)/dz1;
	 else K = 0.;
         if (!nearly_equal(x1+K*dx1, x2+L*dx2, 3e-4))
	 {
	     fprintf(dbg,"not3dx\n");
	     return 0; /* not 3d intersect */
	 }
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
     if (!nearly_equal(dx1, 0., 1e-7)) K = (x2-x1+L*dx2)/dx1;
     else if (!nearly_equal(dz1, 0., 1e-9)) K = (z2-z1+L*dz2)/dz1;
     else K = 0.;
     if (!nearly_equal(y1+K*dy1, y2+L*dy2, 3e-4))
	 {
	     fprintf(dbg,"not3dy\n");
	     return 0; /* not 3d intersect */
	 }
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
 if (!nearly_equal(dx1, 0., 1e-7)) K = (x2-x1+L*dx2)/dx1;
 else if (!nearly_equal(dy1, 0., 1e-9)) K = (y2-y1+L*dy2)/dy1;
 else K = 0.;
 if (!nearly_equal(z1+K*dz1, z2+L*dz2, 3e-4))
	 {
/*	     printf("difference: %0.13lf\n", (z1+K*dz1) - (z2+L*dz2));
	     printf("x can be: %f or %f diff: %f\n", x1+K*dx1, x2+L*dx2, (x1+K*dx1) - (x2+L*dx2));
	     printf("y can be: %f or %f diff: %f\n", y1+K*dy1, y2+L*dy2, (y1+K*dy1) - (y2+L*dy2));
	     printf("z can be: %f or %f diff: %f\n", z1+K*dz1, z2+L*dz2, (z1+K*dz1) - (z2+L*dz2));*/
	     fprintf(dbg,"not3dz\n");
	     return 0; /* not 3d intersect */
	 }
 p->x = x1 + K*dx1;
 p->y = y1 + K*dy1;
 p->z = z1 + K*dz1;
 return 1;	/* intersection OK */
 p_case:
 /* FIXME: seldom met */
/* printf("p_case!\n");*/
 xr = (x2 - x1) / dx1;
 yr = (y2 - y1) / dy1;
 zr = (z2 - z1) / dz1;
 if (!nearly_equal(xr, yr, 1e-8)) return 0; /* parallel disjoint */
 if (!nearly_equal(xr, zr, 1e-8)) return 0; /* parallel disjoint */
 if (!nearly_equal(yr, zr, 1e-8)) return 0; /* parallel disjoint */
 p->x = x1;
 p->y = y1;
 p->z = z1;
 return 2;	/* parallel the same*/
}


void make_ray(Ray* l, Vertex* a, Vertex* b)
{
/* if (!a || !b || !l) panic(HERE);*/
 l->P.x = a->x;
 l->P.y = a->y;
 l->P.z = a->z;
 l->d.x = b->x - a->x;
 l->d.y = b->y - a->y;
 l->d.z = b->z - a->z;
 l->r = 0;
 l->c = RED;
}


double pseudo_distance(Vertex* a, Vertex* b)
{
 double x,y,z;
 x = b->x - a->x;
 y = b->y - a->y;
 z = b->z - a->z;
 return x*x+y*y+z*z;
}


double distance(Vertex* a, Vertex* b)
{
 double x,y,z;
 x = b->x - a->x;
 y = b->y - a->y;
 z = b->z - a->z;
 return sqrt(x*x+y*y+z*z);
}


void normalize(Vector* v)
{
 double len;
 len = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
 /* FIXME: for speed: dangerous */
/* if (nearly_equal(len, 0., 1e-9)) panic(HERE);*/
 v->x /= len;
 v->y /= len;
 v->z /= len;
}


void materialize(Material* t)
{
 double f;
 f = t->c + t->s + t->t;
 if (f <= 1.) return;
 t->c /= f;
 t->s /= f;
 t->t /= f;
}


void get_normal_new(Vector* n, Triangle* t, Vertex* p, Material* R, Material* G, Material* B, TexCoord* tex)
{
 Ray ab,bc,ca;
 Ray ap, bp, cp;
 Vertex pa,pb,pc;
 TexCoord ttex;
 register double la,lb,lc;
 double sum,fact;
 int err,sol;
 make_ray(&ab, &t->a, &t->b);
 make_ray(&bc, &t->b, &t->c);
 make_ray(&ca, &t->c, &t->a);
 make_ray(&ap, &t->a, p);
 make_ray(&bp, &t->b, p);
 make_ray(&cp, &t->c, p);
 err = 0;
 la = lb = lc = 0.;
 /* FIXME: ignore unsolved triangles, just use values from a unsolved vetex */
 if ((sol = line_intersect(&pa, &bc, &ap)) != 1) { la = 1.; lb = lc = 0.; err = 1; /*printf("+");*/ goto solved; }
 if ((sol = line_intersect(&pb, &ca, &bp)) != 1) { lb = 1.; la = lc = 0.; err = 1; /*printf("+");*/ goto solved; }
 if ((sol = line_intersect(&pc, &ab, &cp)) != 1) { lc = 1.; la = lb = 0.; err = 1; /*printf("+");*/ goto solved; }
 solved:
 if (!err)
   {
    la = distance(&pa, p);
    lb = distance(&pb, p);
    lc = distance(&pc, p);
    sum = la+lb+lc;
    la /= sum;
    lb /= sum;
    lc /= sum;
   }
 else
   {
    if (la > .5) { tex->x = t->ta.x; tex->y = t->ta.y; }
    if (lb > .5) { tex->x = t->tb.x; tex->y = t->tb.y; }
    if (lc > .5) { tex->x = t->tc.x; tex->y = t->tc.y; }
   }
 n->x = la*t->na.x + lb*t->nb.x + lc*t->nc.x;
 n->y = la*t->na.y + lb*t->nb.y + lc*t->nc.y;
 n->z = la*t->na.z + lb*t->nb.z + lc*t->nc.z;
 R->c = la*t->mra.c + lb*t->mrb.c + lc*t->mrc.c;
 R->s = la*t->mra.s + lb*t->mrb.s + lc*t->mrc.s;
 R->t = la*t->mra.t + lb*t->mrb.t + lc*t->mrc.t;
 G->c = la*t->mga.c + lb*t->mgb.c + lc*t->mgc.c;
 G->s = la*t->mga.s + lb*t->mgb.s + lc*t->mgc.s;
 G->t = la*t->mga.t + lb*t->mgb.t + lc*t->mgc.t;
 B->c = la*t->mba.c + lb*t->mbb.c + lc*t->mbc.c;
 B->s = la*t->mba.s + lb*t->mbb.s + lc*t->mbc.s;
 B->t = la*t->mba.t + lb*t->mbb.t + lc*t->mbc.t;
 /*tex->x = la*t->ta.x + lb*t->tb.x + lc*t->tc.x;
 tex->y = la*t->ta.y + lb*t->tb.y + lc*t->tc.y;*/
 materialize(R);
 materialize(G);
 materialize(B);
 normalize(n);
 /***********************/
 if (err) return;
 if (nearly_equal(t->a.x, pa.x, 1e-7) &&
	 nearly_equal(t->a.y, pa.y, 1e-7) &&
	 nearly_equal(t->a.y, pa.y, 1e-7))
   {
    tex->x = t->ta.x;
    tex->y = t->ta.y;
    return;
   }
 la = pseudo_distance(&t->c, &t->b);
 lb = pseudo_distance(&pa, &t->b);
 fact = lb / la;
 fact = sqrt(fact);
 ttex.x = (1.-fact)* t->tb.x + fact*t->tc.x;
 ttex.y = (1.-fact)* t->tb.y + fact*t->tc.y;
 la = pseudo_distance(&pa, &t->a);
 lb = pseudo_distance(p, &t->a);
 fact = lb / la;
 fact = sqrt(fact);
 tex->x = (1.-fact)* t->ta.x + fact*ttex.x;
 tex->y = (1.-fact)* t->ta.y + fact*ttex.y;
 /***********************/
}


void get_normal_old(Vector* n, Triangle* t, Vertex* p, Material* R, Material* G, Material* B, TexCoord* tex)
{
 Ray bc;
 Ray ap;
 int sol;
 Vertex d;
 Vector tn;
 Material tR, tG, tB;
 TexCoord ttex;
 double lenBC, lenBD;
 double fact;
 make_ray(&bc, &t->b, &t->c);
 make_ray(&ap, &t->a, p);
 /*printf("P = (%f,%f,%f)\n", p->x, p->y, p->z);
 print_ray(&bc);
 print_ray(&ap);*/;
 if ((sol = line_intersect(&d, &bc, &ap)) != 1)
   {
/*    printf("Line intersect returns: %d\n", sol);*/
/*       printf("+");*/
       fflush(stdout);
    goto get_a;
/*    panic(HERE);*/
    /*FIXME: this is a stub */
   }
 if (nearly_equal(t->a.x, d.x, 1e-7) &&
	 nearly_equal(t->a.y, d.y, 1e-7) &&
	 nearly_equal(t->a.y, d.y, 1e-7))
   {
    get_a:
/*       printf("-");*/
    n->x = t->na.x;
    n->y = t->na.y;
    n->z = t->na.z;
    R->c = t->mra.c;
    R->s = t->mra.s;
    R->t = t->mra.t;
    G->c = t->mga.c;
    G->s = t->mga.s;
    G->t = t->mga.t;
    B->c = t->mba.c;
    B->s = t->mba.s;
    B->t = t->mba.t;
    tex->x = t->ta.x;
    tex->y = t->ta.y;
    return ;
   }
 lenBC = pseudo_distance(&t->c, &t->b);
/* if (nearly_equal(lenBC, 0., 1e-9)) panic(HERE);*/
 lenBD = pseudo_distance(&d, &t->b);
 fact = lenBD / lenBC;
 fact = sqrt(fact);
 tn.x = (1.-fact)* t->nb.x + fact*t->nc.x;
 tn.y = (1.-fact)* t->nb.y + fact*t->nc.y;
 tn.z = (1.-fact)* t->nb.z + fact*t->nc.z;
 tR.c = (1.-fact)* t->mrb.c + fact*t->mrc.c;
 tR.s = (1.-fact)* t->mrb.s + fact*t->mrc.s;
 tR.t = (1.-fact)* t->mrb.t + fact*t->mrc.t;
 tG.c = (1.-fact)* t->mgb.c + fact*t->mgc.c;
 tG.s = (1.-fact)* t->mgb.s + fact*t->mgc.s;
 tG.t = (1.-fact)* t->mgb.t + fact*t->mgc.t;
 tB.c = (1.-fact)* t->mbb.c + fact*t->mbc.c;
 tB.s = (1.-fact)* t->mbb.s + fact*t->mbc.s;
 tB.t = (1.-fact)* t->mbb.t + fact*t->mbc.t;
 ttex.x = (1.-fact)* t->tb.x + fact*t->tc.x;
 ttex.y = (1.-fact)* t->tb.y + fact*t->tc.y;
 normalize(&tn);
 /* FIXME: materialize is needed here? */
 lenBC = pseudo_distance(&d, &t->a);
/* if (nearly_equal(lenBC, 0., 1e-9)) panic(HERE);*/
 lenBD = pseudo_distance(p, &t->a);
 fact = lenBD / lenBC;
 fact = sqrt(fact);
 n->x = (1.-fact)* t->na.x + fact*tn.x;
 n->y = (1.-fact)* t->na.y + fact*tn.y;
 n->z = (1.-fact)* t->na.z + fact*tn.z;
 R->c = (1.-fact)* t->mra.c + fact*tR.c;
 R->s = (1.-fact)* t->mra.s + fact*tR.s;
 R->t = (1.-fact)* t->mra.t + fact*tR.t;
 G->c = (1.-fact)* t->mga.c + fact*tG.c;
 G->s = (1.-fact)* t->mga.s + fact*tG.s;
 G->t = (1.-fact)* t->mga.t + fact*tG.t;
 B->c = (1.-fact)* t->mba.c + fact*tB.c;
 B->s = (1.-fact)* t->mba.s + fact*tB.s;
 B->t = (1.-fact)* t->mba.t + fact*tB.t;
 tex->x = (1.-fact)* t->ta.x + fact*ttex.x;
 tex->y = (1.-fact)* t->ta.y + fact*ttex.y;
 materialize(R);
 materialize(G);
 materialize(B);
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


double pseudo_length(Vector* v)
{
 return v->x*v->x+v->y*v->y+v->z*v->z;
}

int get_red(Screen* s, int x, int y);
int get_green(Screen* s, int x, int y);
int get_blue(Screen* s, int x, int y);

void get_texture(Texture* t, double x, double y, double *r, double* g, double* b)
{
 int i,j;
 i = (int)(x*(double)t->x);
 j = (int)(y*(double)t->y);
/* printf("(%f,%f) from [%d,%d] --> (%d,%d)\n", x,y,t->x,t->y,i,j);*/
 *r = (double)get_red(t, i, j)/255.;
 *g = (double)get_green(t, i, j)/255.;
 *b = (double)get_blue(t, i, j)/255.;
}


double recurse_color(Ray* r, Triangle* ptrlist, Vertex* v, Vector* n, Material* R, Material* G, Material* B, TexCoord* tC, int id, double f)
{
 double mU,mD;
 double rv;
 double ca;
 double fact, fact2;
 double a2;
 double nf;
 double add;
 double rt,gt,bt;
 double beta,arg;
 double l1,l2,shadow;
 int idx;
 Triangle* t;
 Vector obsV, nn, tr, pn, pv, in;
 Vertex nv;
 Material nR, nG, nB;
 TexCoord ntC;
 Ray vi;
 Ray ti;
 Ray li;
 fact = fact2 = 0.;
 mU = mD = 1.;
 t = &ptrlist[id];
/* printf("%f,%f)\n", tC->x, tC->y);*/
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
 if (pseudo_length(&obsV) < 3e-4) return 0.;
 if (pseudo_length(&tr) < 3e-4) return 0.;
 li.P.x = v->x;
 li.P.y = v->y;
 li.P.z = v->z;
 li.d.x = tr.x;
 li.d.y = tr.y;
 li.d.z = tr.z;
 li.r = r->r + 1;
 li.c = r->c;
 normalize(&tr);
 normalize(&obsV);
 shadow = 1.;
 if (intersection(ptrlist, &li, &nv, &idx))
   {
    l1 = pseudo_distance(v, &light);
    l2 = pseudo_distance(v, &nv);
    if (l2 < l1)
      {
       get_normal(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
       if (r->c == RED)   shadow = nR.t;
       if (r->c == GREEN) shadow = nG.t;
       if (r->c == BLUE)  shadow = nB.t;
       if (shadow > .9) shadow = .9;
      }
   }
 if (shadow < .25) shadow = .25;
/* if (shadow > 1.) shadow = 1.;*/
/* printf("nL: (%f,%f,%f)\n", tr.x, tr.y, tr.z);*/
/* printf("nV: (%f,%f,%f)\n", obsV.x, obsV.y, obsV.z);*/
 ca = scalar_prod(&tr, n);
/* printf("ca = %f\n", ca);*/
 if (ca < 0. && t->faces == 1) ca = 0.;
 if (ca < 0. && t->faces == 2) ca *= -1.;
/* printf("nca = %f\n", ca);*/
 if (r->c == RED)   fact = R->c;
 if (r->c == GREEN) fact = G->c;
 if (r->c == BLUE)  fact = B->c;
 if (t->tex)
 {
  get_texture(t->tex, tC->x, tC->y, &rt, &gt, &bt);
  if (r->c == RED)   fact *= rt;
  if (r->c == GREEN) fact *= gt;
  if (r->c == BLUE)  fact *= bt;
 }
/* else printf(" no texture.\n");*/
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
/* if (ca < 0.) ca = 0.;*/
 if (ca > 1.) ca = 1.;
/* printf("ca = %f\n", ca);*/
/* printf("diffuse adds: %f\n", fact * ca);*/
 rv += fact * ca * shadow;
 /* FIXME: is this specular flash OK? */
 if (t->ca > 0.)
   {
    ca = scalar_prod(&vi.d, &tr);
    if (ca > 0.) rv += shadow * pow(ca, t->ca);
   }
 if (rv > 1.) rv = 1.;
/* printf("rv = %f\n", rv);*/
 if (r->c == RED)   fact = R->s /*+ t->tr*/;
 if (r->c == GREEN) fact = G->s /*+ t->tg*/;
 if (r->c == BLUE)  fact = B->s /*+ t->tb*/;
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
 if (nf > step && (rv < 1.-step))
   {
/*       return rv;*/
    /* FIXME: continue here! */
    /*oblicz odbicie...*/
    if (!intersection(ptrlist, &vi, &nv, &idx))
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
         get_normal(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
	 add = nf * recurse_color(&vi, ptrlist, &nv, &nn, &nR, &nG, &nB, &ntC, idx, nf);
	 rv += add;
       }
   }
 if (r->c == RED)   { fact = R->t; mU = t->mUR; mD = t->mDR; }
 if (r->c == GREEN) { fact = G->t; mU = t->mUG; mD = t->mDG; }
 if (r->c == BLUE)  { fact = B->t; mU = t->mUB; mD = t->mDB; }
 nf = f * fact;
 if (nf > step && (rv <1.-step))
   {
    ti.P.x = v->x;
    ti.P.y = v->y;
    ti.P.z = v->z;
    a2 = scalar_prod(&obsV, n);
    if (nearly_equal(a2, 0., 1e-8)) goto ret_val;
    if (a2 > 0.)
      {
    /*zalamanie z osrodka mU do mD*/		/* FIXME: dalej nie jest pewne czy zalamanie dziala dobrze. */
       arg = ((sin(acos(a2)))*mU)/mD;
       if (arg < .999999) beta = asin(arg);
       else beta = 1.5705;
/*    printf("alfa = %f, beta = %f\n", (acos(a2)*180.)/(3.1415926), (beta*180.)/(3.1415929));*/
       arg = cos(beta);
      pn.x = -(n->x * arg);
      pn.y = -(n->y * arg);
      pn.z = -(n->z * arg);
      arg = mU/mD;
      pv.x = (a2*n->x - obsV.x)*arg;
      pv.y = (a2*n->y - obsV.y)*arg;
      pv.z = (a2*n->z - obsV.z)*arg;
    /*printf("spr = %f\n", scalar_prod(&pn, &pv));
    printf("n = (%f,%f,%f)\n", n->x, n->y, n->z);
    printf("pn = (%f,%f,%f)\n", pn.x, pn.y, pn.z);
    printf("lpn = (%f,%f,%f)\n", n->x*a2, n->y*a2, n->z*a2);
    printf("v = (%f,%f,%f)\n", obsV.x-a2*n->x, obsV.y-a2*n->y,obsV.z-a2*n->z);
    printf("pv = (%f,%f,%f)\n", pv.x, pv.y, pv.z);*/
      ti.d.x = pn.x + pv.x;
      ti.d.y = pn.y + pv.y;
      ti.d.z = pn.z + pv.z;
      normalize(&ti.d);
     }
   if (a2 < 0.)
     {
/*       printf("from denser!\n");*/
      in.x = -n->x;
      in.y = -n->y;
      in.z = -n->z;
      a2 = scalar_prod(&obsV, &in);
/*    printf("a2 = %f\n", a2);*/
    /*zalamanie z osrodka mD do mU*/
      arg = ((sin(acos(a2)))*mD)/mU;
/*    printf("arg = %f\n", arg);*/
      if (arg < .999999) beta = asin(arg);
      else beta = 1.5705;
/*    if (beta < 3.14 / 4.) printf("alfa = %f, beta = %f\n", (acos(a2)*180.)/(3.1415926), (beta*180.)/(3.1415929));*/
      arg = cos(beta);
      pn.x = n->x * cos(beta);
      pn.y = n->y * cos(beta);
      pn.z = n->z * cos(beta);
      arg = mD/mU;
      pv.x = (a2*in.x - obsV.x)*arg;
      pv.y = (a2*in.y - obsV.y)*arg;
      pv.z = (a2*in.z - obsV.z)*arg;
    		/* FIXME: is pv correctly computed ? */
      /*printf("spr = %f\n", scalar_prod(&pn, &pv));
      printf("n = (%f,%f,%f)\n", n->x, n->y, n->z);
      printf("pn = (%f,%f,%f)\n", pn.x, pn.y, pn.z);
      printf("lpn = (%f,%f,%f)\n", n->x*a2, n->y*a2, n->z*a2);
      printf("v = (%f,%f,%f)\n", obsV.x-a2*n->x, obsV.y-a2*n->y,obsV.z-a2*n->z);
      printf("pv = (%f,%f,%f)\n", pv.x, pv.y, pv.z);*/
      ti.d.x = pn.x + pv.x;
      ti.d.y = pn.y + pv.y;
      ti.d.z = pn.z + pv.z;
      normalize(&ti.d);
     }
   /*printf("t:(%1.2f/%1.2f)  (%2.2f) (%3.2f -> %3.2f) (%1.2f,%1.2f,%1.2f) -> (%1.2f,%1.2f,%1.2f)\n",
	t->mU, t->mD, 180.*acos(a2)/3.1415 - 180.*beta/3.1415, 180.*acos(a2)/3.1415, 180.*beta/3.1415, -obsV.x, -obsV.y, -obsV.z, ti.d.x, ti.d.y, ti.d.z);*/
    ti.r = r->r + 1;
    ti.c = r->c;
/*       return rv;*/
    /* FIXME: continue here! */
    /*oblicz odbicie...*/
    if (!intersection(ptrlist, &ti, &nv, &idx))
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
         get_normal(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
	 add = nf * recurse_color(&ti, ptrlist, &nv, &nn, &nR, &nG, &nB, &ntC, idx, nf);
	 rv += add;
       }
   }
 ret_val:
/* else printf("skipped!\n");*/
 if (rv <= step) return 0.;
 if (rv > 1.) return 1.;
/* printf("FINAL rv = %f\n", rv);*/
 return rv;
}


void calculate_color(Screen* s, Ray* r, Triangle* t, int x, int y)
{
 Vertex v;
 Vector n;
 Material Rm, Gm, Bm;
 TexCoord tC;
 double rr,gg,bb;
 int idx;
 /* moze byc wiele przeciec z obiektami, posortowac w/g kolejnosci */
 if (!intersection(t, r, &v, &idx))
   {
/*     printf("-");*/
/*       printf("error with:\n");*/
/*       print_ray(r);*/
/*       fprintf(dbg,"no intr at: (%d,%d)\n", x,y);*/
/*       fprintf(dbg, "(%f,%f,%f) --> (%f,%f,%f)\n", r->P.x, r->P.y, r->P.z, r->d.x, r->d.y, r->d.z);*/
       color_background(s, x, y);
/*       set_color(s, x, y, 0,0,0);*/
   }
 else
   {
/*       printf("intersect: (%d,%d)\n", x,y);*/
/*       printf("V = (%f,%f,%f)\n", v.x, v.y, v.z);*/
       get_normal(&n, &t[idx], &v, &Rm, &Gm, &Bm, &tC);
       r->c = RED;
       rr = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1.);
       r->c = GREEN;
       gg = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1.);
       r->c = BLUE;
       bb = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1.);
/*       if (rr > .99)*/
/*       fprintf(dbg,"COMPUTED: RGB: (%f,%f,%f) at (%d,%d)\n", rr,gg,bb,x,y);*/
       set_color(s, x, y, (int)(rr*255.), (int)(gg*255.), (int)(bb*255.));
/*       panic(HERE);*/		/* BIG_PANIC */
/*       printf("*");*/
   }
}


double calc_surface(Triangle* t)
{
 double alfa;
 double la,lb;
 Vector a,b;
 a.x = t->c.x - t->a.x;
 a.y = t->c.y - t->a.y;
 a.z = t->c.z - t->a.z;
 b.x = t->b.x - t->a.x;
 b.y = t->b.y - t->a.y;
 b.z = t->b.z - t->a.z;
 la = length(&a);
 lb = length(&b);
 if (la < 1e-7 || lb < 1e-7)
   {
    printf("Warning: 0 length edge\n");
    return 0.;
   }
 normalize(&a);
 normalize(&b);
 alfa = acos(scalar_prod(&a, &b));
 return .5 * sin(alfa) * la * lb;
}


PtrList* get_head(PtrList* ptr)
{
 while (ptr->prev) ptr = ptr->prev;
 return ptr;
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


double calc_box(Triangle* t1, Triangle* t2)
{
 double minx,miny,minz;
 double maxx,maxy,maxz;
 double dx,dy,dz;
/* if (!t1 || !t2) panic(HERE);*/
 minx = t1->a.x;
 if (t1->b.x < minx) minx = t1->b.x;
 if (t1->c.x < minx) minx = t1->c.x;
 if (t2->a.x < minx) minx = t2->a.x;
 if (t2->b.x < minx) minx = t2->b.x;
 if (t2->c.x < minx) minx = t2->c.x;
 miny = t1->a.y;
 if (t1->b.y < miny) miny = t1->b.y;
 if (t1->c.y < miny) miny = t1->c.y;
 if (t2->a.y < miny) miny = t2->a.y;
 if (t2->b.y < miny) miny = t2->b.y;
 if (t2->c.y < miny) miny = t2->c.y;
 minz = t1->a.z;
 if (t1->b.z < minz) minz = t1->b.z;
 if (t1->c.z < minz) minz = t1->c.z;
 if (t2->a.z < minz) minz = t2->a.z;
 if (t2->b.z < minz) minz = t2->b.z;
 if (t2->c.z < minz) minz = t2->c.z;
 maxx = t1->a.x;
 if (t1->b.x > maxx) maxx = t1->b.x;
 if (t1->c.x > maxx) maxx = t1->c.x;
 if (t2->a.x > maxx) maxx = t2->a.x;
 if (t2->b.x > maxx) maxx = t2->b.x;
 if (t2->c.x > maxx) maxx = t2->c.x;
 maxy = t1->a.y;
 if (t1->b.y > maxy) maxy = t1->b.y;
 if (t1->c.y > maxy) maxy = t1->c.y;
 if (t2->a.y > maxy) maxy = t2->a.y;
 if (t2->b.y > maxy) maxy = t2->b.y;
 if (t2->c.y > maxy) maxy = t2->c.y;
 maxz = t1->a.z;
 if (t1->b.z > maxz) maxz = t1->b.z;
 if (t1->c.z > maxz) maxz = t1->c.z;
 if (t2->a.z > maxz) maxz = t2->a.z;
 if (t2->b.z > maxz) maxz = t2->b.z;
 if (t2->c.z > maxz) maxz = t2->c.z;
 dx = maxx - minx;
 dy = maxy - miny;
 dz = maxz - minz;
/* return sqrt(dx*dx + dy*dy + dz*dz);*/
 return dx*dx + dy*dy + dz*dz;
}


void add_box(BList** head, Triangle* t1, Triangle* t2)
{
 double minx,miny,minz;
 double maxx,maxy,maxz;
/* if (!head || ! t1) panic(HERE);*/
 minx = t1->a.x;
 if (t1->b.x < minx) minx = t1->b.x;
 if (t1->c.x < minx) minx = t1->c.x;
 if (t2)
   {
    if (t2->a.x < minx) minx = t2->a.x;
    if (t2->b.x < minx) minx = t2->b.x;
    if (t2->c.x < minx) minx = t2->c.x;
   }
 miny = t1->a.y;
 if (t1->b.y < miny) miny = t1->b.y;
 if (t1->c.y < miny) miny = t1->c.y;
 if (t2)
   {
    if (t2->a.y < miny) miny = t2->a.y;
    if (t2->b.y < miny) miny = t2->b.y;
    if (t2->c.y < miny) miny = t2->c.y;
   }
 minz = t1->a.z;
 if (t1->b.z < minz) minz = t1->b.z;
 if (t1->c.z < minz) minz = t1->c.z;
 if (t2)
   {
    if (t2->a.z < minz) minz = t2->a.z;
    if (t2->b.z < minz) minz = t2->b.z;
    if (t2->c.z < minz) minz = t2->c.z;
   }
 maxx = t1->a.x;
 if (t1->b.x > maxx) maxx = t1->b.x;
 if (t1->c.x > maxx) maxx = t1->c.x;
 if (t2)
   {
    if (t2->a.x > maxx) maxx = t2->a.x;
    if (t2->b.x > maxx) maxx = t2->b.x;
    if (t2->c.x > maxx) maxx = t2->c.x;
   }
 maxy = t1->a.y;
 if (t1->b.y > maxy) maxy = t1->b.y;
 if (t1->c.y > maxy) maxy = t1->c.y;
 if (t2)
   {
    if (t2->a.y > maxy) maxy = t2->a.y;
    if (t2->b.y > maxy) maxy = t2->b.y;
    if (t2->c.y > maxy) maxy = t2->c.y;
   }
 maxz = t1->a.z;
 if (t1->b.z > maxz) maxz = t1->b.z;
 if (t1->c.z > maxz) maxz = t1->c.z;
 if (t2)
   {
    if (t2->a.z > maxz) maxz = t2->a.z;
    if (t2->b.z > maxz) maxz = t2->b.z;
    if (t2->c.z > maxz) maxz = t2->c.z;
   }
 blist_add(head, minx, miny, minz, maxx, maxy, maxz, t1, t2);
/* printf("BOX: (%f,%f,%f) --> (%f,%f,%f)\n", minx,miny,minz,maxx,maxy,maxz);*/
}


double calc_boxes(BTree* b1, BTree* b2)
{
 double minx,miny,minz;
 double maxx,maxy,maxz;
 double dx,dy,dz;
/* if (!b1 || !b2) panic(HERE);*/
 dx = dy = dz = 0.;
 minx = b1->b->minx;
 if (b2->b->minx < minx) minx = b2->b->minx;
 miny = b1->b->miny;
 if (b2->b->miny < miny) miny = b2->b->miny;
 minz = b1->b->minz;
 if (b2->b->minz < minz) minz = b2->b->minz;
 maxx = b1->b->maxx;
 if (b2->b->maxx > maxx) maxx = b2->b->maxx;
 maxy = b1->b->maxy;
 if (b2->b->maxy > maxy) maxy = b2->b->maxy;
 maxz = b1->b->maxz;
 if (b2->b->maxz > maxz) maxz = b2->b->maxz;
/* return sqrt(dx*dx+dy*dy+dz*dz);*/
 return dx*dx+dy*dy+dz*dz;
}


void create_btree_box(Box** b, BTree* b1, BTree* b2)
{
 Box* temp;
 double minx,miny,minz;
 double maxx,maxy,maxz;
/* if (!b || !b1 || !b2) panic(HERE);*/
 minx = b1->b->minx;
 if (b2->b->minx < minx) minx = b2->b->minx;
 miny = b1->b->miny;
 if (b2->b->miny < miny) miny = b2->b->miny;
 minz = b1->b->minz;
 if (b2->b->minz < minz) minz = b2->b->minz;
 maxx = b1->b->maxx;
 if (b2->b->maxx > maxx) maxx = b2->b->maxx;
 maxy = b1->b->maxy;
 if (b2->b->maxy > maxy) maxy = b2->b->maxy;
 maxz = b1->b->maxz;
 if (b2->b->maxz > maxz) maxz = b2->b->maxz;
 temp = (Box*)malloc(sizeof(Box));
 temp->t1 = temp->t2 = NULL;	/* abstract box-> node */
 temp->minx = minx;
 temp->miny = miny;
 temp->minz = minz;
 temp->maxx = maxx;
 temp->maxy = maxy;
 temp->maxz = maxz;
 *b = temp;
}


BTree* find_nearest_boxes(PtrList* head, PtrList** b1, PtrList** b2, int swap)
{
 BTree* nbox;
 PtrList *p1, *p2;
 PtrList *bl, *br;
 double minD,d;
/* if (!head || !b1 || !b2) panic(HERE);*/
 p1 = head;
 minD = 1e12;
 *b1 = *b2 = NULL;
 bl = br = NULL;
 while (p1)
   {
    p2 = p1->next;
    while (p2)
      {
       if ( (d = calc_boxes((BTree*)p1->ptr, (BTree*)p2->ptr)) < minD)
         {
          bl = p1;
	  br = p2;
         }
       p2 = p2->next;
      }
    p1 = p1->next;
   }
/* if (!bl || !br) return NULL;*/
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


int count_items(PtrList* head)
{
 int item;
 item = 0;
 while (head) { item++; head = head->next; }
 return item;
}


void create_btree(BList* head)
{
 BList *hd;
 PtrList* phead;
 BTree* temp;
 PtrList *bt1, *bt2;
 int no;
 hd = head;
 phead = NULL;
 bt1 = bt2 = NULL;
/* if (!head) panic(HERE);*/
 no = 0;
 while (head)
   {
    temp = (BTree*)malloc(sizeof(BTree));
    temp->l = temp->r = NULL;	/* real box: leaf */
    temp->b = &head->b;		/* uses a real box containing triangle(s) */
    ptrlist_add(&phead, (void*)temp);
    head = head->next;
   }
 while ((no = count_items(phead)) > 1)
   {
/*       printf("creating btree... items: %d\n", no);*/
    if (!(no % 32)) { printf("."); fflush(stdout); }
    temp = find_nearest_boxes(phead, &bt1, &bt2, no % 2);
    ptrlist_delete(&phead, bt1);
    ptrlist_delete(&phead, bt2);
    ptrlist_add(&phead, (void*)temp);
   }
 btree = (BTree*)phead->ptr;
 ptrlist_free(&phead);
}


void preprocess_scene(Triangle* tl)
{
 int i;
 double minS,s;
 PtrList *head, *hd;
 PtrList *minT;
 BList * bhead;
 Triangle *t1, *t2;
 printf("creating btree...");
 fflush(stdout);
 head = NULL;
 bhead = NULL;
 for (i=0;i<nTriangles;i++) ptrlist_add(&head, &tl[i]);
 while (head)
 {
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
 /* printf("minT: (%f,%f,%f), (%f,%f,%f), (%f,%f,%f) %f\n",
	 minT->t->a.x, minT->t->a.y, minT->t->a.z,
	 minT->t->b.x, minT->t->b.y, minT->t->b.z,
	 minT->t->c.x, minT->t->c.y, minT->t->c.z, minS);*/
/*  printf("minS: %f ", minS);*/
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
     /*printf("boxD: (%f,%f,%f), (%f,%f,%f), (%f,%f,%f) %f\n",
	 minT->t->a.x, minT->t->a.y, minT->t->a.z,
	 minT->t->b.x, minT->t->b.y, minT->t->b.z,
	 minT->t->c.x, minT->t->c.y, minT->t->c.z, minS);*/
     t2 = (Triangle*)minT->ptr;
/*     printf("boxD: %f\n", minS);*/
     head = hd;
     ptrlist_delete(&head, minT);
     hd = head;
    }
/*  else printf("just one triangle\n");*/
  add_box(&bhead, t1, t2);
 }
 ptrlist_free(&head);
 boxes = bhead;
 create_btree(bhead);
 printf("done.\n");
 /*while (btree)
 {
 printf("(%p [%p (%3.0f,%3.0f,%3.0f) (%3.0f,%3.0f,%3.0f) %p] %p)\n",
	 btree->l,btree->b->t1,
	 btree->b->minx, btree->b->miny, btree->b->minz,
	 btree->b->maxx, btree->b->maxy, btree->b->maxz,
	 btree->b->t2,btree->r);
 btree = btree->l;
 }*/
}


void raytrace(Screen* s, char* scenefile, int bkup)
{
 Ray r;
 Triangle* ts;
 int i,j;
 dbg = fopen("debug.out", "w");
 if (!dbg) panic(HERE);
 step = 2./256.;
/* nmalloc = nfree = 0;*/
 load_scene(&ts, scenefile);
 preprocess_scene(ts);
 get_normal = get_normal_old;
 if (nTriangles <= 8) intersection = intersection_old;
 else intersection = intersection_new;
 r.P.x = observer.x;
 r.P.y = observer.y;
 r.P.z = observer.z;;
 r.r = 0;
 r.d.z = (s->x+s->y)/3.;	/* FIXME: mozna zmieniac aby dokladniej zobaczyc */
 printf("\n");
 /*calculate_color(s, &r, ts, 400, 300);*/
 /*printf("TEST1\n");
 r.d.x = 255 - s->x/2.;
 r.d.y = 191 - s->y/2.;
 calculate_color(s, &r, ts, 255, 191);
 panic(HERE);*/
 /*printf("TEST\n");
 r.d.x = 256 - s->x/2.;
 r.d.y = 256 - s->y/2.;
 calculate_color(s, &r, ts, 256, 256);
 panic(HERE);*/
 for (i=0;i<s->x;i++)
   {
    r.d.x = i - s->x/2.;
    printf("%04d/%04d", i+1,s->x);
    proc_tr = proc_tr2 = 0;
    for (j=0;j<s->y;j++)
      {
       if (!(j % 24)) { printf("."); fflush(stdout); }
       r.d.y = j - s->y/2.;
       calculate_color(s, &r, ts, i, j);
      }
    printf("%3.1f%%,%3.1f%%,%dK/%dK\n", ((double)(i+1)*100.)/(double)s->x,(double)(proc_tr*100)/(double)proc_tr2,proc_tr/1000,proc_tr2/1000);
    if (i > 0 && !(i % bkup))
      {
       wrt_bmp(&screen, "partial.bmp");
      }
/*    printf("\n");*/
   }
 printf("\n");
 free_scene(&ts);
 fclose(dbg);
}


int get_red(Screen* s, int x, int y)
{
    /* FIXME: fast version */
    return s->pixels[3*(s->y * x + y)];
 if (x == s->x) x--;
 if (y == s->y) y--;
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)];
   }
 else panic(HERE);
 return -1;
}


int get_green(Screen* s, int x, int y)
{
    /* FIXME: fast version */
    return s->pixels[3*(s->y * x + y)+1];
 if (x == s->x) x--;
 if (y == s->y) y--;
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)+1];
   }
 else panic(HERE);
 return -1;
}


int get_blue(Screen* s, int x, int y)
{
    /* FIXME: fast version */
    return s->pixels[3*(s->y * x + y)+2];
 if (x == s->x) x--;
 if (y == s->y) y--;
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
 if (!plik) { printf("Error writing BMP: %s\n", out_f); exit(0); }
 fprintf(plik,"%c%c",'B', 'M');
 bm_handle.bm_y = s->y;
 bm_handle.bm_x = s->x;
 bm_handle.fsize = sizeof(BMPTag)+(bm_handle.bm_y*bm_handle.bm_x*3);
 fwrite(&bm_handle.fsize,4,1,plik);
 fwrite(&bm_handle.dummy,4,1,plik);
 bm_handle.offset=sizeof(BMPTag);
 bm_handle.planes=1;
 bm_handle.bpp=24;
/* bm_handle.nbytes = s->x * s->y * 3;*/
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
 printf("Bitmap: %s written.\n", out_f);
}


int main(int lb, char** par)
{
 int bkup;
 bkup = 64;
 if (lb < 2) { printf("usage: %s scenefile [nlines_per_backup]\n",par[0]); return 1; }
 if (lb > 2) bkup = atoi(par[2]);
 raytrace(&screen, par[1], bkup);
 wrt_bmp(&screen, "screen.bmp");
 free_screen(&screen);
 blist_free(&boxes);
 return 0;
}

