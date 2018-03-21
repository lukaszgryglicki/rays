#ifndef __RT_LIB_H__
#define __RT_LIB_H__
#define HERE __FILE__,__LINE__
#define OPEN_TEXTURE 0
#define OPEN_PARTIAL 1
#define PI 3.14159265

#define __REAL_80bit  long double
#define __REAL_64bit  double
#define __REAL_32bit  float
#define sREAL double

#ifndef REAL
#define REAL __REAL_80bit
#endif

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
#define Smin 12

#include <stdio.h>
#include <stdlib.h>

enum { RED = 0x0606, GREEN, BLUE};

typedef struct _Vertex
{
 REAL x,y,z;
#ifdef DIM_4D
 REAL v;
#endif
} Vertex;
typedef struct _Vertex Vector;

typedef struct _Material
{
 REAL c,s,t;
} Material;

typedef struct _TexCoord
{
 REAL x,y;
#ifdef DIM_4D
 REAL v;
#endif
} TexCoord;

typedef struct _Surface
{
 REAL A,B,C,D;
#ifdef DIM_4D
 REAL E;
#endif
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
 REAL caR, caG, caB;
 REAL n_dist;
 int faces;
 int idx;
 int tid;
 int nidx;
 NURBS* nurbs;
 TexCoord nca;
 TexCoord ncb;
 TexCoord ncc;
#ifdef DIM_4D
 Vertex d;
 Vector nd;
 Material mrd, mgd, mbd;
 TexCoord td;
 TexCoord ncd;
#endif
} Triangle;

typedef struct _Screen
{
 unsigned char* pixels;
 int x,y;
#ifdef DIM_4D
 int z;
#endif
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
#ifdef DIM_4D
 REAL minv, maxv;
#endif
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

typedef struct _MaterialTransform
{
 int on_texture;

 int t_id;

 int i_mra, i_mrb, i_mrc;
 int i_mga, i_mgb, i_mgc;
 int i_mba, i_mbb, i_mbc;
 Material mra, mrb, mrc;
 Material mga, mgb, mgc;
 Material mba, mbb, mbc;

 int i_mur, i_mdr;
 int i_mug, i_mdg;
 int i_mub, i_mdb;
 REAL mUR,mDR;
 REAL mUG,mDG;
 REAL mUB,mDB;

 int i_ca;
 REAL caR, caG, caB;
 REAL n_dist;
#ifdef DIM_4D
 int i_mrd, i_mgd, i_mbd;
 Material mrd, mgd, mbd;
#endif

} MaterialTransform;

typedef struct _ListTransform
{
 REAL ** M;
 REAL ** MN;
 int i1, i2;
 int ni1, ni2;
 MaterialTransform mt;
} ListTransform;

typedef struct _BList
{
 struct _BList *next, *prev;
 Box b;
} BList;

typedef struct _Voxel
{
 PtrList* tris;
 Box box;
 int ntbox;
 struct _Voxel* varr;
} Voxel;

typedef struct _IdxTree
{
 struct _IdxTree *s, *h, *t;
 int idx;
} IdxTree;

typedef struct _ThreadRuntime
{
 int idx;
 Ray r;
 Screen* s;
 int thrnum;
} ThreadRuntime;

int run_RT(int lb, char** par);
int load_jpeg_file(unsigned long*** bits, int* x, int* y, FILE* infile);
int save_jpeg_file(unsigned char* bits, int x, int y, FILE* outfile);
int save_gray_jpeg_file(unsigned char* bits, int x, int y, FILE* outfile);
void clear_material_transform(MaterialTransform*);
int debug(char* f, int l, char* fmt, ...);
void panic(char* f, int l);
void spanic(char* f, int l, char* why, ...);
REAL* vector(int siz);
REAL** matrix(int siz);
REAL** matrix2(int siz1, int siz2);
REAL*** matrix3(int siz1, int siz2, int siz3);
void I_matrix(REAL** dst, int siz);
int debugs(char* fmt, ...);
int debug(char* f, int l, char* fmt, ...);
void matrix_mul_vector(REAL* dst, REAL** m, REAL* v, int len);
void free_matrix(REAL** m, int siz);
void free_matrix3(REAL*** m, int siz1, int siz2);
void copy_matrix(REAL** dst, REAL** src, int siz);
void copy_matrix3(REAL*** dst, REAL*** src, int siz1, int siz2, int siz3);
REAL** invert_matrix(REAL** srcC, int dim);
void rotatex(REAL** m, REAL ang);
void rotatey(REAL** m, REAL ang);
void rotatez(REAL** m, REAL ang);
void translatex(REAL** m, REAL arg);
void translatey(REAL** m, REAL arg);
void translatez(REAL** m, REAL arg);
void translate(REAL** m, REAL x, REAL y, REAL z);
void scalex(REAL** m, REAL arg);
void scaley(REAL** m, REAL arg);
void scalez(REAL** m, REAL arg);
void scale(REAL** m, REAL x, REAL y, REAL z);
void init_bmp(BMPTag* b);
void compute_normals(Triangle* t);
void materialize(Material* t);
void translatex(REAL** m, REAL arg);
REAL length(Vector* v);
void transform_triangle(Triangle* t, REAL** m, REAL** nm);
void transform_triangle_material(Triangle* t, MaterialTransform* lt);
void translate_jpeg_to_uchar_format(unsigned long** bits, Texture* t);
void load_jpeg_texture(Texture* t, FILE* texfile);
void free_texture(Texture* t);
void ptrlist_add(PtrList** head, void* t);
void prenurbs_free(REAL**** ptr, int dim, int npts);
REAL*** precompute_nurbs(REAL* knot, REAL t, int dim, int npts);
REAL b0(REAL* knot, REAL t, int i);
REAL fastnurbs(NURBS* nurb, REAL u, REAL v, int xyz);
void fastnurbs_array(NURBS* nurb, REAL u, REAL v, REAL* tab, int siz);
void calc_normal(NURBS* nu, REAL u, REAL v, REAL* x, REAL* y, REAL* z);
int nearly_equal(REAL x, REAL y, REAL eps);
void triangulate_nurbs(NURBS* n, Triangle* tts);
void set_color_uchar(unsigned char* s, int Y, int x, int y, int r, int g, int b);
void set_color(Screen* s, int x, int y, int r, int g, int b);
void get_texturei(Texture* t, REAL x, REAL y, int *r, int* g, int* b);
REAL compute_angle(Vertex* a, Vertex* b, Vertex* c);
int vertex_in_triangle(Triangle* t, Vertex* w);
void ptrlist_add_on_tail(PtrList** tail, void* t);
void ptrlist_add(PtrList** head, void* t);
int intersection_triangle(Triangle* tr, Ray* r, Vertex* ret);
REAL distance(Vertex*, Vertex*);
REAL pseudo_distance(Vertex*, Vertex*);
void normalize_2d(Vector* v);
void ptrlist_free(PtrList** head);
PtrList* ptrlist_get_tail(PtrList* head);
int line_intersect(Vertex* p, Ray* l1, Ray* l2);
void make_ray(Ray* l, Vertex* a, Vertex* b);
REAL pseudo_distance(Vertex* a, Vertex* b);
void normalize(Vector* v);
REAL scalar_prod(Vector* v, Vector* w);
REAL length(Vector* v);
REAL pseudo_length(Vector* v);
int get_red(Screen* s, int x, int y);
int get_green(Screen* s, int x, int y);
int get_blue(Screen* s, int x, int y);
void get_texture(Texture* t, REAL x, REAL y, REAL *r, REAL* g, REAL* b);
void get_texturei(Texture* t, REAL x, REAL y, int *r, int* g, int* b);
REAL calc_surface(Triangle* t);
PtrList* get_head(PtrList* ptr);
void ptrlist_delete_with_tail(PtrList** head, PtrList** tail, PtrList* ptr);
void ptrlist_delete(PtrList** head, PtrList* ptr);
void wrt_bmp(Screen* s, char* out_f);
void downcase(char* b);
void help(char* fn);
void check_options();
void write_config();
void glCube(float x, float y, float z, float v);
void set_light_col(char* pcarg);
void check_light();
void dump_variables(FILE*);
void intersection_btree(Ray* r, BTree* tree, int* l_idx, int* l_glob_sol, Vertex* l_last_int);
void intersection_btree_old_algorithm(VList** head, Ray* r, BTree* tree, int* l_idx, int* l_glob_sol, Vertex* l_last_int);
void compute_slimits(int n);
int args_fprintf(FILE* f, char* fmt, ...);
int args_sprintf(char* f, char* fmt, ...);
int args_fscanf(FILE* f, char* fmt, ...);
int args_sscanf(char* f, char* fmt, ...);
int v_fscanf(FILE* f, char* fmt, ...);
void thr_compute_glpixels_all(unsigned char* glpix);

#endif
