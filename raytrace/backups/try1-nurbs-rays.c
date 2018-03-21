/* Opensource software, license: BSD: author: morgothdbma */
/* morgothdbma@o2.pl, +48693582014 */
/*#define VS*/
#ifdef VS
#define NOGL
#define NOJPEG
#define NOSIGNALS
#define NOGETOPT
#define mfscanf fscanf
#include "getopt.h"
#else
#define mfscanf v_fscanf
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#ifndef VS
#include <unistd.h>
#endif
#include <stdarg.h>
#include <time.h>
#ifdef LINUX
#include <getopt.h>
#include <sys/types.h>
#endif
#ifndef NOGL
#include <GL/glut.h>
#include <pthread.h>
#include <unistd.h>
#endif
#ifndef NOSIGNALS
#include <signal.h>
#endif
#ifndef NOJPEG
#include <jpeglib.h>
#include <setjmp.h>
#define ERR_CANNOTREAD 1
#define ERR_BADJPEG    2
#define ERR_GRAYJPEG   3
#define ERR_256CJPEG   4
#endif
#define HERE __FILE__,__LINE__
#define BACK_R 0x1F
#define BACK_G 0x4F
#define BACK_B 0x7F
#define OPEN_TEXTURE 0
#define OPEN_PARTIAL 1
#define PI 3.14159265
#define REAL long double
#define sREAL double
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
 REAL n_dist;
 int t_id;
} ListTransform;

typedef struct _BList
{
 struct _BList *next, *prev;
 Box b;
} BList;
/* variables */
#ifndef NOGL
pthread_t thread;
#endif
time_t t1,t2;
char demandbmp[1024];
char screenbmp[1024];
char panicbmp[1024];
char debugfile[1024];
char partialbmp[1024];
char signalledbmp[1024];
char textureDir[1024];
char scenef[1024];
char tprefix[1024];
Vertex light, observer;	/* FIXME: more lights ? */
Screen screen;
unsigned int* tids;
#ifndef NOGL
GLfloat rotX,rotY,rotZ;
GLfloat scaX,scaY,scaZ;
GLfloat traX,traY,traZ;
GLfloat ltX,ltY,ltZ;
REAL **pv_m, **pv_mn;
Triangle* pv_tbuf;
int preview_ps;
int pv_usetm;
#endif
int l_enable;
int l_disabled;
int t_disabled;
int sh_disabled;
int invN;
int fps;
int skip_minimalize_algorithm;
int apply_steps;
int g_idx;
int loaded;
int no_sig;
int timeout;
int nTriangles;
int n_pl;
int nNURBS;
int nTriNURBS;
int aa;
int want_bin_pps;
int want_save_pps;
int want_load_pps;
int rt_hlt;
int jqual;
int thr_sx,thr_sy;
int proc_signal;
int want_gjpeg;
int double_res;
int aa_bkup;
int n_dist;
int prev_li;
int hlevel;
unsigned char* glpixels;
int use_jpeg;
int use_gl;
REAL pzoomx, pzoomy;
REAL global_dist;
REAL ambient,minshadow;
REAL maxshadow;
REAL glob_u, glob_v;
REAL glob_rv[Smin];
FILE* dbg;
Texture* texture;
Vertex old_p;
Triangle* g_ts;
NURBS* g_nurbs;
int nTex;
int rseed;
REAL step;
BList* boxes;
BTree* btree;
PtrList* tlist;
int vlight;
int bkup;
int glob_sol;
REAL proc_tr;
REAL proc_tr2;
REAL proc_bt;
int max_rec;
int line_idx;
int ovr_x;
int ovr_y;
int want_bin;
int ovr_b;
int ovr_m;
int ovr_ms;
int ovr_r;
int ovr_a;
int ovr_no;
int norm_dist;
REAL apply_steps_perc;
REAL lookz;
REAL** world;
REAL** worldn;
REAL distorber1, distorber2;
int trans_used;
/*end variables*/
int (*intersection)(Triangle*, Ray*, Vertex*, int*);
void (*get_normal)(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);
/*int nmalloc,nfree;*/
void wrt_bmp(Screen* s, char* out_f);
int intersection_new(Triangle*, Ray*, Vertex*, int*);
int intersection_old(Triangle*, Ray*, Vertex*, int*);
void get_normal_new(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);
void get_normal_old(Vector*, Triangle*, Vertex*, Material*, Material*, Material*, TexCoord*);
void set_color(Screen*, int, int, int, int, int);
#ifndef NOJPEG

struct my_error_mgr
{
 struct jpeg_error_mgr pub;
 jmp_buf setjmp_buffer;
};
typedef struct my_error_mgr * my_error_ptr;
static my_error_ptr errptr = NULL;

static void my_error_exit (j_common_ptr cinfo)
{
 my_error_ptr myerr = (my_error_ptr) cinfo->err;
 (*cinfo->err->output_message) (cinfo);
 longjmp(myerr->setjmp_buffer, 1);
}


void set_rows(int cur_row, JSAMPARRAY pixel_data, int width, unsigned long*** bits)
{
 JSAMPROW ptr;
 int j;
 int r, g, b;
 ptr = pixel_data[0];
 for (j=0;j<width;j++)
   {
    r=GETJSAMPLE(*ptr++);
    g=GETJSAMPLE(*ptr++);
    b=GETJSAMPLE(*ptr++);
    (*bits)[cur_row][j] = (b << 0x10) + (g << 0x8) + r;
    }
}


int load_jpeg_file(unsigned long*** bits, int* x, int* y, FILE* infile)
{
    struct jpeg_decompress_struct cinfo;
    struct my_error_mgr jerr;
    JSAMPARRAY buffer;
    int row_stride;
    int i;
    *x = *y = 0;
    *bits = NULL;
    if (infile == NULL) return ERR_CANNOTREAD;
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    errptr = &jerr;
    if (setjmp(jerr.setjmp_buffer))
       {
	jpeg_destroy_decompress(&cinfo);
/*	fclose(infile);*/
	return ERR_BADJPEG;
       }
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);
    if (cinfo.jpeg_color_space == JCS_GRAYSCALE)
       {
	jpeg_destroy_decompress(&cinfo);
/*	fclose(infile);*/
	return ERR_GRAYJPEG;
       }
    else cinfo.quantize_colors = FALSE;
    jpeg_start_decompress(&cinfo);
    if (cinfo.output_components == 1)
       {
	jpeg_destroy_decompress(&cinfo);
/*	fclose(infile);*/
	return ERR_256CJPEG;
       }
    *x = cinfo.output_width;
    *y = cinfo.output_height;
    *bits = calloc(*y, sizeof(unsigned long*));
    for (i=0;i<*y;i++) (*bits)[i] = calloc(*x, sizeof(unsigned long));
    row_stride = cinfo.output_width * cinfo.output_components;
    buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    i = (*y)-1;
    while (cinfo.output_scanline < cinfo.output_height)
       {
	jpeg_read_scanlines(&cinfo, buffer, 1);
	set_rows(i,buffer,*x, bits);
	i--;
       }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    return 0;
}


void jpeg_compute_pixels(unsigned char* pix)
{
 int i,j,k;
 int sx,sy;
 sx = screen.x;
 sy = screen.y;
 if (!aa)
   {
    for (i=0;i<sy;i++)
    for (j=0;j<sx;j++)
    for (k=0;k<3;k++)
       pix[3*(sx*((sy-i)-1)+j)+k] = screen.pixels[3*(sy*j+i)+k];
   }
 else
   {
    sx /= 2;
    sy /= 2;
    for (i=0;i<sy;i++)
    for (j=0;j<sx;j++)
    for (k=0;k<3;k++)
       pix[3*(sx*((sy-i)-1)+j)+k] = (
	   screen.pixels[3*(2*sy*(2*j)+(2*i))+k]+
	   screen.pixels[3*(2*sy*(2*j+1)+(2*i))+k]+
	   screen.pixels[3*(2*sy*(2*j)+(2*i+1))+k]+
	   screen.pixels[3*(2*sy*(2*j+1)+(2*i+1))+k])/4;
   }
}


void jpeg_compute_gray_pixels(unsigned char* pix)
{
 int i,j;
 int sx,sy;
 sx = screen.x;
 sy = screen.y;
 if (!aa)
   {
    for (i=0;i<sy;i++)
    for (j=0;j<sx;j++)
      {
       pix[sx*((sy-i)-1)+j] =
	   (int)(0.299*(REAL)screen.pixels[3*(sy*j+i)])+
	   (int)(0.587*(REAL)screen.pixels[3*(sy*j+i)+1])+
	   (int)(0.114*(REAL)screen.pixels[3*(sy*j+i)+2]);
      }
   }
 else
   {
    sx /= 2;
    sy /= 2;
    for (i=0;i<sy;i++)
    for (j=0;j<sx;j++)
       pix[sx*((sy-i)-1)+j] = 	/* THAT WAS HARD... */
	   (((int)(0.299*(REAL)screen.pixels[3*(2*sy*(2*j)+(2*i))])+
	   (int)(0.587*(REAL)screen.pixels[3*(2*sy*(2*j)+(2*i))+1])+
	   (int)(0.114*(REAL)screen.pixels[3*(2*sy*(2*j)+(2*i))+2]))+
	   ((int)(0.299*(REAL)screen.pixels[3*(2*sy*(2*j+1)+(2*i))])+
	   (int)(0.587*(REAL)screen.pixels[3*(2*sy*(2*j+1)+(2*i))+1])+
	   (int)(0.114*(REAL)screen.pixels[3*(2*sy*(2*j+1)+(2*i))+2]))+
	   ((int)(0.299*(REAL)screen.pixels[3*(2*sy*(2*j)+(2*i+1))])+
	   (int)(0.587*(REAL)screen.pixels[3*(2*sy*(2*j)+(2*i+1))+1])+
	   (int)(0.114*(REAL)screen.pixels[3*(2*sy*(2*j)+(2*i+1))+2]))+
	   ((int)(0.299*(REAL)screen.pixels[3*(2*sy*(2*j+1)+(2*i+1))])+
	   (int)(0.587*(REAL)screen.pixels[3*(2*sy*(2*j+1)+(2*i+1))+1])+
	   (int)(0.114*(REAL)screen.pixels[3*(2*sy*(2*j)+(2*i))+2])))/4;
   }
}


void init_pixels(unsigned char** pix, int x, int y, int gray)
{
 if (!gray)
   {
    *pix = (unsigned char*)malloc((x*y*3)*sizeof(unsigned char));
    jpeg_compute_pixels(*pix);
   }
 else
   {
    *pix = (unsigned char*)malloc((x*y)*sizeof(unsigned char));
    jpeg_compute_gray_pixels(*pix);
   }
}


void free_pixels(unsigned char** pix)
{
 free(*pix);
 *pix = NULL;
}


int save_jpeg_file(unsigned char* bits, int x, int y, FILE* outfile)
{
 struct jpeg_compress_struct cinfo;
 struct my_error_mgr jerr;
 unsigned char* pixels;
 JSAMPROW row_pointer[1];
 int row_stride;
 pixels = NULL;
 if (aa) { x /= 2; y /= 2; }
/* printf("x:y=%d:%d\n", x,y);*/
 init_pixels(&pixels, x, y, 0);
 cinfo.err = jpeg_std_error(&jerr.pub);
 jerr.pub.error_exit = my_error_exit;
 errptr = &jerr;
 if (setjmp(jerr.setjmp_buffer))
       {
	jpeg_destroy_compress(&cinfo);
	return ERR_BADJPEG;
       }
 jpeg_create_compress(&cinfo);
 jpeg_stdio_dest(&cinfo, outfile);
 cinfo.image_width = x;
 cinfo.image_height = y;
 cinfo.input_components = 3;
 cinfo.in_color_space = JCS_RGB;
/* jpeg_set_default_colorspace(&cinfo);*/
 jpeg_set_defaults(&cinfo);
 jpeg_set_quality(&cinfo, jqual, FALSE);
 jpeg_start_compress(&cinfo, TRUE);
/* jpeg_write_m_header(&cinfo, TRUE);*/
 row_stride = x*3;
 while (cinfo.next_scanline < cinfo.image_height)
   {
    row_pointer[0] = &pixels[cinfo.next_scanline * row_stride];
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
   }
 jpeg_finish_compress(&cinfo);
 jpeg_destroy_compress(&cinfo);
 free_pixels(&pixels);
 return 0;
}


int save_gray_jpeg_file(unsigned char* bits, int x, int y, FILE* outfile)
{
 struct jpeg_compress_struct cinfo;
 struct my_error_mgr jerr;
 unsigned char* pixels;
 JSAMPROW row_pointer[1];
 pixels = NULL;
 if (aa) { x /= 2; y /= 2; }
/* printf("x:y=%d:%d\n", x,y);*/
 init_pixels(&pixels, x, y, 1);
 cinfo.err = jpeg_std_error(&jerr.pub);
 jerr.pub.error_exit = my_error_exit;
 errptr = &jerr;
 if (setjmp(jerr.setjmp_buffer))
       {
	jpeg_destroy_compress(&cinfo);
	return ERR_BADJPEG;
       }
 jpeg_create_compress(&cinfo);
 jpeg_stdio_dest(&cinfo, outfile);
 cinfo.image_width = x;
 cinfo.image_height = y;
 cinfo.input_components = 1;
 cinfo.in_color_space = JCS_GRAYSCALE;
 jpeg_set_defaults(&cinfo);
 jpeg_set_quality(&cinfo, jqual, FALSE);
 jpeg_start_compress(&cinfo, TRUE);
 while (cinfo.next_scanline < cinfo.image_height)
   {
    row_pointer[0] = &pixels[cinfo.next_scanline*x];
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
   }
 jpeg_finish_compress(&cinfo);
 jpeg_destroy_compress(&cinfo);
 free_pixels(&pixels);
 return 0;
}

#endif

int debug(char* f, int l, char* fmt, ...);

void blank_line(Screen* s, int idx)
{
 int i;
 if (!s || idx < 0 || idx>= s->x) 
   {
    debug(HERE, "bad idx: %d", idx);
   }
 for (i=0;i<s->y;i++)
   {
    s->pixels[3*(s->y * idx + i)   ] = 0;
    s->pixels[3*(s->y * idx + i) +1] = 0;
    s->pixels[3*(s->y * idx + i) +2] = 0;
   }
}



void panic(char* f, int l)
{
 printf("RayTracer Engine Panic: file: %s, line: %d\n", f,l);
 blank_line(&screen, line_idx);
 if (screen.x > 0 && screen.y > 0) wrt_bmp(&screen, panicbmp);
#ifdef DEBUG
 debug(HERE, "RayTracer Engine Panic: file: %s, line: %d\n", f,l);
 printf("Entering dead loop, wait for GDB to attach, my pid is: %d\n", getpid());
 if (dbg) fclose(dbg);
 dbg = NULL;
 while (1) sleep(1);
#else
 exit(1);
#endif
}


void spanic(char* why, char* f, int l)
{
 printf("\nFatal error occured: %s\n", why);
 panic(f,l);
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
 if (!dst || siz<= 0) spanic("I_matrix: no dst or bad size",HERE);
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
 if (!dst || !m || !v || len <=0) spanic("matrix_mul_vector: bad parameters",HERE);
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


void copy_matrix(REAL** dst, REAL** src, int siz)
{
 int i,j;
 if (!src || !dst || siz<= 0) return;
 for (i=0;i<siz;i++) for (j=0;j<siz;j++) dst[i][j] = src[i][j];
}


int try_swap(REAL** m, int idx, int dim)
{
 int x;
 for (x=idx;x<dim;x++) if (m[idx][x]) return x;
 return -1;
}


REAL** invert_matrix(REAL** srcC, int dim)
{
 REAL** src, **dst;
 REAL div, pom;
 register int x,k;
 int i,swit;
 REAL* vectmp;
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
       if (swit < 0) spanic("invert_matrix: cant do that",HERE); /*uninvertable matrix*/
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
 free_matrix(src, dim);
 free(vectmp);
 return dst;
}


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


void print_matrix(REAL** m, int siz)
{
 int i,j;
 if (!m || siz <=0) spanic("print_matrix: is null or dim < 0",HERE);
 printf("matrix %p\n", (void*)m);
 for (i=0;i<siz;i++)
   {
     for (j=0;j<siz;j++) printf("%3.4Lf\t", m[i][j]);
    printf("\n");
   }
}


void init_bmp(BMPTag* b)	
{
 int i;
 if (!b) spanic("init_bmp: BMPTag is null",HERE);
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
 int i,n;
 if (!scr) spanic("init_screen: scr is null",HERE);
 scr->pixels = (unsigned char*)malloc((x*y*3+1)*sizeof(unsigned char));
 if (!scr->pixels) spanic("init_screen: malloc pixels failed",HERE);
 scr->x = x;
 scr->y = y;
 n = x * y * 3;
 for (i=0;i<=n;i++) scr->pixels[i] = 0;
}


void free_screen(Screen* scr)
{
 if (!scr) spanic("free_screen: no screen to be freed",HERE);
 free(scr->pixels);
 scr->x = scr->y = 0;
}


void print_ray(Ray* r)
{
 printf("FROM: (%Lf,%Lf,%Lf), DIR: (%Lf,%Lf,%Lf), REC(%d), COL(%d)\n",
	 r->P.x, r->P.y, r->P.z, r->d.x, r->d.y, r->d.z, r->r, r->c);
}


int read_vector(FILE* f, Vector* v)
{
 int nr;
 if (!f || !v) spanic("read_vector: bad params",HERE);
 nr = fscanf(f, "Vector: (%Lf,%Lf,%Lf)", &v->x, &v->y, &v->z);
 if (nr == 3) return 1;
 else return 0;
}


int read_vertex(FILE* f, Vertex* v)
{
 int nr;
 if (!f || !v) spanic("read_vertex: bad params",HERE);
 nr = fscanf(f, "Vertex: (%Lf,%Lf,%Lf)", &v->x, &v->y, &v->z);
/* printf("nr = %d\n", nr); */
 if (nr == 3) return 1;
 else return 0;
}


int read_texcoord(FILE* f, TexCoord* t)
{
 int nr;
 if (!f || !t) spanic("read_texcoord: bad params",HERE);
 nr = fscanf(f, "TexCoord: (%Lf,%Lf)", &t->x, &t->y);
/* printf("%d %f,%f\n", nr, t->x, t->y);*/
 if (nr != 2) return 0;
 if (t->x < 0. || t->x > 1.) return 0;
 if (t->y < 0. || t->y > 1.) return 0;
 return 1;
}

void normalize(Vector*);

void normalize_color_factors(Triangle* t)
{
 REAL f;
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

REAL length(Vector* v);

void transform_observer(Vertex* t, REAL** m, REAL** nm)
{
 REAL *x, *rx;
 if (!t || !m || !nm) spanic("transform_observer: no observer or matrix", HERE);
 x  = vector(4);
 rx = vector(4);
 x[0] = t->x;
 x[1] = t->y;
 x[2] = t->z;
 x[3] = 1.;
 matrix_mul_vector(rx, m, x, 4);
 t->x = rx[0];
 t->y = rx[1];
 t->z = rx[2];
 free(x);
 free(rx);
}


void transform_light(void* t, REAL** m, REAL** nm, int is_vec)
{
 REAL *x, *rx;
 REAL **M;
 Vector* v;
 if (is_vec) { M = nm; v = (Vector*)t; }
 else        { M = m;  v = (Vertex*)t; }
 if (!t || !m || !nm) spanic("transform_light: no light or matrix", HERE);
 /*print_matrix(M, 4);
 print_matrix(m, 4);
 print_matrix(nm, 4);*/
 x  = vector(4);
 rx = vector(4);
 x[0] = v->x;
 x[1] = v->y;
 x[2] = v->z;
 x[3] = 1.;
 matrix_mul_vector(rx, M, x, 4);
 v->x = rx[0];
 v->y = rx[1];
 v->z = rx[2];
/* printf("%Lf,%Lf,%Lf\n", v->x, v->y, v->z);*/
 free(x);
 free(rx);
}


void roll_back_transform_triangle(Triangle* t, REAL** mi, REAL** nmi)
{
 REAL *x, *rx;
 REAL **m, **nm;
 m = invert_matrix(mi, 4);
 nm = invert_matrix(nmi, 4);
 if (!t || !m || !nm) spanic("roll_back_transform_triangle: no triangle or matrix", HERE);
 x  = vector(4);
 rx = vector(4);
 x[0] = t->a.x;
 x[1] = t->a.y;
 x[2] = t->a.z;
 x[3] = 1.;
 matrix_mul_vector(rx, m, x, 4);
 t->a.x = rx[0];
 t->a.y = rx[1];
 t->a.z = rx[2];
 x[0] = t->b.x;
 x[1] = t->b.y;
 x[2] = t->b.z;
 x[3] = 1.;
 matrix_mul_vector(rx, m, x, 4);
 t->b.x = rx[0];
 t->b.y = rx[1];
 t->b.z = rx[2];
 x[0] = t->c.x;
 x[1] = t->c.y;
 x[2] = t->c.z;
 x[3] = 1.;
 matrix_mul_vector(rx, m, x, 4);
 t->c.x = rx[0];
 t->c.y = rx[1];
 t->c.z = rx[2];
 x[0] = t->na.x;
 x[1] = t->na.y;
 x[2] = t->na.z;
 x[3] = 1.;
 matrix_mul_vector(rx, nm, x, 4);
 t->na.x = rx[0];
 t->na.y = rx[1];
 t->na.z = rx[2];
 x[0] = t->nb.x;
 x[1] = t->nb.y;
 x[2] = t->nb.z;
 x[3] = 1.;
 matrix_mul_vector(rx, nm, x, 4);
 t->nb.x = rx[0];
 t->nb.y = rx[1];
 t->nb.z = rx[2];
 x[0] = t->nc.x;
 x[1] = t->nc.y;
 x[2] = t->nc.z;
 x[3] = 1.;
 matrix_mul_vector(rx, nm, x, 4);
 t->nc.x = rx[0];
 t->nc.y = rx[1];
 t->nc.z = rx[2];
 normalize(&t->na);
 normalize(&t->nb);
 normalize(&t->nc);
 compute_surface(t);
 free(x);
 free(rx);
 free_matrix(m, 4);
 free_matrix(nm, 4);
}


void transform_triangle(Triangle* t, REAL** m, REAL** nm)
{
 REAL *x, *rx;
 if (!t || !m || !nm) spanic("transform_triangle: no triangle or matrix", HERE);
 x  = vector(4);
 rx = vector(4);
 x[0] = t->a.x;
 x[1] = t->a.y;
 x[2] = t->a.z;
 x[3] = 1.;
 matrix_mul_vector(rx, m, x, 4);
 t->a.x = rx[0];
 t->a.y = rx[1];
 t->a.z = rx[2];
 x[0] = t->b.x;
 x[1] = t->b.y;
 x[2] = t->b.z;
 x[3] = 1.;
 matrix_mul_vector(rx, m, x, 4);
 t->b.x = rx[0];
 t->b.y = rx[1];
 t->b.z = rx[2];
 x[0] = t->c.x;
 x[1] = t->c.y;
 x[2] = t->c.z;
 x[3] = 1.;
 matrix_mul_vector(rx, m, x, 4);
 t->c.x = rx[0];
 t->c.y = rx[1];
 t->c.z = rx[2];
 x[0] = t->na.x;
 x[1] = t->na.y;
 x[2] = t->na.z;
 x[3] = 1.;
 matrix_mul_vector(rx, nm, x, 4);
 t->na.x = rx[0];
 t->na.y = rx[1];
 t->na.z = rx[2];
 x[0] = t->nb.x;
 x[1] = t->nb.y;
 x[2] = t->nb.z;
 x[3] = 1.;
 matrix_mul_vector(rx, nm, x, 4);
 t->nb.x = rx[0];
 t->nb.y = rx[1];
 t->nb.z = rx[2];
 x[0] = t->nc.x;
 x[1] = t->nc.y;
 x[2] = t->nc.z;
 x[3] = 1.;
 matrix_mul_vector(rx, nm, x, 4);
 t->nc.x = rx[0];
 t->nc.y = rx[1];
 t->nc.z = rx[2];
 normalize(&t->na);
 normalize(&t->nb);
 normalize(&t->nc);
 compute_surface(t);
 free(x);
 free(rx);
}

void read_transformation_long(FILE* f, REAL*** mm, REAL*** mn,int* tuse, REAL* nor_dist, int* texid);

void read_transformation(FILE* f, REAL*** mm, REAL*** mn, int* tuse, REAL* nor_dist)
{
 int tmpl;
 read_transformation_long(f, mm, mn, tuse, nor_dist, &tmpl);
}


void roll_back_list_transform_triangle(Triangle* t)
{
 PtrList* tmp;
 ListTransform* lt;
/* printf("roll back list transform\n");*/
 tmp = tlist;
 while (tmp->next) tmp = tmp->next;
 while (tmp)
   {
    lt = (ListTransform*)tmp->ptr;
    if (t->idx >= lt->i1 && t->idx <= lt->i2)
      {
/*	  printf("removing: %p\n", (void*)tmp);*/
 /*      printf("roll_back transform list: \n");
       print_matrix(lt->M, 4);
       print_matrix(lt->MN, 4);*/
       roll_back_transform_triangle(t, lt->M, lt->MN);
       if (lt->n_dist >= 0.)
         {
	  t->n_dist = lt->n_dist;
/*	  printf("New n_dist is: %Lf\n", t->n_dist);*/
	 }
       if (lt->t_id >= 0)
         {
	  t->tid = lt->t_id;
	 }
      }
    tmp = tmp->prev;
   }
}


void list_transform_triangle(Triangle* t)
{
 PtrList* tmp;
 ListTransform* lt;
 tmp = tlist;
/* printf("list transform\n");*/
 while (tmp)
   {
    lt = (ListTransform*)tmp->ptr;
    if (t->idx >= lt->i1 && t->idx <= lt->i2)
      {
/*	  printf("applying: %p\n", (void*)tmp);*/
/*       printf("transform list: \n");
       print_matrix(lt->M, 4);
       print_matrix(lt->MN, 4);*/
       transform_triangle(t, lt->M, lt->MN);
       if (lt->n_dist >= 0.)
         {
	  t->n_dist = lt->n_dist;
/*	  printf("New n_dist is: %Lf\n", t->n_dist);*/
	 }
       if (lt->t_id >= 0)
         {
	  t->tid = lt->t_id;
	 }
      }
    tmp = tmp->next;
   }
}

int v_fscanf(FILE* f, char* fmt, ...);

int read_triangle(FILE* f, Triangle* t, int* idx, int* t_use, char** tmap, Triangle* all_t, int rd)
{
 int nr,n,ti;
 int t_used;
 int i,id,ok;
 int rw,rlt,cpt;
 REAL **trmat, **trmatn;
 REAL dd;
 char str[512];
 long pos;
 if (!f || !t) spanic("read_triangle: no triangle or file not opened.",HERE);
 t_used = 0;
 trmat = trmatn = NULL;
 if (!rd)  goto finalize;
 pos = ftell(f);
 rw = rlt = 1;
 nr = mfscanf(f, "CopyTriangles: dst=%d,src=%d,num=%d\n", &i, &ti,&n);
 cpt = 0;
 if (nr != 3)
   {
    fseek(f, pos, SEEK_SET);
    nr = mfscanf(f, "CopyTrianglesAdv: dst=%d,src=%d,num=%d,rollW=%d,rollLT=%d\n", &i, &ti,&n, &rw, &rlt);
    if (nr != 5) fseek(f, pos, SEEK_SET);
    else cpt = 1;
   }
 else cpt = 1;
 if (cpt)
   {
/*       printf("i = %d, ti = %d, idx = %d, nt = %d\n", i, ti, idx, nTriangles);*/
    if (n < 0) spanic("read_triange: mcopy: bad numeric count", HERE);
    if (i != *idx) spanic("read_triangle: mcopy: destination index bad", HERE);
    if (ti + (n-1) >= i) spanic("read_triangle: mcopy: would overlap indices", HERE);
    if (ti < 0 || ti >= nTriangles-(n-1)) spanic("read_triangle: source index bad", HERE);
/*    printf("running mcopy... src = %d, dst = %d, num = %d\n", ti, i, n);*/
    id = *idx;
/*    printf("id = %d\n", id);*/
    for (i=id;i<id+n;i++)
      {
/*       printf(" i = %d, ti = %d, id = %d, memcpy %d < %d\n", i, ti, id,  i, ti+i-id);*/
       memcpy(&all_t[i], &all_t[ti+i-id], sizeof(Triangle));
/*       printf("setting temprary index to rollback: %d\n", ti+i-id);*/
       all_t[i].idx = ti+i-id;
       if (trans_used && rw)
         {
          roll_back_transform_triangle(&all_t[i], world, worldn);
         }
       if (rlt) roll_back_list_transform_triangle(&all_t[i]);
/*       printf("setting to final idx to apply trans: %d\n", i);*/
       all_t[i].idx = i;
       list_transform_triangle(&all_t[i]);
       if (trans_used)
         {
          transform_triangle(&all_t[i], world, worldn);
         }
       if (all_t[i].tid != 0)
         {
          t_use[all_t[i].tid-1] ++;
         }
      }
/*    printf("increasing idx from %d to %d\n", *idx, (*idx)+(n-1));*/
    *idx += (n-1);
    return 1;
   }
 nr = mfscanf(f, "CopyTriangle: %d<-%d\n", &i, &ti);
 if (nr != 2) fseek(f, pos, SEEK_SET);
 else
   {
/*       printf("i = %d, ti = %d, idx = %d, nt = %d\n", i, ti, idx, nTriangles);*/
    if (ti == i) spanic("read_triangle: copy indices the same", HERE);
    if (ti > i)  spanic("read_triangle: cannot copy from unloaded triangle", HERE);
    if (i != *idx) spanic("read_triangle: destination index bad", HERE);
    if (ti < 0 || ti >= nTriangles) spanic("read_triangle: source index bad", HERE);
    memcpy(t, &all_t[ti], sizeof(Triangle));
    t->idx = ti;
    if (trans_used)
      {
/*       printf("world roll back\n");*/
       roll_back_transform_triangle(t, world, worldn);
      }
    roll_back_list_transform_triangle(t);
    t->idx = *idx;
    list_transform_triangle(t);
    if (trans_used)
      {
/*       printf("world again\n");*/
       transform_triangle(t, world, worldn);
      }
    if (t->tid != 0)
      {
       t_use[t->tid-1] ++;
      }
/*    printf("normal (%Lf,%Lf,%Lf)\n", t->na.x, t->na.y, t->na.z);*/
    return 1;
   }
 nr = mfscanf(f, "Triangle: %d\n", &n);
/* printf("nr = %d\n", nr);*/
 if (nr != 1) spanic("read_triangle: cant read idx",HERE);
 if (n < 0) spanic("read_triangle: idx < 0", HERE);
 mfscanf(f, "{\n");
 t->idx = *idx;
 t->nurbs = NULL;
 t->nidx = -1;
 t->nca.x = t->nca.y = 0.;
 t->ncb.x = t->ncb.y = 0.;
 t->ncc.x = t->ncc.y = 0.;
 mfscanf(f, " a: ");
 if (!read_vertex(f, &t->a)) spanic("read_triangle: cant read a vertex",HERE);
 mfscanf(f,"\n");
/* printf("(%Lf,%Lf,%Lf)\n", t->a.x, t->a.y, t->a.z);*/
 /*t->a.x += distorber1;
 t->a.y += distorber2;*/
 mfscanf(f, " b: ");
 if (!read_vertex(f, &t->b)) spanic("read_triangle: cant read b vertex",HERE);
 mfscanf(f,"\n");
 /*t->b.x += distorber1;
 t->b.y += distorber2;*/
 mfscanf(f, " c: ");
 if (!read_vertex(f, &t->c)) spanic("read_triangle: cant read c vertex", HERE);
 mfscanf(f,"\n");
 /*t->c.x += distorber1;
 t->c.y += distorber2;*/
 mfscanf(f, " texA: ");
 if (!read_texcoord(f, &t->ta)) spanic("read_triangle: cant read a texcoord", HERE);
 mfscanf(f,"\n");
 mfscanf(f, " texB: ");
 if (!read_texcoord(f, &t->tb)) spanic("read_triangle: cant read b texcoord", HERE);
 mfscanf(f,"\n");
 mfscanf(f, " texC: ");
 if (!read_texcoord(f, &t->tc)) spanic("read_triangle: cant read c texcoord", HERE);
 mfscanf(f,"\n");
 mfscanf(f, " na: ");
 if (!read_vector(f, &t->na)) spanic("read_triangle: cant read a normal", HERE);
 mfscanf(f,"\n");
 mfscanf(f, " nb: ");
 if (!read_vector(f, &t->nb)) spanic("read_triangle: cant read b normal", HERE);
 mfscanf(f,"\n");
 mfscanf(f, " nc: ");
 if (!read_vector(f, &t->nc)) spanic("read_triangle: cant read c normal", HERE);
 finalize:
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
 if (!rd) 
  {
   t->s.A = t->s.B = t->s.C = t->s.D = 0.;
   t->idx = *idx;
   goto finalize2;
  }
 mfscanf(f,"\n");
 nr = mfscanf(f," %s RGB: (%Lf,%Lf,%Lf)\n", str, &t->mra.t, &t->mga.t, &t->mba.t);
/* printf("nr = %d\n", nr);*/
 if (nr != 4) spanic("read_triangle: cant read material transparency (all or a only)", HERE);
 if (!strcmp(str, "transparency:"))
   {
    nr = mfscanf(f," specular: RGB: (%Lf,%Lf,%Lf)\n", &t->mra.s, &t->mga.s, &t->mba.s);
    if (nr != 3) spanic("read_triangle: cant read material specular", HERE);
    nr = mfscanf(f," diffuse: RGB: (%Lf,%Lf,%Lf)\n", &t->mra.c, &t->mga.c, &t->mba.c);
    if (nr != 3) spanic("read_triangle: cant read material diffuse", HERE);
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
    nr = mfscanf(f," specularA: RGB: (%Lf,%Lf,%Lf)\n", &t->mra.s, &t->mga.s, &t->mba.s);
/*    printf("nr = %d\n", nr);*/
    if (nr != 3) spanic("read_triangle: cant read a material specular", HERE);
    nr = mfscanf(f," diffuseA: RGB: (%Lf,%Lf,%Lf)\n", &t->mra.c, &t->mga.c, &t->mba.c);
    if (nr != 3) spanic("read_triangle: cant read a material diffuse", HERE);
    nr = mfscanf(f," transparencyB: RGB: (%Lf,%Lf,%Lf)\n", &t->mrb.t, &t->mgb.t, &t->mbb.t);
    if (nr != 3) spanic("read_triangle: cant read b material transparency", HERE);
    nr = mfscanf(f," specularB: RGB: (%Lf,%Lf,%Lf)\n", &t->mrb.s, &t->mgb.s, &t->mbb.s);
    if (nr != 3) spanic("read_triangle: cant read b material specular", HERE);
    nr = mfscanf(f," diffuseB: RGB: (%Lf,%Lf,%Lf)\n", &t->mrb.c, &t->mgb.c, &t->mbb.c);
    if (nr != 3) spanic("read_triangle: cant read b material diffuse", HERE);
    nr = mfscanf(f," transparencyC: RGB: (%Lf,%Lf,%Lf)\n", &t->mrc.t, &t->mgc.t, &t->mbc.t);
    if (nr != 3) spanic("read_triangle: cant read c material transparency", HERE);
    nr = mfscanf(f," specularC: RGB: (%Lf,%Lf,%Lf)\n", &t->mrc.s, &t->mgc.s, &t->mbc.s);
    if (nr != 3) spanic("read_triangle: cant read c material specular", HERE);
    nr = mfscanf(f," diffuseC: RGB: (%Lf,%Lf,%Lf)\n", &t->mrc.c, &t->mgc.c, &t->mbc.c);
    if (nr != 3) spanic("read_triangle: cant read c material diffuse", HERE);
   }
 else spanic("read_triangle: bad keyword: expected: transparency(A):", HERE);
 pos = ftell(f);
 nr = mfscanf(f," surface: ABCD: (%Lf,%Lf,%Lf,%Lf)\n", &t->s.A, &t->s.B, &t->s.C, &t->s.D);
 if (nr != 4)
   {
/*    spanic("read_triangle: cant read surface", HERE);*/
    t->s.A = t->s.B = t->s.C = t->s.D = 0.;
    fseek(f, pos, SEEK_SET);
   }
 finalize2:
 normalize_color_factors(t);
 if (t->s.A == 0. && t->s.B == 0. && t->s.C == 0. && t->s.D == 0.) compute_surface(t);
 if (!rd) goto finalize3;
 nr = mfscanf(f, " %s (%Lf,%Lf)\n", str, &t->mUR, &t->mDR);
 if (nr != 3) spanic("read_triangle: cant read mU/mD factors (all or r)", HERE);
 if (!strcmp(str, "transparencyFact:"))
   {
    t->mUG = t->mUR;
    t->mUB = t->mUR;
    t->mDG = t->mDR;
    t->mDB = t->mDR;
   }
 else if (!strcmp(str, "transparencyFactR:"))
   {
    nr = mfscanf(f, " transparencyFactG: (%Lf,%Lf)\n", &t->mUG, &t->mDG);
    if (nr != 2) spanic("read_triangle: cant read g mU/mD factors", HERE);
    nr = mfscanf(f, " transparencyFactB: (%Lf,%Lf)\n", &t->mUB, &t->mDB);
    if (nr != 2) spanic("read_triangle: cant read b mU/mD factors", HERE);
   }
 else spanic("read_triangle: bad keyword: expected: transparencyFact(R):", HERE);
 if (t->mUR < .05 || t->mUR > 20.) spanic("read_triangle: bad mU/mD r value", HERE);
 if (t->mDR < .05 || t->mDR > 20.) spanic("read_triangle: bad mU/mD r value", HERE);
 if (t->mUG < .05 || t->mUG > 20.) spanic("read_triangle: bad mU/mD g value", HERE);
 if (t->mDG < .05 || t->mDG > 20.) spanic("read_triangle: bad mU/mD g value", HERE);
 if (t->mUB < .05 || t->mUB > 20.) spanic("read_triangle: bad mU/mD b value", HERE);
 if (t->mDB < .05 || t->mDB > 20.) spanic("read_triangle: bad mU/mD b value", HERE);
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr !=1) spanic("raed_triangle: cant read specularFact or normalDist", HERE);
 if (!strcmp(str, "normalDist:"))
   {
    nr = mfscanf(f, " %Lf%%\n", &t->n_dist);
    t->n_dist /= 100.;
    if (t->n_dist < 0.) spanic("read_triangle: negative n_dist", HERE);
   }
 else { t->n_dist = global_dist; fseek(f, pos, SEEK_SET); }
/* printf("ndist= %Lf\n", t->n_dist);*/
 nr = mfscanf(f, " specularFact: %Lf\n", &t->ca);
 if (nr != 1) spanic("read_triangle: cant read specularFact", HERE);
 if (t->ca < -1.) spanic("read_triangle: specularFact < 0", HERE);
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr !=1) spanic("raed_triangle: truncated triangle: cant read faces", HERE);
 if (!strcmp(str, "faces:"))
   {
    nr = mfscanf(f, " %d\n", &t->faces);
   }
 else { t->faces = 1; fseek(f, pos, SEEK_SET); }
 if (nr != 1) spanic("read_triangle: cant read faces", HERE);
 if (t->faces != 1 && t->faces != 2) spanic("read_triangle: bad faces count: expected 1 or 2", HERE);
 if (!tmap)
   {
    nr = mfscanf(f, " texture: %d\n", &ti);
    if (nr != 1) spanic("read_triangle: cant read texture", HERE);
   }
 else
   {
    pos = ftell(f);
    nr = mfscanf(f, " texture: %d\n", &ti);
    if (nr <= 0)
      {
       fseek(f, pos, SEEK_SET);
       nr = mfscanf(f, " texture: %s\n", str);
       if (nr != 1) spanic("read_triangle: cant read texture (with name mapping)", HERE);
       ok = 0;
       ti = 0;
       for (i=0;i<nTex;i++) if (tmap[i] && !strcmp(str, tmap[i]))
         {
	  ti = i+1;
	  ok = 1;
	  break;
	 }
/*       printf("str = %s\n", str);*/
       if (!ok) spanic("read_triangle: unknown texture mapping name", HERE);
      }
   }
/* if (n < 0) spanic("read_triangle: texture id < 0", HERE);*/
 finalize3:
 if (!rd) ti = t->tid;
/* printf("tid = %d\n", ti);*/
 if (ti < 0 || ti > nTex) spanic("read_triangle: texture id out of range", HERE);
 if (ti != 0)
   {
    /*t->tex = &texture[ti-1];*/
    t->tid = ti;
    t_use[ti-1] ++;
   }
 else
   {
    t->tex = NULL;
    t->tid = 0;
   }
 if (!rd) goto finalize4;
 pos = ftell(f);
 nr = mfscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "Transform:"))
    {
     trmat = matrix(4);
     I_matrix(trmat, 4);
     trmatn = matrix(4);
     I_matrix(trmatn, 4);
     read_transformation_long(f, &trmat, &trmatn, &t_used, &dd, &id);
     if (dd >= 0.) t->n_dist = dd;
     if (id >= 0.) t->tid = id;
    }
 else fseek(f, pos, SEEK_SET);
 mfscanf(f,"}\n");
/* printf("transformations:\n");*/
/* printf("normal (%Lf,%Lf,%Lf)\n", t->na.x, t->na.y, t->na.z);*/
 if (t_used)
   {
    /*printf("transform internal: \n");
    print_matrix(trmat, 4);
    print_matrix(trmatn, 4);*/
    transform_triangle(t, trmat, trmatn);
    free_matrix(trmat, 4);
    free_matrix(trmatn, 4);
   }
 finalize4:
 list_transform_triangle(t);
 if (trans_used)
   {
/*       printf("world transform\n");*/
    /*printf("transform world: \n");
    print_matrix(world, 4);
    print_matrix(worldn, 4);*/
    transform_triangle(t, world, worldn);
   }
/* printf("Triangle %d (%Lf,%Lf,%Lf) (%Lf,%Lf,%Lf) (%Lf,%Lf,%Lf)\n", *idx, t->a.x, t->a.y, t->a.z, t->b.x, t->b.y, t->b.z, t->c.x, t->c.y, t->c.z);*/
 /*printf("texA: (%Lf,%Lf)\n", (t->a.x + 175.) / 350., (t->a.y + 150.) / 350.);
 printf("texB: (%Lf,%Lf)\n", (t->b.x + 175.) / 350., (t->b.y + 150.) / 350.);
 printf("texC: (%Lf,%Lf)\n", (t->c.x + 175.) / 350., (t->c.y + 150.) / 350.);*/
/* printf(" idx = %d\n", idx);*/
/* printf("normal (%Lf,%Lf,%Lf)\n", t->na.x, t->na.y, t->na.z);*/
 t->a.x += distorber1;
 t->a.y += distorber2;
 t->b.x += distorber1;
 t->b.y += distorber2;
 t->c.x += distorber1;
 t->c.y += distorber2;
 t->s.D -= t->s.A*distorber1 + t->s.B*distorber2;
 if (t->tid != 0)
   {
    t_use[t->tid-1] ++;
   }
/* if (t->tex) printf("texture = %p(%d,%d)\n", t->tex->pixels, t->tex->x, t->tex->y);*/
 return 1;
}

#ifndef NOJPEG

void translate_jpeg_to_uchar_format(unsigned long** bits, Texture* t)
{
 int i,j;
 if (!bits) spanic("translate_jpeg_to_uchar_format: no jbits", HERE);
/* printf("x:y = %d:%d\n", t->x, t->y);*/
 for (i=0;i<t->x;i++)
 for (j=0;j<t->y;j++)
     {
      t->pixels[3*(t->y*i+j)  ] = bits[j][i] & 0xFF;
      t->pixels[3*(t->y*i+j)+1] = (bits[j][i] & 0xFF00) >> 0x8;
      t->pixels[3*(t->y*i+j)+2] = (bits[j][i] & 0xFF0000) >> 0x10;
     }
}


void load_jpeg_texture(Texture* t, FILE* texfile)
{
 unsigned long** tex_data;
 int err,i;
 int tex_x, tex_y;
 tex_data = NULL;
 err = load_jpeg_file(&tex_data, &tex_x, &tex_y, texfile);
 fclose(texfile);
 if (err) spanic("load_jpeg_texture: jpeg decompress error", HERE);
 t->x = tex_x;
 t->y = tex_y;
 t->pixels = (unsigned char*)malloc(3*tex_x*tex_y*sizeof(unsigned char));
 translate_jpeg_to_uchar_format(tex_data, t);
 for (i=0;i<tex_y;i++)
   {
    if (tex_data[i]) free((void*)(tex_data[i]));
    (tex_data[i]) = 0;
   }
 if (tex_data) free((void*)tex_data);
/* *tex_data = 0;*/
}

#endif

void create_texture(Texture* t, char* dir, int n, int mode)
{
 BMPTag bm_handle;
 FILE* plik;
 char fn[1024];
 int i,j;
 char r,g,b,m;
 plik = NULL;
 init_bmp(&bm_handle);
 if (mode == OPEN_TEXTURE)
  {
   sprintf(fn, "%s./%s/%d.bmp", tprefix, dir, n+1);
   printf("Opening texture file: %s\n", fn);
  }
 else if (mode == OPEN_PARTIAL)
  {
   strcpy(fn, dir);
   printf("Recovering from partial file: %s\n", fn);
  }
/* printf("%s %d && %d && (%p || %p)\n", fn, mode==OPEN_PARTIAL,use_jpeg, strstr(fn,"jpg"), strstr(fn,"jpeg"));*/
 if (mode == OPEN_PARTIAL && use_jpeg && (strstr(fn, "jpg") || strstr(fn, "jpeg")))
 {
     printf("Opening BACKUP from JPEG is REALLY NOT recomended\nSpecial marks in file can be destroyed by compress/decompress process\n");
     printf("YOU MAY EXPECT WRONG RESULTS; YOU HAVE BEEN WARNED!\n");
     plik = NULL;
 }
 else plik = fopen(fn,"rb");
 if (!plik)
   {
    if (!use_jpeg)
      {
       printf("Error opening BMP: %s\n", fn);
       spanic("create_texture: cant open file", HERE);
      }
    else
      {
       if (mode == OPEN_TEXTURE)
         {
          sprintf(fn, "./%s/%d.jpeg", dir, n+1);
          printf("Opening texture file: %s\n", fn);
         }
      plik = fopen(fn,"rb");
      if (!plik) spanic("create_texture: cant open jpeg file", HERE);
      }
#ifndef NOJPEG
    load_jpeg_texture(t, plik);
#endif
 /*   fclose(plik);*/
    return;
   }
 i = fscanf(plik,"%c%c",&b,&m);
 if (i != 2) spanic("create_texture: truncated file", HERE);
 if (b != 'B' || m != 'M') spanic("create_texture: not a BMP file", HERE);
 fread(&bm_handle.fsize,4,1,plik);
 fread(&bm_handle.dummy,4,1,plik);
 fread(&bm_handle.offset,4,1,plik);
 fread(&bm_handle.dummy2,4,1,plik);
 fread(&bm_handle.bm_x,4,1,plik);
 fread(&bm_handle.bm_y,4,1,plik);
 fread(&bm_handle.planes,2,1,plik);
 fread(&bm_handle.bpp,2,1,plik);
 if (bm_handle.bpp != 24) spanic("read_triangle: only 24BPP BMPs suported", HERE);
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
 /*int i;*/
 if (!strcmp(dir, "")) spanic("create_textures: no texture directory", HERE);
 if (nTex == 0) { texture = NULL; return; }
 texture = (Texture*)malloc(nTex*sizeof(Texture));
 if (!texture) spanic("create_textures: malloc failed", HERE);
 /*for (i=0;i<nTex;i++) create_texture(&texture[i], dir, i, OPEN_TEXTURE);*/
}


int is_binary(FILE* f)
{
 char str[64];
 fread(str,4,1,f);
 str[4] = 0;
 fseek(f, 0, SEEK_SET);
 if (!strcmp(str,"BINT")) return 1;
 else return 0;
}


void postprocess_textures(int* t_use, Triangle* t, char* dir)
{
 int i;
 for (i=0;i<nTex;i++) if (t_use[i] > 0) create_texture(&texture[i], dir, i, OPEN_TEXTURE);
 else texture[i].pixels = NULL;
 for (i=0;i<nTriangles;i++)
   {
    if (t[i].tid > 0) t[i].tex = &texture[t[i].tid - 1];
    else t[i].tex = NULL;
   }
}

void reconstruct_nurbs_addr(Triangle* t)
{
 int i;
 for (i=0;i<nNURBS;i++)
   {
    if (g_nurbs[i].idx == t->nidx) 
      {
       t->nurbs = &g_nurbs[i];
       return;
      }
   }
 spanic("cannot find NURBS idx for current triangle", HERE);
}

void load_binary_nurb(FILE* f, NURBS* n)
{
 sREAL tmp;
 int i,j,k;
 fread(&n->idx, sizeof(int), 1, f);
 fread(&n->ntri, sizeof(int), 1, f);
 fread(&n->itex, sizeof(int), 1, f);
 fread(&n->faces, sizeof(int), 1, f);
 fread(&n->tid, sizeof(int), 1, f);
 fread(&n->n1, sizeof(int), 1, f);
 fread(&n->p1, sizeof(int), 1, f);
 fread(&n->m1, sizeof(int), 1, f);
 fread(&n->d1, sizeof(int), 1, f);
 fread(&n->n2, sizeof(int), 1, f);
 fread(&n->p2, sizeof(int), 1, f);
 fread(&n->m2, sizeof(int), 1, f);
 fread(&n->d2, sizeof(int), 1, f);
 n->knot1 = (REAL*)malloc((n->m1+1)*sizeof(REAL));
 n->knot2 = (REAL*)malloc((n->m2+1)*sizeof(REAL));
 n->D = matrix3(n->n1, n->n2, S);
 n->w = matrix2(n->n1, n->n2);
 for (i=0;i<=n->m1;i++)
   {
    fread(&tmp, sizeof(sREAL), 1, f);
    n->knot1[i] = tmp;
   }
 for (i=0;i<=n->m2;i++)
   {
    fread(&tmp, sizeof(sREAL), 1, f);
    n->knot2[i] = tmp;
   }
 for (i=0;i<n->n1;i++)
 for (j=0;j<n->n2;j++)
   {
    fread(&tmp, sizeof(sREAL), 1, f);
    n->w[i][j] = tmp;
   }
 for (i=0;i<n->n1;i++)
 for (j=0;j<n->n2;j++)
 for (k=0;k<S;k++)
   {
    fread(&tmp, sizeof(sREAL), 1, f);
    n->D[i][j][k] = tmp;
   }
}

void save_binary_nurb(FILE* f, NURBS* n)
{
 sREAL tmp;
 int i,j,k;
 fwrite(&n->idx, sizeof(int), 1, f);
 fwrite(&n->ntri, sizeof(int), 1, f);
 fwrite(&n->itex, sizeof(int), 1, f);
 fwrite(&n->faces, sizeof(int), 1, f);
 fwrite(&n->tid, sizeof(int), 1, f);
 fwrite(&n->n1, sizeof(int), 1, f);
 fwrite(&n->p1, sizeof(int), 1, f);
 fwrite(&n->m1, sizeof(int), 1, f);
 fwrite(&n->d1, sizeof(int), 1, f);
 fwrite(&n->n2, sizeof(int), 1, f);
 fwrite(&n->p2, sizeof(int), 1, f);
 fwrite(&n->m2, sizeof(int), 1, f);
 fwrite(&n->d2, sizeof(int), 1, f);
 for (i=0;i<=n->m1;i++)
   {
    tmp = n->knot1[i];
    fwrite(&tmp, sizeof(sREAL), 1, f);
   }
 for (i=0;i<=n->m2;i++)
   {
    tmp = n->knot2[i];
    fwrite(&tmp, sizeof(sREAL), 1, f);
   }
 for (i=0;i<n->n1;i++)
 for (j=0;j<n->n2;j++)
   {
    tmp = n->w[i][j];
    fwrite(&tmp, sizeof(sREAL), 1, f);
   }
 for (i=0;i<n->n1;i++)
 for (j=0;j<n->n2;j++)
 for (k=0;k<S;k++)
   {
    tmp = n->D[i][j][k];
    fwrite(&tmp, sizeof(sREAL), 1, f);
   }
}

void save_binary_nurbs(FILE* f)
{
 int i;
 for (i=0;i<nNURBS;i++) save_binary_nurb(f, &g_nurbs[i]);
}

void load_binary_nurbs(FILE* f)
{
 int i;
 printf("Loading pretriangulated binary NURBS surfaces...\n");
 g_nurbs = (NURBS*)malloc(nNURBS*sizeof(NURBS));
 for (i=0;i<nNURBS;i++) load_binary_nurb(f, &g_nurbs[i]);
}


void load_binary_scene(FILE* f, Triangle** ts)
{
 int tlen,ti,i;
 int sx,sy;
 char texDir[1024];
 int* t_use;
 sREAL tmp;
 Triangle* t;
/* int nTex;*/
 fread(texDir, 4, 1, f);
 strcpy(texDir, "");
 fread(&sx,sizeof(int),1,f);
 fread(&sy,sizeof(int),1,f);
 if (sx <= 0 || sy <= 0) spanic("load_binary_scene: bad screen size", HERE);
 if (ovr_x > 0) sx = ovr_x;
 if (ovr_y > 0) sy = ovr_y;
 if (double_res) { sx *= 2; sy *= 2; }
 init_screen(&screen, sx, sy);
 fread(&tmp,sizeof(sREAL),1,f);
 observer.x = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 observer.y = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 observer.z = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 light.x = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 light.y = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 light.z = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 lookz = tmp;
 fread(&i,sizeof(int),1,f);
 if (!ovr_b) bkup = i;
 fread(&i,sizeof(int),1,f);
 if (!ovr_r) max_rec = i;
 fread(&tmp,sizeof(sREAL),1,f);
 if (!ovr_m) minshadow = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 if (!ovr_ms) maxshadow = tmp;
 fread(&tmp,sizeof(sREAL),1,f);
 if (!ovr_a) ambient = tmp;
 fread(&vlight,sizeof(int),1,f);
 fread(&i,sizeof(int),1,f);
 free_matrix(world, 4);
 free_matrix(worldn, 4);
 world = matrix(4);
 I_matrix(world, 4);
 worldn = matrix(4);
 I_matrix(worldn, 4);
 fread(&tlen,sizeof(int),1,f);
/* printf("tlen: %d\n", tlen);*/
 fread(texDir,tlen,1,f);
 texDir[tlen] = 0;
 fread(&nTex,sizeof(int),1,f);
 if (nTex < 0) spanic("load_binary_scene: negative textures count", HERE);
 create_textures(texDir);
 t_use = (int*)malloc(nTex*sizeof(int));
 for (i=0;i<nTex;i++) t_use[i] = 0;
 fread(&nTriangles,sizeof(int),1,f);
 if (nTriangles < 0) spanic("load_binary_scene: negative triangles count", HERE);
 if (nTriangles > 0) *ts = (Triangle*)malloc(nTriangles*sizeof(Triangle));
 t = *ts;
 fread(&nNURBS,sizeof(int),1,f);
 if (nNURBS > 0) load_binary_nurbs(f);
 if (nTriangles <= 0 && nNURBS <= 0) spanic("nothing to RT", HERE);
 for (i=0;i<nTriangles;i++)
   {
    fread(&ti,sizeof(int),1,f);
    if (ti < 0) spanic("load_binary_scene: bad triangle idx", HERE);
    t[i].idx = i;
    fread(&ti,sizeof(int),1,f);
    if (ti < -2) spanic("load_binary_scene: bad nurbs idx", HERE);
    t[i].nidx = ti;
    if (t[i].nidx >= 0) reconstruct_nurbs_addr(&t[i]);
    fread(&tmp, sizeof(sREAL),1,f);
    t[i].nca.x = tmp;
    fread(&tmp, sizeof(sREAL),1,f);
    t[i].nca.y = tmp;
    fread(&tmp, sizeof(sREAL),1,f);
    t[i].ncb.x = tmp;
    fread(&tmp, sizeof(sREAL),1,f);
    t[i].ncb.y = tmp;
    fread(&tmp, sizeof(sREAL),1,f);
    t[i].ncc.x = tmp;
    fread(&tmp, sizeof(sREAL),1,f);
    t[i].ncc.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].a.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].a.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].a.z = tmp;
    t[i].a.x += distorber1;
    t[i].a.y += distorber2;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].b.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].b.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].b.z = tmp;
    t[i].b.x += distorber1;
    t[i].b.y += distorber2;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].c.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].c.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f);
    t[i].c.z = tmp;
    t[i].c.x += distorber1;
    t[i].c.y += distorber2;
    fread(&tmp,sizeof(sREAL),1,f); t[i].ta.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].ta.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].tb.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].tb.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].tc.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].tc.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].na.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].na.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].na.z = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].nb.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].nb.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].nb.z = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].nc.x = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].nc.y = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].nc.z = tmp;
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
    fread(&tmp,sizeof(sREAL),1,f); t[i].mra.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mga.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mba.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mra.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mga.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mba.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mra.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mga.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mba.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mrb.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mgb.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mbb.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mrb.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mgb.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mbb.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mrb.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mgb.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mbb.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mrc.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mgc.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mbc.t = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mrc.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mgc.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mbc.s = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mrc.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mgc.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mbc.c = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].s.A = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].s.B = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].s.C = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].s.D = tmp;
    normalize_color_factors(t);
    if (t[i].s.A == 0. && t[i].s.B == 0. && t[i].s.C == 0. && t[i].s.D == 0.) compute_surface(&t[i]);
    fread(&tmp,sizeof(sREAL),1,f); t[i].mUR = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mDR = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mUG = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mDG = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mUB = tmp;
    fread(&tmp,sizeof(sREAL),1,f); t[i].mDB = tmp;
    if (t[i].mUR < .05 || t[i].mUR > 20.) spanic("load_binary_scene: bad mU/mD r value", HERE);
    if (t[i].mDR < .05 || t[i].mDR > 20.) spanic("load_binary_scene: bad mU/mD r value", HERE);
    if (t[i].mUG < .05 || t[i].mUG > 20.) spanic("load_binary_scene: bad mU/mD g value", HERE);
    if (t[i].mDG < .05 || t[i].mDG > 20.) spanic("load_binary_scene: bad mU/mD g value", HERE);
    if (t[i].mUB < .05 || t[i].mUB > 20.) spanic("load_binary_scene: bad mU/mD b value", HERE);
    if (t[i].mDB < .05 || t[i].mDB > 20.) spanic("load_binary_scene: bad mU/mD b value", HERE);
    fread(&tmp,sizeof(sREAL),1,f); t[i].n_dist = tmp;
    if (t[i].n_dist < 0.) spanic("load_binary_scene: negative n_dist", HERE);
    fread(&tmp,sizeof(sREAL),1,f); t[i].ca = tmp;
    if (t[i].ca < -1.) spanic("load_binary_scene: negative specularFact", HERE);
    fread(&t[i].faces,sizeof(int),1,f);
    if (t[i].faces != 1 && t[i].faces != 2) spanic("load_binary_scene: bad faces count (should be 1 or 2)", HERE);
    fread(&ti,sizeof(int),1,f);
    if (ti > nTex) spanic("load_binary_scene: texture out of range", HERE);
    if (ti > 0 && ti <= nTex)
      {
	/*t[i].tex = &texture[ti-1];*/
	t[i].tid = ti;
        t_use[ti-1] ++;
      }
    else t[i].tex = NULL;
    t[i].tid = ti;
/*    if (trans_used) transform_triangle(&t[i], world, worldn);*/
   }
 strcpy(textureDir, texDir);
 fclose(f);
 postprocess_textures(t_use, t, texDir);
 free(t_use);
}


void save_binary_scene(char* name, Triangle* t)
{
 int tlen,i,tu;
 char out[1024];
 char sig[5];
 FILE* f;
 sREAL tmp;
 tu = 0;
 strcpy(out,name);
 for (i=0;i<(int)strlen(out);i++) if (out[i] == '.') out[i] = 0;
 strcat(out,".bin");
 printf("Saving binary scene to: %s\n", out);
 f = fopen(out, "wb");
 if (!f) spanic("save_binary_scene: file not open for writing", HERE);
/* int nTex;*/
 strcpy(sig, "BINT");
 fwrite(sig,4,1,f);
 fwrite(&screen.x,sizeof(int),1,f);
 fwrite(&screen.y,sizeof(int),1,f);
 tmp = observer.x;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = observer.y;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = observer.z;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = light.x;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = light.y;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = light.z;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = lookz;
 fwrite(&tmp,sizeof(sREAL),1,f);
 fwrite(&bkup,sizeof(int),1,f);
 fwrite(&max_rec,sizeof(int),1,f);
 tmp = minshadow;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = maxshadow;
 fwrite(&tmp,sizeof(sREAL),1,f);
 tmp = ambient;
 fwrite(&tmp,sizeof(sREAL),1,f);
 fwrite(&vlight,sizeof(int),1,f);
 fwrite(&tu,sizeof(int),1,f);
 tlen = strlen(textureDir);
/* printf("textureDir: %s, len = %d\n", textureDir, tlen);*/
 fwrite(&tlen,sizeof(int),1,f);
 fwrite(textureDir,tlen,1,f);
 fwrite(&nTex,sizeof(int),1,f);
 fwrite(&nTriangles,sizeof(int),1,f);
 fwrite(&nNURBS,sizeof(int),1,f);
 if (nNURBS > 0) save_binary_nurbs(f);
 for (i=0;i<nTriangles;i++)
   {
    fwrite(&t[i].idx,sizeof(int),1,f);
    fwrite(&t[i].nidx,sizeof(int),1,f);
    tmp = t[i].nca.x;
    fwrite(&tmp, sizeof(sREAL),1,f);
    tmp = t[i].nca.y;
    fwrite(&tmp, sizeof(sREAL),1,f);
    tmp = t[i].ncb.x;
    fwrite(&tmp, sizeof(sREAL),1,f);
    tmp = t[i].ncb.y;
    fwrite(&tmp, sizeof(sREAL),1,f);
    tmp = t[i].ncc.x;
    fwrite(&tmp, sizeof(sREAL),1,f);
    tmp = t[i].ncc.y;
    fwrite(&tmp, sizeof(sREAL),1,f);
    t[i].a.x -= distorber1;
    t[i].a.y -= distorber2;
    tmp = t[i].a.x;
    fwrite(&tmp,sizeof(sREAL),1,f);
    tmp = t[i].a.y;
    fwrite(&tmp,sizeof(sREAL),1,f);
    tmp = t[i].a.z;
    fwrite(&tmp,sizeof(sREAL),1,f);
    t[i].a.x += distorber1;
    t[i].a.y += distorber2;
    t[i].b.x -= distorber1;
    t[i].b.y -= distorber2;
    tmp = t[i].b.x;
    fwrite(&tmp,sizeof(sREAL),1,f);
    tmp = t[i].b.y;
    fwrite(&tmp,sizeof(sREAL),1,f);
    tmp = t[i].b.z;
    fwrite(&tmp,sizeof(sREAL),1,f);
    t[i].b.x += distorber1;
    t[i].b.y += distorber2;
    t[i].c.x -= distorber1;
    t[i].c.y -= distorber2;
    tmp = t[i].c.x;
    fwrite(&tmp,sizeof(sREAL),1,f);
    tmp = t[i].c.y;
    fwrite(&tmp,sizeof(sREAL),1,f);
    tmp = t[i].c.z;
    fwrite(&tmp,sizeof(sREAL),1,f);
    t[i].c.x += distorber1;
    t[i].c.y += distorber2;
     tmp = t[i].ta.x ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].ta.y ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].tb.x ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].tb.y ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].tc.x ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].tc.y ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].na.x ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].na.y ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].na.z ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].nb.x ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].nb.y ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].nb.z ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].nc.x ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].nc.y ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].nc.z ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mra.t;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mga.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mba.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mra.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mga.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mba.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mra.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mga.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mba.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mrb.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mgb.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mbb.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mrb.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mgb.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mbb.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mrb.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mgb.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mbb.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mrc.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mgc.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mbc.t ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mrc.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mgc.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mbc.s ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mrc.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mgc.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mbc.c ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].s.A ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].s.B ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].s.C ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].s.D ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mUR ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mDR ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mUG ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mDG ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mUB ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].mDB ;    fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].n_dist;  fwrite(&tmp,sizeof(sREAL),1,f);
     tmp = t[i].ca ;     fwrite(&tmp,sizeof(sREAL),1,f);
     fwrite(&t[i].faces,sizeof(int),1,f);
     fwrite(&t[i].tid,sizeof(int),1,f);
   }
 fclose(f);
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


int is_blank(int zn)
{
 if (zn == ' ' || zn == '\t' || zn == '\n' || zn == '\r') return 1;
 else return 0;
}


void skip_to_eol(FILE* f)
{
 while (fgetc(f) != '\n') ;
}


void skip_comments(FILE* f)
{
 long pos;
 int zn,comment;
 int flag;
 while (1)
 {
  while (is_blank(fgetc(f))) ;
  fseek(f, -1, SEEK_CUR);
hash_comm:
  while (is_blank(fgetc(f))) ;
  fseek(f, -1, SEEK_CUR);
  pos = ftell(f);
  zn = fgetc(f);
  if (zn == '#') { skip_to_eol(f); pos = ftell(f); goto hash_comm; }
/*  printf("zn = '%c'\n", zn);*/
  if (zn != '/') { fseek(f, pos, SEEK_SET); return; }
  zn = fgetc(f);
  if (zn != '*') { fseek(f, pos, SEEK_SET); return; }
  flag = '*';
  comment = 1;
  while (comment)
   {
/*       printf("in comment. flag = %c\n", flag);*/
    zn = fgetc(f);
    if (zn == EOF) return;
    if (is_blank(zn)) continue;
    if (zn == '#') skip_to_eol(f);
    if (zn == '*') flag = '/';
    else if (zn == '/' && flag == '/') comment = 0;
    else flag = '*';
   }
  while (is_blank(fgetc(f))) ;
  fseek(f, -1, SEEK_CUR);
 }
}


#ifndef VS
int debug(char* f, int l, char* fmt, ...)
{
#ifdef DEBUG
 va_list ap;
 int err;
 if (!dbg) return 1;
 err = fprintf(dbg, "%s, line: %d: ", f, l);
 va_start(ap,fmt);
 err = vfprintf(dbg,fmt,ap);
 va_end(ap);
 fflush(dbg);
 return err;
#else
 return 0;
#endif
}
int v_fscanf(FILE* f, char* fmt, ...)
{
 va_list ap;
 int err;
 va_start(ap,fmt);
 skip_comments(f);
 err = vfscanf(f,fmt,ap);
 skip_comments(f);
 va_end(ap);
 return err;
}
#endif


void read_transformation_long(FILE* f, REAL*** mm, REAL*** mn, int* tuse, REAL* nor_dist, int* texid)
{
 int done;
 int nr;
 int ipar;
 REAL par,par2,par3;
 REAL** mat;
 char str[1024];
 mfscanf(f,"{\n");
 done = 0;
 *nor_dist = -1;
 *texid = -1;
 while (!done)
       {
        nr = mfscanf(f,"%s", str);
	if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
        if (!strcmp(str,"}")) { mfscanf(f,"\n"); done = 1; }
	else if (!strcmp(str,"NormalDistorber:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf%%\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   if (par < 0.) spanic("read_transformation: negative normal_distorber", HERE);
	   *nor_dist = par / 100.;
	  }
	else if (!strcmp(str,"NewTexture:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %d\n", &ipar);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   if (ipar < 0.) spanic("read_transformation: negative texture_id", HERE);
	   if (ipar > nTex) spanic("read_transformation: bad new texture id", HERE);
	   *texid = ipar;
	  }
	else if (!strcmp(str,"RotateX:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatex(mm,par);
	   m_rotatex(mn,par);
	  }
	else if (!strcmp(str,"RotateY:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatey(mm,par);
	   m_rotatey(mn,par);
	  }
	else if (!strcmp(str,"RotateZ:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatez(mm,par);
	   m_rotatez(mn,par);
	  }
	else if (!strcmp(str,"RotateNX:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatex(mn,par);
	  }
	else if (!strcmp(str,"RotateNY:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatey(mn,par);
	  }
	else if (!strcmp(str,"RotateNZ:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatez(mn,par);
	  }
	else if (!strcmp(str,"RotateWX:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatex(mm,par);
	  }
	else if (!strcmp(str,"RotateWY:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatey(mm,par);
	  }
	else if (!strcmp(str,"RotateWZ:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," %Lf\n", &par);
	   if (nr != 1) spanic("read_transformation: cant read parameter", HERE);
	   m_rotatez(mm,par);
	  }
	else if (!strcmp(str,"Translate:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," (%Lf,%Lf,%Lf)\n", &par, &par2, &par3);
	   if (nr != 3) spanic("read_transformation: cant read parameters: 3", HERE);
	   m_translate(mm,par,par2,par3);
	  }
	else if (!strcmp(str,"Scale:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," (%Lf,%Lf,%Lf)\n", &par, &par2, &par3);
/*	   printf("nr = %d\n", nr);*/
	   if (nr != 3) spanic("read_transformation: cant read parameters: 3", HERE);
	   m_scale(mm,par,par2,par3);
	   m_scale(mn,par,par2,par3);
	  }
	else if (!strcmp(str,"ScaleN:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," (%Lf,%Lf,%Lf)\n", &par, &par2, &par3);
/*	   printf("nr = %d\n", nr);*/
	   if (nr != 3) spanic("read_transformation: cant read parameters: 3", HERE);
	   m_scale(mn,par,par2,par3);
	  }
	else if (!strcmp(str,"ScaleW:"))
	  {
	   *tuse = (*tuse) + 1;
	   nr = mfscanf(f," (%Lf,%Lf,%Lf)\n", &par, &par2, &par3);
/*	   printf("nr = %d\n", nr);*/
	   if (nr != 3) spanic("read_transformation: cant read parameters: 3", HERE);
	   m_scale(mm,par,par2,par3);
	  }
	else if (!strcmp(str,"Matrix:") || !strcmp(str,"MatrixW:"))
	  {
	   *tuse = (*tuse) + 1;
	   mat = matrix(4);
	   nr = mfscanf(f," [(%Lf,%Lf,%Lf,%Lf),(%Lf,%Lf,%Lf,%Lf),(%Lf,%Lf,%Lf,%Lf),(%Lf,%Lf,%Lf,%Lf)]\n",
	   &mat[0][0], &mat[0][1], &mat[0][2], &mat[0][3],
	   &mat[1][0], &mat[1][1], &mat[1][2], &mat[1][3],
	   &mat[2][0], &mat[2][1], &mat[2][2], &mat[2][3],
	   &mat[3][0], &mat[3][1], &mat[3][2], &mat[3][3]);
	   if (nr != 16) spanic("read_transformation: cant read parameters: 16", HERE);
	   m_ownmatrix(mm,mat);
	   free_matrix(mat, 4);
	  }
	else if (!strcmp(str,"MatrixN:"))
	  {
	   *tuse = (*tuse) + 1;
	   mat = matrix(4);
	   nr = mfscanf(f," [(%Lf,%Lf,%Lf,%Lf),(%Lf,%Lf,%Lf,%Lf),(%Lf,%Lf,%Lf,%Lf),(%Lf,%Lf,%Lf,%Lf)]\n",
	   &mat[0][0], &mat[0][1], &mat[0][2], &mat[0][3],
	   &mat[1][0], &mat[1][1], &mat[1][2], &mat[1][3],
	   &mat[2][0], &mat[2][1], &mat[2][2], &mat[2][3],
	   &mat[3][0], &mat[3][1], &mat[3][2], &mat[3][3]);
	   if (nr != 16) spanic("read_transformation: cant read parameteris: 16", HERE);
	   m_ownmatrix(mn,mat);
	   free_matrix(mat, 4);
	  }
	else if (!strcmp(str,"Identity:"))
	  {
	   *tuse = 0;
	   mfscanf(f,"\n");
	   I_matrix(*mm, 4);
	   I_matrix(*mn, 4);
	  }
	else if (!strcmp(str,"IdentityW:"))
	  {
	   *tuse = 0;
	   mfscanf(f,"\n");
	   I_matrix(*mm, 4);
	  }
	else if (!strcmp(str,"IdentityN:"))
	  {
	   *tuse = 0;
	   mfscanf(f,"\n");
	   I_matrix(*mn, 4);
	  }
	else if (!strcmp(str,"NegateN:"))
	  {
	   *tuse = (*tuse) + 1;
	   mfscanf(f,"\n");
	   m_scale(mn,-1.,-1.,-1.);
	  }
	else if (!strcmp(str,"Invert:"))
	  {
	   *tuse = (*tuse) + 1;
	   mfscanf(f,"\n");
	   mat = invert_matrix(*mm, 4);
	   free_matrix(*mm, 4);
	   *mm = mat;
	   mat = invert_matrix(*mn, 4);
	   free_matrix(*mn, 4);
	   *mn = mat;
	  }
	else if (!strcmp(str,"InvertW:"))
	  {
	   *tuse = (*tuse) + 1;
	   mfscanf(f,"\n");
	   mat = invert_matrix(*mm, 4);
	   free_matrix(*mm, 4);
	   *mm = mat;
	  }
	else if (!strcmp(str,"InvertN:"))
	  {
	   *tuse = (*tuse) + 1;
	   mfscanf(f,"\n");
	   mat = invert_matrix(*mn, 4);
	   free_matrix(*mn, 4);
	   *mn = mat;
	  }
        else { printf("%s\n", str); spanic("read_transformation: unknown transformation", HERE); }
       }
}

void ptrlist_add(PtrList** head, void* t);
int nearly_equal(REAL x, REAL y, REAL eps);

void read_nurbs(FILE* f, NURBS* n)
{
 int idx,nr;
 int uK, uN;
 int i,j,k;
 nr = mfscanf(f, "NURBS: %d\n", &idx);
 if (nr != 1 || idx < 0) spanic("cannot read nurbs header", HERE);
 mfscanf(f, "{\n");
 nr = mfscanf(f, "DimU: %d\n", &n->p1);
 if (nr != 1) spanic("cannot read DimU", HERE);
 nr = mfscanf(f, "DimV: %d\n", &n->p2);
 if (nr != 1) spanic("cannot read DimV", HERE);
 nr = mfscanf(f, "nptsU: %d\n", &n->n1);
 if (nr != 1) spanic("cannot read nptsU", HERE);
 nr = mfscanf(f, "nptsV: %d\n", &n->n2);
 if (nr != 1) spanic("cannot read nptsV", HERE);
 nr = mfscanf(f, "divU: %d\n", &n->d1);
 if (nr != 1) spanic("cannot read divU", HERE);
 nr = mfscanf(f, "divV: %d\n", &n->d2);
 if (nr != 1) spanic("cannot read divV", HERE);
 nr = mfscanf(f, "ITex: %d\n", &n->itex);
 if (nr != 1) spanic("cannot read ITex", HERE);
 nr = mfscanf(f, "Faces: %d\n", &n->faces);
 if (nr != 1) spanic("cannot read Faces", HERE);
 nr = mfscanf(f, "Tid: %d\n", &n->tid);
 if (nr != 1) spanic("cannot read Tid", HERE);
 if (n->n1 <= n->p1) spanic("not enough U points", HERE);
 if (n->n2 <= n->p2) spanic("not enough V points", HERE);
 if (n->n1 < 2 || n->n2 < 2) spanic("not enough points", HERE);
 if (n->p1 < 1 || n->p2 < 1) spanic("too low degree", HERE);
 if (n->p1 >= 10 || n->p2 >= 10) spanic("too high degree", HERE);
 if (n->n1 > 32 || n->n2 > 32) spanic("too high points number", HERE);
 nr = mfscanf(f, "userNodes: %d\n", &uN);
 if (nr != 1) spanic("cannot read uN", HERE);
 nr = mfscanf(f, "userKnots: %d\n", &uK);
 if (nr != 1) spanic("cannot read uK", HERE);
 n->t1 = (REAL*)malloc(n->n1*sizeof(REAL));
 n->t2 = (REAL*)malloc(n->n2*sizeof(REAL));
 if (uN)
   {
    mfscanf(f, "UNodes:");
    for (i=0;i<n->n1;i++) 
      {
       nr = mfscanf(f, " %Lf", &n->t1[i]);
       if (nr != 1) spanic("cannot read Unode", HERE);
      }
    mfscanf(f, "\n");
    mfscanf(f, "VNodes:");
    for (i=0;i<n->n2;i++) 
      {
       nr = mfscanf(f, " %Lf", &n->t2[i]);
       if (nr != 1) spanic("cannot read Vnode", HERE);
      }
    mfscanf(f, "\n");
   }
 else
   {
    for (i=0;i<n->n1;i++) n->t1[i] = (REAL)i / (REAL)(n->n1 - 1);
    for (i=0;i<n->n2;i++) n->t2[i] = (REAL)i / (REAL)(n->n2 - 1);
   }
 /*printf("NodesU:");
 for (i=0;i<n->n1;i++) printf(" %f", n->t1[i]);
 printf("\n");
 printf("NodesV:");
 for (i=0;i<n->n2;i++) printf(" %f", n->t2[i]);
 printf("\n");*/
 n->m1 = n->n1 + n->p1;
 n->m2 = n->n2 + n->p2;
 n->knot1 = (REAL*)malloc((n->m1+1)*sizeof(REAL));
 n->knot2 = (REAL*)malloc((n->m2+1)*sizeof(REAL));
 if (uK)
   {
    mfscanf(f, "UKnots:");
    for (i=0;i<=n->m1;i++) 
      {
       nr = mfscanf(f, " %Lf", &n->knot1[i]);
       if (nr != 1) spanic("cannot read Uknot", HERE);
      }
    mfscanf(f, "\n");
    mfscanf(f, "VKnots:");
    for (i=0;i<=n->m2;i++) 
      {
       nr = mfscanf(f, " %Lf", &n->knot2[i]);
       if (nr != 1) spanic("cannot read Vknot", HERE);
      }
    mfscanf(f, "\n");
   }
 else
   {
    for (i=0;i<=n->p1;i++) n->knot1[i] = 0.;
    for (i=0;i<=n->p2;i++) n->knot2[i] = 0.;
    for (j=1;j<n->n1-n->p1;j++)
      {
	n->knot1[j+n->p1] = 0.;
	for (i=j;i<j+n->p1;i++) n->knot1[j+n->p1] += n->t1[i];
	n->knot1[j+n->p1] /= (REAL)n->p1;
      }
    for (j=1;j<n->n2-n->p2;j++)
      {
	n->knot2[j+n->p2] = 0.;
	for (i=j;i<j+n->p2;i++) n->knot2[j+n->p2] += n->t2[i];
	n->knot2[j+n->p2] /= (REAL)n->p2;
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
 free(n->t1);
 free(n->t2);
 n->t1 = n->t2 = NULL;
 n->D = matrix3(n->n1, n->n2, S);
 n->w = matrix2(n->n1, n->n2);
 for (i=0;i<n->n1;i++)
   {
    nr = mfscanf(f, "Line=%d\n", &k);
/*    printf(" i = %d\n", i);*/
    if (nr != 1 || k != i) spanic("No line definition", HERE);
    for (j=0;j<n->n2;j++)
      {
       nr = mfscanf(f, 
       "{v=(%Lf,%Lf,%Lf),t=(%Lf,%Lf,%Lf),s=(%Lf,%Lf,%Lf),c=(%Lf,%Lf,%Lf),"
       "f=(%Lf,%Lf,%Lf),sf=%Lf,nd=%Lf%%,tc=(%Lf,%Lf),w=%Lf}\n",&n->D[i][j][0],&n->D[i][j][1],&n->D[i][j][2],
        &n->D[i][j][3],&n->D[i][j][4],&n->D[i][j][5],&n->D[i][j][6],&n->D[i][j][7],&n->D[i][j][8],
        &n->D[i][j][9],&n->D[i][j][10],&n->D[i][j][11],&n->D[i][j][12],&n->D[i][j][13],&n->D[i][j][14],
	&n->D[i][j][15],&n->D[i][j][16],&n->D[i][j][17],&n->D[i][j][18], &n->w[i][j]);
        if (nr != S+1) 
	  { 
	   printf("data scanf error at: (%d,%d)\n", i,j); 
	   spanic("error", HERE); 
	  }
      }
     mfscanf(f, "\n");
    }
 n->ntri = 2*n->d1*n->d2;
 mfscanf(f, "}\n");
}

REAL b0(REAL* knot, REAL t, int i)
{
/* printf("checking: [%d,%d] %f in [%f,%f]\n", i,i+1,t,knot[i],knot[i+1]);*/
 if (t >= knot[i] && t <= knot[i+1]) return 1.;
 else return 0;
}

REAL*** precompute_nurbs(REAL* knot, REAL t, int dim, int npts)
{
 int i,crd,idx;
 REAL*** prep;
 REAL up;
 REAL down;
 REAL f1, f2;
/* printf("t = %f, dim = %d, idx = %d\n", t, dim, idx);*/
 prep = (REAL***)malloc(npts*sizeof(REAL**));
 for (idx=0;idx<npts;idx++)
 {
  prep[idx] = (REAL**)malloc((dim+1)*sizeof(REAL*));
  for (i=0;i<=dim;i++) prep[idx][i] = (REAL*)malloc(((dim-i)+1)*sizeof(REAL));
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

void prenurbs_free(REAL**** ptr, int dim, int npts)
{
 int i,j;
 for (i=0;i<npts;i++)
   for (j=0;j<=dim;j++) free((*ptr)[i][j]);
 for (i=0;i<npts;i++) free((*ptr)[i]);
 free(*ptr);
 *ptr = NULL;
}

REAL fastnurbs(NURBS* nurb, REAL u, REAL v, int xyz)
{
 register int i,j;
 REAL up,down;
 REAL factor;
 REAL b1,b2;
 REAL*** pre1;
 REAL*** pre2;
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

void fastnurbs_array(NURBS* nurb, REAL u, REAL v, REAL* tab, int siz)
{
 register int i,j,k;
 REAL up,down;
 REAL factor;
 REAL b1,b2;
 REAL*** pre1;
 REAL*** pre2;
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

void calc_normal(NURBS* nu, REAL u, REAL v, REAL* x, REAL* y, REAL* z)
{
 REAL x1,x2,x3;
 REAL y1,y2,y3;
 REAL z1,z2,z3;
 REAL nx,ny,nz;
 REAL dx1,dx2;
 REAL dy1,dy2;
 REAL dz1,dz2;
 REAL epsu;
 REAL epsv;
 REAL len;
/* REAL** pre;*/
 
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
/* printf("%Lf,%Lf,%Lf\n", dx1,dy1,dz1);*/
 
 nx = dy1*dz2 - dz1*dy2;
 ny = dz1*dx2 - dx1*dz2;
 nz = dx1*dy2 - dy1*dx2;
 
 len = sqrt(nx*nx + ny*ny + nz*nz);
/* printf("len = %Lf\n", len);*/
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

void triangulate_nurbs(NURBS* n, Triangle* tts)
{
 double u,U,v,V; 
 int i,j, k;
 n->values = matrix3(n->d1+1, n->d2+1, S);
 n->normals = matrix3(n->d1+1, n->d2+1, 3);
 for (i=0;i<=n->d1;i++)
   {
    u = (double)i/(double)n->d1;
    for (j=0;j<=n->d2;j++)
     {
      v = (double)j/(double)n->d2;
      calc_normal(n, u, v, &n->normals[i][j][0], &n->normals[i][j][1], &n->normals[i][j][2]);
      fastnurbs_array(n, u, v, n->values[i][j], S);
     }
   }
 k = n->idx;
 for (i=0;i<n->d1;i++)
  {
   u = (double)i/(double)n->d1;
   U = (double)(i+1)/(double)n->d1;
   for (j=0;j<n->d2;j++)
    {
/*	    printf("Filling triangle: %d\n", k);*/
     v = (double)j/(double)n->d2;
     V = (double)(j+1)/(double)n->d2;
     tts[k].nca.x = u;
     tts[k].nca.y = v;
     tts[k].ncb.x = U;
     tts[k].ncb.y = V;
     tts[k].ncc.x = U;
     tts[k].ncc.y = v;
     tts[k].nurbs = n;
     tts[k].nidx  = n->idx;
     tts[k].a.x = n->values[i][j][COORD_LX];
     tts[k].a.y = n->values[i][j][COORD_LY];
     tts[k].a.z = n->values[i][j][COORD_LZ];
     tts[k].b.x = n->values[i+1][j+1][COORD_LX];
     tts[k].b.y = n->values[i+1][j+1][COORD_LY];
     tts[k].b.z = n->values[i+1][j+1][COORD_LZ];
     tts[k].c.x = n->values[i+1][j][COORD_LX];
     tts[k].c.y = n->values[i+1][j][COORD_LY];
     tts[k].c.z = n->values[i+1][j][COORD_LZ];
     if (!n->itex)
       {
/*	   printf("no tex interpol.\n");*/
        tts[k].ta.x = u;
        tts[k].ta.y = v;
        tts[k].tb.x = U;
        tts[k].tb.y = V;
        tts[k].tc.x = U;
        tts[k].tc.y = v;
       }
     else
       {
        tts[k].ta.x = n->values[i][j][COORD_TX];
        tts[k].ta.y = n->values[i][j][COORD_TY];
        tts[k].tb.x = n->values[i+1][j+1][COORD_TX];
        tts[k].tb.y = n->values[i+1][j+1][COORD_TY];
        tts[k].tc.x = n->values[i+1][j][COORD_TX];
        tts[k].tc.y = n->values[i+1][j][COORD_TY];
       }
     tts[k].na.x = n->normals[i][j][0];
     tts[k].na.y = n->normals[i][j][1];
     tts[k].na.z = n->normals[i][j][2];
     tts[k].nb.x = n->normals[i+1][j+1][0];
     tts[k].nb.y = n->normals[i+1][j+1][1];
     tts[k].nb.z = n->normals[i+1][j+1][2];
     tts[k].nc.x = n->normals[i+1][j][0];
     tts[k].nc.y = n->normals[i+1][j][1];
     tts[k].nc.z = n->normals[i+1][j][2];
     tts[k].mra.t = n->values[i][j][COORD_TR];
     tts[k].mga.t = n->values[i][j][COORD_TG];
     tts[k].mba.t = n->values[i][j][COORD_TB];
     tts[k].mra.s = n->values[i][j][COORD_SR];
     tts[k].mga.s = n->values[i][j][COORD_SG];
     tts[k].mba.s = n->values[i][j][COORD_SB];
     tts[k].mra.c = n->values[i][j][COORD_CR];
     tts[k].mga.c = n->values[i][j][COORD_CG];
     tts[k].mba.c = n->values[i][j][COORD_CB];
     tts[k].mrb.t = n->values[i+1][j+1][COORD_TR];
     tts[k].mgb.t = n->values[i+1][j+1][COORD_TG];
     tts[k].mbb.t = n->values[i+1][j+1][COORD_TB];
     tts[k].mrb.s = n->values[i+1][j+1][COORD_SR];
     tts[k].mgb.s = n->values[i+1][j+1][COORD_SG];
     tts[k].mbb.s = n->values[i+1][j+1][COORD_SB];
     tts[k].mrb.c = n->values[i+1][j+1][COORD_CR];
     tts[k].mgb.c = n->values[i+1][j+1][COORD_CG];
     tts[k].mbb.c = n->values[i+1][j+1][COORD_CB];
     tts[k].mrc.t = n->values[i+1][j][COORD_TR];
     tts[k].mgc.t = n->values[i+1][j][COORD_TG];
     tts[k].mbc.t = n->values[i+1][j][COORD_TB];
     tts[k].mrc.s = n->values[i+1][j][COORD_SR];
     tts[k].mgc.s = n->values[i+1][j][COORD_SG];
     tts[k].mbc.s = n->values[i+1][j][COORD_SB];
     tts[k].mrc.c = n->values[i+1][j][COORD_CR];
     tts[k].mgc.c = n->values[i+1][j][COORD_CG];
     tts[k].mbc.c = n->values[i+1][j][COORD_CB];
     tts[k].mUR = 1.;
     tts[k].mUG = 1.;
     tts[k].mUB = 1.;
     tts[k].mDR = n->values[i][j][COORD_FR];
     tts[k].mDG = n->values[i][j][COORD_FG];
     tts[k].mDB = n->values[i][j][COORD_FB];
     tts[k].n_dist = n->values[i][j][COORD_ND] / 100.;
     tts[k].ca     = n->values[i][j][COORD_SF];
     tts[k].faces  = n->faces;
     tts[k].tid    = n->tid;
     k++;
/*	    printf("Filling triangle: %d\n", k);*/
     tts[k].nca.x = u;
     tts[k].nca.y = v;
     tts[k].ncb.x = u;
     tts[k].ncb.y = V;
     tts[k].ncc.x = U;
     tts[k].ncc.y = V;
     tts[k].nurbs = n;
     tts[k].nidx  = n->idx;
     tts[k].a.x = n->values[i][j][COORD_LX];
     tts[k].a.y = n->values[i][j][COORD_LY];
     tts[k].a.z = n->values[i][j][COORD_LZ];
     tts[k].b.x = n->values[i][j+1][COORD_LX];
     tts[k].b.y = n->values[i][j+1][COORD_LY];
     tts[k].b.z = n->values[i][j+1][COORD_LZ];
     tts[k].c.x = n->values[i+1][j+1][COORD_LX];
     tts[k].c.y = n->values[i+1][j+1][COORD_LY];
     tts[k].c.z = n->values[i+1][j+1][COORD_LZ];
     if (!n->itex)
       {
        tts[k].ta.x = u;
        tts[k].ta.y = v;
        tts[k].tb.x = u;
        tts[k].tb.y = V;
        tts[k].tc.x = U;
        tts[k].tc.y = V;
       }
     else
       {
        tts[k].ta.x = n->values[i][j][COORD_TX];
        tts[k].ta.y = n->values[i][j][COORD_TY];
        tts[k].tb.x = n->values[i][j+1][COORD_TX];
        tts[k].tb.y = n->values[i][j+1][COORD_TY];
        tts[k].tc.x = n->values[i+1][j+1][COORD_TX];
        tts[k].tc.y = n->values[i+1][j+1][COORD_TY];
       }
     tts[k].na.x = n->normals[i][j][0];
     tts[k].na.y = n->normals[i][j][1];
     tts[k].na.z = n->normals[i][j][2];
     tts[k].nb.x = n->normals[i][j+1][0];
     tts[k].nb.y = n->normals[i][j+1][1];
     tts[k].nb.z = n->normals[i][j+1][2];
     tts[k].nc.x = n->normals[i+1][j+1][0];
     tts[k].nc.y = n->normals[i+1][j+1][1];
     tts[k].nc.z = n->normals[i+1][j+1][2];
     tts[k].mra.t = n->values[i][j][COORD_TR];
     tts[k].mga.t = n->values[i][j][COORD_TG];
     tts[k].mba.t = n->values[i][j][COORD_TB];
     tts[k].mra.s = n->values[i][j][COORD_SR];
     tts[k].mga.s = n->values[i][j][COORD_SG];
     tts[k].mba.s = n->values[i][j][COORD_SB];
     tts[k].mra.c = n->values[i][j][COORD_CR];
     tts[k].mga.c = n->values[i][j][COORD_CG];
     tts[k].mba.c = n->values[i][j][COORD_CB];
     tts[k].mrb.t = n->values[i][j+1][COORD_TR];
     tts[k].mgb.t = n->values[i][j+1][COORD_TG];
     tts[k].mbb.t = n->values[i][j+1][COORD_TB];
     tts[k].mrb.s = n->values[i][j+1][COORD_SR];
     tts[k].mgb.s = n->values[i][j+1][COORD_SG];
     tts[k].mbb.s = n->values[i][j+1][COORD_SB];
     tts[k].mrb.c = n->values[i][j+1][COORD_CR];
     tts[k].mgb.c = n->values[i][j+1][COORD_CG];
     tts[k].mbb.c = n->values[i][j+1][COORD_CB];
     tts[k].mrc.t = n->values[i+1][j+1][COORD_TR];
     tts[k].mgc.t = n->values[i+1][j+1][COORD_TG];
     tts[k].mbc.t = n->values[i+1][j+1][COORD_TB];
     tts[k].mrc.s = n->values[i+1][j+1][COORD_SR];
     tts[k].mgc.s = n->values[i+1][j+1][COORD_SG];
     tts[k].mbc.s = n->values[i+1][j+1][COORD_SB];
     tts[k].mrc.c = n->values[i+1][j+1][COORD_CR];
     tts[k].mgc.c = n->values[i+1][j+1][COORD_CG];
     tts[k].mbc.c = n->values[i+1][j+1][COORD_CB];
     tts[k].mUR = 1.;
     tts[k].mUG = 1.;
     tts[k].mUB = 1.;
     tts[k].mDR = n->values[i][j][COORD_FR];
     tts[k].mDG = n->values[i][j][COORD_FG];
     tts[k].mDB = n->values[i][j][COORD_FB];
     tts[k].n_dist = n->values[i][j][COORD_ND] / 100.;
     tts[k].ca     = n->values[i][j][COORD_SF];
     tts[k].faces  = n->faces;
     tts[k].tid    = n->tid;
     k++;
    }
  }
 for (i=0;i<=n->d1;i++)
 for (j=0;j<=n->d2;j++)
   {
    free(n->values[i][j]);
    free(n->normals[i][j]);
   }
 for (i=0;i<=n->d1;i++)
   {
    free(n->values[i]);
    free(n->normals[i]);
   }
 free(n->values);
 free(n->normals);
 n->values = n->normals = NULL;
}

int read_nurbs_from_file(FILE* f, int n, int nt)
{
 int i;
 int ntris;
 nNURBS = n;
 g_nurbs = (NURBS*)malloc(nNURBS*sizeof(NURBS));
 for (i=0;i<n;i++) read_nurbs(f, &g_nurbs[i]);
 ntris = 0;
 for (i=0;i<n;i++) 
   {
    printf("NURBS surface: %d, uses triangles: [%d,%d]\n"
	, i, ntris + nt, ntris + g_nurbs[i].ntri + nt);
    g_nurbs[i].idx = ntris + nt;
    ntris += g_nurbs[i].ntri;
   }
 return ntris;
}

void load_scene(Triangle** ts, char* scenefile)
{
 FILE* f;
 char texDir[1024];
 char str[1024];
 REAL dd;
 int texid;
 int nr,sx,sy,n,i,i1,i2;
 long pos;
 int gtrans;
 REAL **lm, **lmn;
 int *t_use;
 char** tmap;
 ListTransform* ltemp;
 ltemp = NULL;
 tlist = NULL;
 tmap = NULL;
 trans_used = 0;
 vlight = 0;
 lookz = .333333;
 world = matrix(4);
 I_matrix(world, 4);
 worldn = matrix(4);
 I_matrix(worldn, 4);
 f = fopen(scenefile,"rb");	/* is ok with r? */
 if (!f) spanic("load_scene: cant open file", HERE);
 if (is_binary(f)) { load_binary_scene(f, ts); return; }
 nr = mfscanf(f, "Screen: (%d,%d)\n", &sx, &sy);
 if (nr != 2) spanic("load_scene: cant read screen def", HERE);
 if (sx <= 0 || sy <= 0) spanic("load_scene: bad screen def", HERE);
 if (ovr_x > 0) sx = ovr_x;
 if (ovr_y > 0) sy = ovr_y;
 if (double_res) { sx *= 2; sy *= 2; }
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr != 1) spanic("load_scene: cant read global norm or other option", HERE);
 if (!strcmp(str, "NormalDistorber:"))
   {
    nr = mfscanf(f, " %Lf%%\n", &dd);
    if (nr != 1) spanic("load_scene: cant read normal distorber definition", HERE);
    if (dd < 0.) dd = 0.;
    if (!ovr_no) global_dist = dd/100.;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr != 1) spanic("load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "MinShadow:"))
   {
    nr = mfscanf(f, " %Lf\n", &dd);
    if (nr != 1) spanic("load_scene: cant read minshadow", HERE);
    if (dd < 0.) dd = 0.;
    if (dd > 1.) dd = 1.;
    if (!ovr_m) minshadow = 1. - dd;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr != 1) spanic("load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "MaxShadow:"))
   {
    nr = mfscanf(f, " %Lf\n", &dd);
    if (nr != 1) spanic("load_scene: cant read maxshadow", HERE);
    if (dd < 0.) dd = 0.;
    if (dd > 1.) dd = 1.;
    if (!ovr_ms) maxshadow = 1. - dd;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr != 1) spanic("load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "Ambient:"))
   {
    nr = mfscanf(f, " %Lf\n", &dd);
    if (nr != 1) spanic("load_scene: cant read ambient", HERE);
    if (dd < 0.) dd = 0.;
    if (dd > 1.) dd = 1.;
    if (!ovr_a) ambient = dd;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr != 1) spanic("load_scene: cant read shadow or recurse defs", HERE);
 if (!strcmp(str, "MaxRecurse:"))
   {
    nr = mfscanf(f, " %d\n", &i);
    if (nr != 1) spanic("load_scene: cant read maxrecurse", HERE);
    if (i < 0) i = 0;
    if (!ovr_r) max_rec = i;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f, "%s", str);
 if (nr != 1) spanic("load_scene: cant read backup or observer def", HERE);
 if (!strcmp(str, "Backup:"))
   {
    nr = mfscanf(f, " %d\n", &i);
    if (nr != 1) spanic("load_scene: cant read backup definition", HERE);
    if (i < 0.) i = 0.;
    if (!ovr_b) bkup = i;
   }
 else fseek(f, pos, SEEK_SET);
 init_screen(&screen, sx, sy);
 mfscanf(f, "Observer: ");
 if (!read_vertex(f, &observer)) spanic("load_scene: cant read observer def", HERE);
 mfscanf(f,"\n");
 pos = ftell(f);
 nr = mfscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "ObserverTransform:"))
       {
/*	   printf("Observer transform.\n");*/
	lm = matrix(4);
	lmn = matrix(4);
	I_matrix(lm, 4);
	I_matrix(lmn, 4);
 	read_transformation(f, &lm, &lmn, &i, &dd);
	transform_observer(&observer, lm, lmn);
	free_matrix(lm, 4);
	free_matrix(lmn, 4);
	lm = lmn = NULL;
       }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "LookZ:"))
   {
    nr = mfscanf(f, "%Lf%%\n", &lookz);
    if (nr != 1) spanic("load_scene: cant read lookZ def", HERE);
    lookz /= 100.;
   }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "Light:"))
   {
    if (!read_vertex(f, &light)) spanic("load_scene: cant read light vertex def", HERE);
   }
 else if (nr == 1 && !strcmp(str, "VLight:"))
   {
    if (!read_vector(f, &light)) spanic("load_scene: cant read light vector def", HERE);
    vlight = 1;
    if (nearly_equal(length(&light), 0., 1e-6)) spanic("load_scene: vector light: length too close 0\n", HERE);
    normalize(&light);
   }
 else spanic("load_scene: cant read light def", HERE);
 pos = ftell(f);
 nr = mfscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "LightTransform:"))
       {
	lm = matrix(4);
	lmn = matrix(4);
	I_matrix(lm, 4);
	I_matrix(lmn, 4);
 	read_transformation(f, &lm, &lmn, &i, &dd);
 	/*print_matrix(lm, 4);
	print_matrix(lmn, 4);
	printf("done.\n");*/
	transform_light(&light, lm, lmn, vlight);
/* printf("(%Lf,%Lf,%Lf)\n", light.x, light.y, light.z);*/
	free_matrix(lm, 4);
	free_matrix(lmn, 4);
	lm = lmn = NULL;
       }
 else fseek(f, pos, SEEK_SET);
 pos = ftell(f);
 nr = mfscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "WorldTransform:"))
       {
	read_transformation(f, &world, &worldn, &trans_used, &dd);
        if (dd > 0. && !ovr_no && global_dist <= 0.)
	  {
	   global_dist = dd;
	   printf("world_transform -> set global_dist to: %Lf\n", global_dist);
	  }
       }
 else fseek(f, pos, SEEK_SET);
 nr = mfscanf(f, "TexDirectory: %s\n", texDir);
 if (nr != 1) spanic("load_scene: cant read texture directory", HERE);
 nr = mfscanf(f, "NumTextures: %d\n", &nTex);
 if (nr != 1) spanic("load_scene: cant read num textures", HERE);
 if (nTex < 0) spanic("load_scene: negative textures count", HERE);
 create_textures(texDir);
 if (texture)
  {
   t_use = (int*)malloc(nTex*sizeof(int));
   for (i=0;i<nTex;i++) t_use[i] = 0;
  }
 else t_use = NULL;
 pos = ftell(f);
 nr = mfscanf(f,"%s\n", str);
 if (nr == 1 && !strcmp(str, "TextureMapping:"))
   {
    mfscanf(f,"{\n");
    tmap = (char**)malloc(nTex*sizeof(char*));
    for (i=0;i<nTex;i++) tmap[i] = NULL;
    while (1)
      {
       pos = ftell(f);
       nr = mfscanf(f, "%s\n", str);
       if (nr == 1 && !strcmp(str,"}")) break;
       else fseek(f, pos, SEEK_SET);
       nr = mfscanf(f, "%d:", &i);
       if (nr != 1) spanic("read_triangle: cant read mapping idx", HERE);
       if (i < 1 || i > nTex) spanic("read_triangle: bad mapping idx", HERE);
       i--;
       if (tmap[i]) spanic("read_trianmgle: mapping index already used", HERE);
       tmap[i] = (char*)malloc(512*sizeof(char));
       nr = mfscanf(f,"%s\n", tmap[i]);
       if (nr != 1) spanic("load_scene: cant read texture mapping name\n", HERE);
      }
   }
 else fseek(f, pos, SEEK_SET);
 nr = mfscanf(f, "nTriangles: %d\n", &n);
 if (nr != 1) spanic("load_scene: cant read num triangles", HERE);
 if (n < 0) spanic("load_scene: negative triangles count", HERE);
 nNURBS = 0;
 pos = ftell(f);
 nr = mfscanf(f, "nNURBS: %d\n", &i);
 if (nr != 1) fseek(f, pos, SEEK_SET);
 else
   {
    if (i <= 0) spanic("load_scene: negative NURBS count", HERE);
/*    printf("READING NURBS SURFACE!\n");*/
    nTriNURBS = read_nurbs_from_file(f,i,n);
    printf("Adding %d triangles as pretriangulation of surface..\n", nTriNURBS);
    n += nTriNURBS;
   }
 if (n > 0) *ts = (Triangle*)malloc(n*sizeof(Triangle));
 if (n <= 0 && nNURBS <= 0) spanic("nothing to RT", HERE);
 if (nNURBS > 0)
   {
    for (i=0;i<nNURBS;i++) triangulate_nurbs(&g_nurbs[i], *ts);
   }
 gtrans = 1;
 while (gtrans)
   {
    pos = ftell(f);
    nr = mfscanf(f,"%s\n", str);
    if (nr == 1 && !strcmp(str, "ListTransform:"))
       {
	nr = mfscanf(f, "[%d,%d]", &i1, &i2);
	if (nr != 2) spanic("load_scene: bad list transformation (bad idcices format)", HERE);
	lm = matrix(4);
	I_matrix(lm, 4);
	lmn = matrix(4);
	I_matrix(lmn, 4);
	read_transformation_long(f, &lm, &lmn, &i, &dd, &texid);
	ltemp = (ListTransform*)malloc(sizeof(ListTransform));
	ltemp->M = lm;
	ltemp->MN = lmn;
	ltemp->i1 = i1;
	ltemp->i2 = i2;
	ltemp->t_id = texid;
	if (dd < 0.) ltemp->n_dist = -1.;
	else ltemp->n_dist = dd;
	if (texid < 0.) ltemp->t_id = -1;
	else ltemp->t_id = texid;
	ptrlist_add(&tlist, (void*)ltemp);
       }
   else { fseek(f, pos, SEEK_SET); gtrans = 0; }
   }
 nTriangles = n;
 for (i=0;i<n-nTriNURBS;i++)
   {
    if (!read_triangle(f, &((*ts)[i]), &i, t_use, tmap, *ts, 1)) spanic("load_scene: general triangle read failure", HERE);
/*    printf("(%Lf,%Lf,%Lf)\n", (*ts)[0].a.x, (*ts)[0].a.y, (*ts)[0].a.z);*/
/*    printf("read %d triangles\n", i);*/
   }
 for (i=n-nTriNURBS;i<n;i++) 
   {
    if (!read_triangle(f, &((*ts)[i]), &i, t_use, tmap, *ts, 0)) spanic("load_scene: general triangle finalize failure", HERE);
   }
 if (tmap)
   {
    for (i=0;i<nTex;i++) if (tmap[i]) free(tmap[i]);
    free(tmap);
    tmap = NULL;
   }
 fclose(f);
 postprocess_textures(t_use, *ts, texDir);
 strcpy(textureDir, texDir);
 if (want_bin) save_binary_scene(scenefile, *ts);
 if (t_use) free(t_use);
}


void free_scene(Triangle** ts)
{
 int i;
 free((void*)(*ts));
 for (i=0;i<nTex;i++) free_texture(&texture[i]);
 free(texture);
}


void set_color_uchar(unsigned char* s, int Y, int x, int y, int r, int g, int b)
{
 s[3*(Y * x + y)   ] = r;
 s[3*(Y * x + y) +1] = g;
 s[3*(Y * x + y) +2] = b;
}


void set_color(Screen* s, int x, int y, int r, int g, int b)
{
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    s->pixels[3*(s->y * x + y)   ] = r;
    s->pixels[3*(s->y * x + y) +1] = g;
    s->pixels[3*(s->y * x + y) +2] = b;
   }
 else spanic("set_color: idices out of range", HERE);
}


void color_background(Screen* s, int x, int y)
{
 set_color(s, x, y, BACK_R, BACK_G, BACK_B);
}


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


/*REAL compute_angle2(Vertex* a, Vertex* b, Vertex* c)
{
 REAL x1,x2,y1,y2,z1,z2,len;
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
 REAL kat;
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


void blist_add(BList** head, REAL ix, REAL iy, REAL iz, REAL ax, REAL ay, REAL az, Triangle* t1, Triangle* t2)
{
 BList* temp;
 if (*head == NULL)
   {
/*     nmalloc++;*/
    *head = (BList*)(malloc(sizeof(BList)));
    /* FIXME: test disabled for speed\n"); */
/*    if (!(*head)) spanic(":", HERE);*/
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
/*  if (!temp) spanic(":", HERE);*/
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


void ptrlist_add(PtrList** head, void* t)
{
 PtrList* temp;
 if (*head == NULL)
   {
/*     nmalloc++;*/
    *head = (PtrList*)(malloc(sizeof(PtrList)));
    /* FIXME: test disabled for speed\n"); */
/*    if (!(*head)) spanic(":", HERE);*/
    (*head)->ptr = (void*)t;
    (*head)->next = NULL;
    (*head)->prev = NULL;
    return;
   }
/*  nmalloc++;*/
  temp = (PtrList*)(malloc(sizeof(PtrList)));
/*  if (!temp) spanic(":", HERE);*/
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
/*    if (!(*head)) spanic(":", HERE);*/
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
/*  if (!temp) spanic(":", HERE);*/
  temp->P.x = v->x;
  temp->P.y = v->y;
  temp->P.z = v->z;
  temp->idx = idx;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
}

void get_nurbs_pt(Triangle* t, Vector* p);

int intersection_triangle(Triangle* tr, Ray* r, Vertex* ret)
{
 REAL t, tmp;
 Vertex rett;
 Surface* s;
 /* FIXME: some tests skipped for speed.\n"); */
/* if (!ret) spanic(":", HERE);*/
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
	if (tr->nurbs && n_pl >= 1) get_nurbs_pt(tr, ret);
/*	  printf("check the same: (%f,%f,%f)\n", (rett.x-r->P.x)*r->d.x,(rett.y-r->P.y)*r->d.y,(rett.z-r->P.z)*r->d.z);*/
/*	printf("passed!\n");*/
	return 1;
      }
    else return 0;
   }
 else return 0;
}

REAL distance(Vertex*, Vertex*);
REAL pseudo_distance(Vertex*, Vertex*);

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
/* if (!head) return;*/		/* FIXME: optimized */
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


void travel_tree(BTree* t)
{
 if (t->r) travel_tree(t->r);
 if (t->l) travel_tree(t->l);
 g_idx ++;
}


void assoc_idx_tree(BTree* t, void** pi)
{
 if (t->r) assoc_idx_tree(t->r, pi);
 if (t->l) assoc_idx_tree(t->l, pi);
 pi[g_idx] = t;
 g_idx ++;
}


void count_indices()
{
 BTree* tmp;
 g_idx = 0;
 tmp = btree;
 travel_tree(tmp);
/* printf("Indices in the tree: %d\n", g_idx);*/
}


void create_assoc_table(void** pi)
{
 BTree* tmp;
 g_idx = 0;
 tmp = btree;
 assoc_idx_tree(tmp, pi);
/* printf("Associative table created,\n");*/
}


void save_box_to_file(FILE* f, Box* b)
{
 fprintf(f, " Box:\n {\n");
 if (b->t1) fprintf(f, "  Triangle: %d\n", b->t1->idx);
 else fprintf(f, "  Triangle: -1\n");
 if (b->t2) fprintf(f, "  Triangle: %d\n", b->t2->idx);
 else fprintf(f, "  Triangle: -1\n");
 fprintf(f, "  (%Lf,%Lf,%Lf)\n  (%Lf,%Lf,%Lf)\n",
	 b->minx,b->miny,b->minz,b->maxx,b->maxy,b->maxz);
 fprintf(f, " }\n");
}


void save_bin_box_to_file(FILE* f, Box* b)
{
 int mi;
 sREAL tmp;
 mi = -1;
 if (b->t1) fwrite(&b->t1->idx, sizeof(int), 1, f);
 else fwrite(&mi, sizeof(int), 1, f);
 if (b->t2) fwrite(&b->t2->idx, sizeof(int), 1, f);
 else fwrite(&mi, sizeof(int), 1, f);
  tmp = b->minx; fwrite(&tmp, sizeof(sREAL), 1, f);
  tmp = b->miny; fwrite(&tmp, sizeof(sREAL), 1, f);
  tmp = b->minz; fwrite(&tmp, sizeof(sREAL), 1, f);
  tmp = b->maxx; fwrite(&tmp, sizeof(sREAL), 1, f);
  tmp = b->maxy; fwrite(&tmp, sizeof(sREAL), 1, f);
  tmp = b->maxz; fwrite(&tmp, sizeof(sREAL), 1, f);
}


void write_btree_to_file(BTree* t, FILE* f, void** p_idx, int n)
{
 int i;
 if (t->r) write_btree_to_file(t->r, f, p_idx, n);
 if (t->l) write_btree_to_file(t->l, f, p_idx, n);
 fprintf(f, "Btree: %d\n{\n", g_idx);
 if (t->r)
   {
    for (i=0;i<n;i++) if (p_idx[i] == (void*)t->r)
      {
       fprintf(f, " R: %d\n", i);
       goto r_done;
      }
    spanic("write_btree_to_file: bad index associative table", HERE);
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
    spanic("write_btree_to_file: bad index associative table", HERE);
   }
 else fprintf(f, " L: -1\n");
 l_done:
 save_box_to_file(f, t->b);
 fprintf(f, "}\n");
 g_idx ++;
}


void write_bin_btree_to_file(BTree* t, FILE* f, void** p_idx, int n)
{
 int i;
 int mi;
 if (t->r) write_bin_btree_to_file(t->r, f, p_idx, n);
 if (t->l) write_bin_btree_to_file(t->l, f, p_idx, n);
 mi = -1;
 fwrite(&g_idx, sizeof(int), 1, f);
 if (t->r)
   {
    for (i=0;i<n;i++) if (p_idx[i] == (void*)t->r)
      {
       fwrite(&i, sizeof(int), 1, f);
       goto r_done;
      }
    spanic("write_btree_to_file: bad index associative table", HERE);
   }
 else fwrite(&mi, sizeof(int), 1, f);
 r_done:
 if (t->l)
   {
    for (i=0;i<n;i++) if (p_idx[i] == (void*)t->l)
      {
       fwrite(&i, sizeof(int), 1, f);
       goto l_done;
      }
    spanic("write_btree_to_file: bad index associative table", HERE);
   }
 else fwrite(&mi, sizeof(int), 1, f);
 l_done:
 save_bin_box_to_file(f, t->b);
 g_idx ++;
}


void save_preprocessed()
{
 void **p_idx;
 int n,i;
 FILE* f;
 char fn[2048];
 int bin;
 BTree* tmp;
 if (!want_save_pps) return;
 printf("Saving preprocessed scene...\n");
 count_indices();
 n = g_idx;
 bin = want_bin_pps;
/* if (bin) printf("Using binary format.\n");*/
 p_idx = (void**)malloc(n*sizeof(void*));
 create_assoc_table(p_idx);
 strcpy(fn, scenef);
 if (!strstr(fn, ".")) strcat(fn, ".pps");
 else
   {
    i = 0;
    while (fn[i] != '.') i++;
    fn[i] = 0;
    strcat(fn, ".btree");
   }
 printf("Writing preprocessed scene to: %s\n", fn);
 f = fopen(fn, "w");
 if (!f) spanic("save_preprocessed: cannot write to file", HERE);
 tmp = btree;
 g_idx = 0;
 if (!bin) fprintf(f, "nElem: %d\n", n);
 else { fwrite("BINT", 4, 1, f); fwrite(&n, sizeof(int), 1, f); }
 for (i=0;i<n;i++) if (p_idx[i] == (void*)btree)
   {
    if (!bin) fprintf(f,"Root: %d\n", i);
    else fwrite(&i, sizeof(int), 1, f);
    break;
   }
 if (!bin) write_btree_to_file(tmp, f, p_idx, n);
 else write_bin_btree_to_file(tmp, f, p_idx, n);
 fclose(f);
}


void intersection_btree(VList** head, Ray* r, BTree* tree)
{
 Box* b;
 register REAL K,x,y,z;
 Vertex ret;
/* if ((tree->l && !tree->r) || (!tree->l && tree->r)) spanic(":", HERE);*/
 /* OPTIMIZE: miliardes times called when RT, must be SUPEROPTIMIZED */
 proc_bt += 1.0;
 if (!tree->l && !tree->r)
 {
  proc_tr += 2.;
  /* FIXME: code is maximal optimized and assumes that BTree is correct */
  /* so dont check if triangle1 exist because it must */
  /* if You mess with *.btree files be afraid */
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
          goto lcase;		/* GOTO USED FOR SPEED */
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
 REAL dist,ndist;
#ifdef DEBUG
 int ts;
#endif
 int sol;
 VList* head;
 VList* hd;
 VList* amb;
/* VList* tamb;*/
 head = NULL;
 glob_sol = 0;
/* proc_tr = 0;*/
 intersection_btree(&head, r, btree);
/* printf("%f %% (%d/%d)\n", (REAL)(proc_tr * 100)/(REAL)nTriangles, proc_tr, nTriangles);*/
 sol = glob_sol;
 proc_tr2 += (REAL)nTriangles;
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
#ifndef DEBUG
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
#else
 /* comment 3 lines above to have ambiguous solution checking */
 ts = 0;
 while (head)
   {
    if (nearly_equal(distance(&r->P, &head->P), dist, 1e-9))	/* FIXME: was 1e-9, the smallest ossible */
      {								/* to avoid ambiguous solutions */
       ts++;
       vlist_add(&amb, &head->P, head->idx);
       debug(HERE, "detected ambiguous solutions\n");
       /* FIXME: break here to avoid ambiguous solutions */
/*       break;*/
      }
    head = head->next;
   }
 if (ts > 1)
 {
/*  tamb = amb;*/
/*  printf("Ambiguous solutions.\n");*/
     debug(HERE, "*");
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
/*  spanic(":", HERE);*/
  return sol;
 }
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", *idx,ret->x, ret->y, ret->z);*/
/* printf("big_free\n");*/
 head = hd;
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
#endif
}


int intersection_old(Triangle* tr, Ray* r, Vertex* ret, int* idx)
{
 register int i;
 int sol;
 REAL dist,ndist;
#ifdef DEBUG
 int ts;
#endif
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
 proc_tr += (REAL)nTriangles;
 proc_tr2 += (REAL)nTriangles;
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
#ifndef DEBUG
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
#else
 /* comment 3 lines above to have ambiguous solution checking */
 ts = 0;
 while (head)
   {
    if (nearly_equal(pseudo_distance(&r->P, &head->P), dist, 1e-9))	/* FIXME: was 1e-9, the smallest ossible */
      {								/* to avoid ambiguous solutions */
       ts++;
       vlist_add(&amb, &head->P, head->idx);
       debug(HERE, "detected ambiguous solutions\n");
       /* FIXME: break here to avoid ambiguous solutions */
/*       break;*/
      }
    head = head->next;
   }
 if (ts > 1)
 {
/*  tamb = amb;*/
/*  printf("Ambiguous solutions.\n");*/
     debug(HERE, "*");
/*     fflush(stdout);*/
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
/*  spanic(":", HERE);*/
  return sol;
 }
/* printf("intersection with triangle: %d --> (%f,%f,%f)\n", *idx,ret->x, ret->y, ret->z);*/
/* printf("big_free\n");*/
 head = hd;
 vlist_free(&head);
 vlist_free(&amb);	/* FIXME: there were memory leaks */
 return sol;
#endif
}


int line_intersect(Vertex* p, Ray* l1, Ray* l2)
{
 REAL x1,y1,z1,x2,y2,z2;
 REAL dx1,dy1,dz1,dx2,dy2,dz2;
 register REAL K,L;
 REAL xr,yr,zr;
/* printf("LINE_INTERSECT:\n");*/
/* print_ray(l1);*/
/* print_ray(l2);*/
/* printf("LINE_INTERSECT_ENDS\n");*/
 /* FIXME: skip some check for speed */
/* if (!p || !l1 || !l2) spanic(":", HERE);*/
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
	 /* FIXME: OPTIMIZE */
#ifdef DEBUG
         if (!nearly_equal(x1+K*dx1, x2+L*dx2, 3e-4))
	 {
	     debug(HERE, "not3dx\n");
	     return 0;
	 }
#endif
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
#ifdef DEBUG
     if (!nearly_equal(y1+K*dy1, y2+L*dy2, 3e-4))
	 {
	     debug(HERE, "not3dy\n");
	     return 0;
	 }
#endif
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
#ifdef DEBUG
 if (!nearly_equal(z1+K*dz1, z2+L*dz2, 3e-4))
	 {
	     debug(HERE, "difference: %0.13lf\n", (z1+K*dz1) - (z2+L*dz2));
	     debug(HERE, "x can be: %f or %f diff: %f\n", x1+K*dx1, x2+L*dx2, (x1+K*dx1) - (x2+L*dx2));
	     debug(HERE, "y can be: %f or %f diff: %f\n", y1+K*dy1, y2+L*dy2, (y1+K*dy1) - (y2+L*dy2));
	     debug(HERE, "z can be: %f or %f diff: %f\n", z1+K*dz1, z2+L*dz2, (z1+K*dz1) - (z2+L*dz2));
	     debug(HERE, "not3dz\n");
	     return 0;
	 }
#endif
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
/* if (!a || !b || !l) spanic(":", HERE);*/
 l->P.x = a->x;
 l->P.y = a->y;
 l->P.z = a->z;
 l->d.x = b->x - a->x;
 l->d.y = b->y - a->y;
 l->d.z = b->z - a->z;
/* l->r = 0;*/		/* FIXME: speed optimize */
/* l->c = RED;*/
}


REAL pseudo_distance(Vertex* a, Vertex* b)
{
 REAL x,y,z;
 x = b->x - a->x;
 y = b->y - a->y;
 z = b->z - a->z;
 return x*x+y*y+z*z;
}


REAL distance(Vertex* a, Vertex* b)
{
 REAL x,y,z;
 x = b->x - a->x;
 y = b->y - a->y;
 z = b->z - a->z;
 return sqrt(x*x+y*y+z*z);
}


void normalize(Vector* v)
{
 register REAL len;
 len = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
 /* FIXME: for speed: dangerous */
/* if (nearly_equal(len, 0., 1e-9)) spanic(":", HERE);*/
 v->x /= len;
 v->y /= len;
 v->z /= len;
}


void materialize(Material* t)
{
 REAL f;
 f = t->c + t->s + t->t;
 if (f <= 1.) return;
 t->c /= f;
 t->s /= f;
 t->t /= f;
}


void normal_distorb(REAL factor, Vector* n, TexCoord* tc)
{
 REAL dist;
 dist = (REAL)rand() / (REAL)(RAND_MAX);
 dist = (dist * 2. - 1.) * factor;
 n->x += dist;
 tc->x += dist / 5.;
 dist = (REAL)rand() / (REAL)(RAND_MAX);
 dist = (dist * 2. - 1.) * factor;
 n->y += dist;
 tc->y += dist / 5.;
 dist = (REAL)rand() / (REAL)(RAND_MAX);
 dist = (dist * 2. - 1.) * factor;
 n->z += dist;
 if (tc->x < 0.) tc->x = 0.;
 if (tc->y < 0.) tc->y = 0.;
 if (tc->x > 1.) tc->x = 1.;
 if (tc->y > 1.) tc->y = 1.;
 normalize(n);
}


void get_normal_new(Vector* n, Triangle* t, Vertex* p, Material* R, Material* G, Material* B, TexCoord* tex)
{
 Ray ab,bc,ca;
 Ray ap, bp, cp;
 Vertex pa,pb,pc;
 TexCoord ttex;
 register REAL la,lb,lc;
 REAL sum,fact;
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
    if (norm_dist) normal_distorb(t->n_dist, n, tex);
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
 if (norm_dist) normal_distorb(t->n_dist, n, tex);
 /***********************/
}

void get_nurbs_pt(Triangle* t, Vector* p)
{
 REAL u,v,uu,vv;
 Ray bc;
 Ray ap;
 int sol;
 Vertex d;
 REAL lenBC, lenBD;
 REAL fact;
/* printf("npt!\n");*/
 
 make_ray(&bc, &t->b, &t->c);
 make_ray(&ap, &t->a, p);
 
 if ((sol = line_intersect(&d, &bc, &ap)) != 1)
   {
    u = t->nca.x;
    v = t->nca.y;
   }
 else
   {
    lenBC = pseudo_distance(&t->c, &t->b);
    lenBD = pseudo_distance(&d, &t->b);
    fact = lenBD / lenBC;
    fact = sqrt(fact);
    uu = (1.-fact)* t->ncb.x + fact*t->ncc.x;
    vv = (1.-fact)* t->ncb.y + fact*t->ncc.y;
    lenBC = pseudo_distance(&d, &t->a);
    lenBD = pseudo_distance(p, &t->a);
    fact = lenBD / lenBC;
    fact = sqrt(fact);
    u = (1.-fact)* t->nca.x + fact*uu;
    v = (1.-fact)* t->nca.y + fact*vv;
   }
 fastnurbs_array(t->nurbs, u, v, glob_rv, Smin);
 old_p.x = p->x;
 old_p.y = p->y;
 old_p.z = p->z;
 glob_u = u;
 glob_v = v;
 p->x = glob_rv[COORD_LX];
 p->y = glob_rv[COORD_LY];
 p->z = glob_rv[COORD_LZ];
/* printf("Was: (%Lf, %Lf, %Lf), is: (%Lf, %Lf, %Lf)\n", */
/*	 old_p.x, old_p.y, old_p.z, p->x, p->y, p->z);*/
/* printf("diff: (%Lf, %Lf, %Lf), (%Lf,%Lf)\n",*/
/*	 old_p.x - p->x, old_p.y - p->y, old_p.z - p->z, u, v);*/
}

void get_normal_nurbs_l2(Vector* n, Triangle* t, Vertex* p, Material* R, Material* G, Material* B, TexCoord* tex)
{
 Ray bc;
 Ray ap;
 int sol;
 Vertex d;
 TexCoord tc;
 REAL lenBC, lenBD;
 REAL fact;
 make_ray(&bc, &t->b, &t->c);
 make_ray(&ap, &t->a, &old_p);
 R->t = glob_rv[COORD_TR];
 R->s = glob_rv[COORD_SR];
 R->c = glob_rv[COORD_CR];
 G->t = glob_rv[COORD_TG];
 G->s = glob_rv[COORD_SG];
 G->c = glob_rv[COORD_CG];
 B->t = glob_rv[COORD_TB];
 B->s = glob_rv[COORD_SB];
 B->c = glob_rv[COORD_CB];
 materialize(R);
 materialize(G);
 materialize(B);
 calc_normal(t->nurbs, glob_u, glob_v, &n->x, &n->y, &n->z);
 if ((sol = line_intersect(&d, &bc, &ap)) != 1)
   {
    tex->x = t->ta.x;
    tex->y = t->ta.y;
    return ;
   }
 lenBC = pseudo_distance(&t->c, &t->b);
 lenBD = pseudo_distance(&d, &t->b);
 fact = lenBD / lenBC;
 fact = sqrt(fact);
 tc.x = (1.-fact)* t->tb.x + fact*t->tc.x;
 tc.y = (1.-fact)* t->tb.y + fact*t->tc.y;
 lenBC = pseudo_distance(&d, &t->a);
 lenBD = pseudo_distance(p, &t->a);
 fact = lenBD / lenBC;
 fact = sqrt(fact);
 tex->x = (1.-fact)* t->ta.x + fact*tc.x;
 tex->y = (1.-fact)* t->ta.y + fact*tc.y;
 if (norm_dist) normal_distorb(t->n_dist, n, tex);
}

void get_normal_nurbs_l1(Vector* n, Triangle* t, Vertex* p, Material* R, Material* G, Material* B, TexCoord* tex)
{
 Ray bc;
 Ray ap;
 int sol;
 Vertex d;
 Vector tn;
 TexCoord tc;
 REAL lenBC, lenBD;
 REAL fact;
 make_ray(&bc, &t->b, &t->c);
 make_ray(&ap, &t->a, &old_p);
 R->t = glob_rv[COORD_TR];
 R->s = glob_rv[COORD_SR];
 R->c = glob_rv[COORD_CR];
 G->t = glob_rv[COORD_TG];
 G->s = glob_rv[COORD_SG];
 G->c = glob_rv[COORD_CG];
 B->t = glob_rv[COORD_TB];
 B->s = glob_rv[COORD_SB];
 B->c = glob_rv[COORD_CB];
 materialize(R);
 materialize(G);
 materialize(B);
 if ((sol = line_intersect(&d, &bc, &ap)) != 1)
   {
    n->x = t->na.x;
    n->y = t->na.y;
    n->z = t->na.z;
    tex->x = t->ta.x;
    tex->y = t->ta.y;
    return ;
   }
/* printf("sol = %d\n", sol);*/
 lenBC = pseudo_distance(&t->c, &t->b);
 lenBD = pseudo_distance(&d, &t->b);
 fact = lenBD / lenBC;
 fact = sqrt(fact);
 tn.x = (1.-fact)* t->nb.x + fact*t->nc.x;
 tn.y = (1.-fact)* t->nb.y + fact*t->nc.y;
 tn.z = (1.-fact)* t->nb.z + fact*t->nc.z;
 tc.x = (1.-fact)* t->tb.x + fact*t->tc.x;
 tc.y = (1.-fact)* t->tb.y + fact*t->tc.y;
 normalize(&tn);
 lenBC = pseudo_distance(&d, &t->a);
 lenBD = pseudo_distance(p, &t->a);
 fact = lenBD / lenBC;
 fact = sqrt(fact);
 n->x = (1.-fact)* t->na.x + fact*tn.x;
 n->y = (1.-fact)* t->na.y + fact*tn.y;
 n->z = (1.-fact)* t->na.z + fact*tn.z;
 tex->x = (1.-fact)* t->ta.x + fact*tc.x;
 tex->y = (1.-fact)* t->ta.y + fact*tc.y;
 normalize(n);
/* printf("!");*/
 if (norm_dist) normal_distorb(t->n_dist, n, tex);
}

void get_normal_nurbs(Vector* n, Triangle* t, Vertex* p, Material* R, Material* G, Material* B, TexCoord* tex)
{
 if (n_pl == 1) get_normal_nurbs_l1(n, t, p, R, G, B, tex);
 else get_normal_nurbs_l2(n, t, p, R, G, B, tex);
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
 REAL lenBC, lenBD;
 REAL fact;
 make_ray(&bc, &t->b, &t->c);
 make_ray(&ap, &t->a, p);
 /*printf("P = (%f,%f,%f)\n", p->x, p->y, p->z);
 print_ray(&bc);
 print_ray(&ap);*/;
 if ((sol = line_intersect(&d, &bc, &ap)) != 1)
   {
/* if (nearly_equal(t->a.x, d.x, 1e-7) &&
	 nearly_equal(t->a.y, d.y, 1e-7) &&
	 nearly_equal(t->a.y, d.y, 1e-7))*/
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
    if (norm_dist) normal_distorb(t->n_dist, n, tex);
    return ;
   }
 lenBC = pseudo_distance(&t->c, &t->b);
/* if (nearly_equal(lenBC, 0., 1e-9)) spanic(":", HERE);*/
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
/* if (nearly_equal(lenBC, 0., 1e-9)) spanic(":", HERE);*/
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
 if (norm_dist) normal_distorb(t->n_dist, n, tex);
/* printf("gn: (%f,%f,%f)\n",p->x,p->y,p->z);*/
/* printf("intersect: (%3.0f,%3.0f,%3.0f),( %3.0f,%3.0f,%3.0f), %3.0f,%3.0f; (%1.2f,%1.2f,%1.2f)\n", p->x, p->y, p->z, d.x, d.y, d.z, lenBC, lenBD, n->x, n->y, n->z);*/
}


REAL scalar_prod(Vector* v, Vector* w)
{
 return v->x*w->x+v->y*w->y+v->z*w->z;
}


REAL length(Vector* v)
{
 return sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
}


REAL pseudo_length(Vector* v)
{
 return v->x*v->x+v->y*v->y+v->z*v->z;
}

int get_red(Screen* s, int x, int y);
int get_green(Screen* s, int x, int y);
int get_blue(Screen* s, int x, int y);

void get_texture(Texture* t, REAL x, REAL y, REAL *r, REAL* g, REAL* b)
{
 int i,j;
 i = (int)(x*(REAL)t->x);
 j = (int)(y*(REAL)t->y);
/* printf("(%f,%f) from [%d,%d] --> (%d,%d)\n", x,y,t->x,t->y,i,j);*/
 *r = (REAL)get_red(t, i, j)/255.;
 *g = (REAL)get_green(t, i, j)/255.;
 *b = (REAL)get_blue(t, i, j)/255.;
}


REAL recurse_color(Ray* r, Triangle* ptrlist, Vertex* v, Vector* n, Material* R, Material* G, Material* B, TexCoord* tC, int id, REAL f)
{
 REAL mU,mD;
 REAL rv;
 REAL ca;
 REAL fact, fact2;
 REAL a2;
 REAL nf;
 REAL add;
 REAL rt,gt,bt;
 REAL beta,arg;
 REAL l1,l2,shadow;
 int idx;
 int tmp_npl;
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
 if (r->r >= max_rec)
   {
    if (r->c == RED)   return (REAL)BACK_R/511.;
    if (r->c == GREEN) return (REAL)BACK_G/511.;
    if (r->c == BLUE)  return (REAL)BACK_B/511.;
   }
/* printf("intersect: (%f,%f,%f)\n", v->x, v->y, v->z);*/
/* printf("N: (%f,%f,%f)\n", n->x, n->y, n->z);*/
 rv = 0.;
 if (!vlight)
   {
    tr.x = light.x - v->x;
    tr.y = light.y - v->y;
    tr.z = light.z - v->z;
   }
 else
   {
    tr.x = -light.x;
    tr.y = -light.y;
    tr.z = -light.z;
   }
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
 if (!vlight) normalize(&tr);
 normalize(&obsV);
 shadow = 1.;
 tmp_npl = n_pl;
 n_pl = 0;
 if (intersection(ptrlist, &li, &nv, &idx))	/* FIXME, ptrlist is NOW unneeded because it is in btree */
   {
    if (!vlight) l1 = pseudo_distance(v, &light);
    else l1 = 1e12;	/* light in infinity */
    l2 = pseudo_distance(v, &nv);
    if (l2 < l1)
      {
       get_normal(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
       if (r->c == RED)   shadow = nR.t;
       if (r->c == GREEN) shadow = nG.t;
       if (r->c == BLUE)  shadow = nB.t;
       if (shadow > minshadow) shadow = minshadow;
      }
   }
 n_pl = tmp_npl;
 if (shadow < maxshadow) shadow = maxshadow;
/* if (shadow > 1.) shadow = 1.;*/
/* printf("nL: (%f,%f,%f)\n", tr.x, tr.y, tr.z);*/
/* printf("nV: (%f,%f,%f)\n", obsV.x, obsV.y, obsV.z);*/
 if (!l_disabled) ca = scalar_prod(&tr, n);
 else ca = 1.;
/* printf("ca = %f\n", ca);*/
 if (ca < 0. && t->faces == 1) ca = 0.;
 if (ca < 0. && t->faces == 2) ca *= -1.;
/* printf("nca = %f\n", ca);*/
 if (r->c == RED)   fact = R->c;
 if (r->c == GREEN) fact = G->c;
 if (r->c == BLUE)  fact = B->c;
 if (t->tex && !t_disabled)
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
 if (ca < ambient) ca = ambient;
 if (ca > 1.) ca = 1.;
/* printf("ca = %f\n", ca);*/
/* printf("diffuse adds: %f\n", fact * ca);*/
 rv += fact * ca * shadow;
 /* FIXME: is this specular flash OK? */
 if (t->ca > 0. && !sh_disabled)
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
	 if (ptrlist[idx].nurbs && n_pl >= 1) get_normal_nurbs(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
	 else get_normal(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
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
/*       debug(HERE, "mU --> mD\n");*/
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
/*       debug(HERE, "mD --> mU\n");*/
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
      pn.x = n->x * arg;	/* FIXME: OPTIMIZE: was cos(beta) */
      pn.y = n->y * arg;
      pn.z = n->z * arg;
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
	 if (ptrlist[idx].nurbs && n_pl >= 1) get_normal_nurbs(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
	 else get_normal(&nn, &ptrlist[idx], &nv, &nR, &nG, &nB, &ntC);
	 add = nf * recurse_color(&ti, ptrlist, &nv, &nn, &nR, &nG, &nB, &ntC, idx, nf);
	 rv += add;
       }
   }
 ret_val:
/* else printf("skipped!\n");*/
/* if (rv <= step) return 0.;*/		/* FIXME: may be needed to uncomment this */
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
 REAL rr,gg,bb;
 int idx;
#ifndef NOGL
#ifndef NOSIGNALS
 if (use_gl) while (proc_signal) sleep(1);
#endif
#endif
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
       if (t[idx].nurbs && n_pl >= 1) get_normal_nurbs(&n, &t[idx], &v, &Rm, &Gm, &Bm, &tC);
       else get_normal(&n, &t[idx], &v, &Rm, &Gm, &Bm, &tC);
       r->c = RED;
       rr = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1.);
       r->c = GREEN;
       gg = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1.);
       r->c = BLUE;
       bb = recurse_color(r, t, &v, &n, &Rm, &Gm, &Bm, &tC, idx, 1.);
/*       if (rr > .99)*/
/*       fprintf(dbg,"COMPUTED: RGB: (%f,%f,%f) at (%d,%d)\n", rr,gg,bb,x,y);*/
       set_color(s, x, y, (int)(rr*255.), (int)(gg*255.), (int)(bb*255.));
/*       spanic(":", HERE);*/		/* BIG_PANIC */
/*       printf("*");*/
   }
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
 if (la < 3e-4 || lb < 3e-4)	/* was 1e-7 with length not pseudo_length */
   {
    printf("Warning: 0 length edge\n");
    return 0.;
   }
 normalize(&a);
 normalize(&b);
 alfa = acos(scalar_prod(&a, &b));
/* return .5 * sin(alfa) * la * lb;*/
 return sin(alfa) * la * lb;
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


REAL calc_box(Triangle* t1, Triangle* t2)
{
 REAL minx,miny,minz;
 REAL maxx,maxy,maxz;
 REAL dx,dy,dz;
/* if (!t1 || !t2) spanic(":", HERE);*/
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
 REAL minx,miny,minz;
 REAL maxx,maxy,maxz;
/* if (!head || ! t1) spanic(":", HERE);*/
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


REAL calc_boxes(BTree* b1, BTree* b2)
{
 REAL minx,miny,minz;
 REAL maxx,maxy,maxz;
 REAL dx,dy,dz;
/* if (!b1 || !b2) spanic(":", HERE);*/
/* dx = dy = dz = 0.;*/
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
 dx = maxx - minx;		/* THIS WAS NOT COMPUTED UNTIL 31.09.2005 */
 dy = maxy - miny;		/* and algorithm worked, but used just two */
 dz = maxz - minz;		/* adjacent boxes when crating btree */
 			/* result: bigger boxes in btree --> more intersections */
 			/* to compute, slower algorithm: GRRRRR */
/* return sqrt(dx*dx+dy*dy+dz*dz);*/
/* printf("dx = %Lf, dy = %Lf, dz = %Lf\n", dx, dy, dz);*/
/* return 0.;*/
 return dx*dx+dy*dy+dz*dz;
}


void create_btree_box(Box** b, BTree* b1, BTree* b2)
{
 Box* temp;
 REAL minx,miny,minz;
 REAL maxx,maxy,maxz;
/* if (!b || !b1 || !b2) spanic(":", HERE);*/
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
 REAL minD,d;
 int i,j;
/* if (!head || !b1 || !b2) spanic(":", HERE);*/
 p1 = head;
 minD = 1e15;
 *b1 = *b2 = NULL;
 bl = br = NULL;
/* printf("search minimum...\n");*/
  i = j = 0;
  while (p1)
   {
    p2 = p1->next;
    while (p2)
      {
       if ( (d = calc_boxes((BTree*)p1->ptr, (BTree*)p2->ptr)) < minD)
         {
	  minD = d;
/*          printf("new minD = %Lf\n", minD);*/
          bl = p1;
	  br = p2;
         }
       j++;
       p2 = p2->next;
       if (swap == 2 && j >= apply_steps) p2 = NULL;	/* skipper for speed */
      }
    j = 0;
    i++;
    p1 = p1->next;
    if (swap == 3 && i >= apply_steps) p1 = NULL;
   }
 /*else			if swap > 1 (2 or 3)
 {
  p2 = p1->next;
  bl = p1;
  br = p2;
 }*/
/* printf("final minD = %Lf\n", minD);*/
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


void print_bd(Box* b)
{
 printf("%p ==> %Lf\n",(void*)b, b->maxx-b->minx+b->maxy-b->miny+b->maxz-b->minz);
}


void print_btree(BTree* b)
{
 if (!b) return;
 if (!b->l && !b->r)
   {
    print_bd(b->b);
   }
 else
  {
   printf("L"); print_btree(b->l);
   printf("R"); print_btree(b->r);
  }
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
/* if (!head) spanic(":", HERE);*/
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
/* while ((no = count_items(phead)) > 1)*/
 while (no > 1)
   {
/*       printf("creating btree... items: %d\n", no);*/
    if (!(no % 32)) { printf("."); fflush(stdout); }
    temp = find_nearest_boxes(phead, &bt1, &bt2, (no % 2) + skip_minimalize_algorithm);
    ptrlist_delete(&phead, bt1);
    ptrlist_delete(&phead, bt2);
    ptrlist_add(&phead, (void*)temp);
    no --;
   }
 btree = (BTree*)phead->ptr;
 ptrlist_free(&phead);
/* print_btree(btree);*/
}


void preprocess_scene(Triangle* tl)
{
 int i;
 int ii;
 REAL minS,s;
 PtrList *head, *hd;
 PtrList *minT;
 BList * bhead;
 Triangle *t1, *t2;
 printf("creating btree...phase1.");
 fflush(stdout);
 head = NULL;
 bhead = NULL;
 for (i=0;i<nTriangles;i++) ptrlist_add(&head, &tl[i]);
 ii = 0;
 while (head)
 {
  ii++;
  if (!(ii % 32)) { printf("."); fflush(stdout); }
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
 printf(".phase2");
 if (!skip_minimalize_algorithm) printf("(full minimalize).");
 else
   {
    apply_steps = (int)(apply_steps_perc * (REAL)nTriangles);
    printf("(partial minimalize %3.3Lf%%).", ((REAL)apply_steps*100.)/(REAL)nTriangles);
   }
 fflush(stdout);
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


void free_mem()
{
 PtrList* tmp;
 ListTransform *lt;
 free_matrix(world, 4);
 free_matrix(worldn, 4);
 if (tlist)
   {
    tmp = tlist;
    while (tmp)
      {
       lt = (ListTransform*)tmp->ptr;
/*       printf("lt = %p, (%f) (%f)\n", (void*)lt, lt->M[0][0], lt->MN[0][0]);*/
       free_matrix(lt->M, 4);
       free_matrix(lt->MN, 4);
       tmp = tmp->next;
      }
    ptrlist_free(&tlist);
   }
}


int find_last_ok_line(char* recfr)
{
 Texture t;
 int i,j,r,g,b,ok;
 if (!recfr) return 0;
 create_texture(&t, recfr, 0, OPEN_PARTIAL);
 if (t.x != screen.x || t.y != screen.y) spanic("find_last_ok_line: current screen and recover screen has different sizes", HERE);
 for (i=0;i<t.x;i++)
   {
    ok = 1;
    for (j=0;j<t.y;j++)
      {
       if ((r = get_red(&t, i, j))   != 0) ok = 0;
       if ((g = get_green(&t, i, j)) != 0) ok = 0;
       if ((b = get_blue(&t, i, j))  != 0) ok = 0;
       set_color(&screen, i, j, r, g, b);
      }
    if (ok)
      {
       printf("Recovering from line: %d\n", i);
	   free_texture(&t);
	   /* FIXME: possible error here */
       return i;
      }
   }
 spanic("find_last_ok_line: looks like there is nothing to recover", HERE);
 return 0;
}


void copy_line(unsigned char* line, Screen* s, int idx)
{
 int i,r,g,b;
 if (!line || !s || idx <0 || idx >= s->x) spanic("copy_line: bad parameters", HERE);
 for (i=0;i<s->y;i++)
   {
    r = get_red(s, idx, i);
    g = get_green(s, idx, i);
    b = get_blue(s, idx, i);
    line[3*i]   = r;
    line[3*i+1] = g;
    line[3*i+2] = b;
   }
}


void unblank_line(Screen* s, int idx, unsigned char* l)
{
 int i, r, g, b;
 if (!s || !l || idx < 0 || idx>= s->x) spanic("unblank_line: bad parameters", HERE);
 for (i=0;i<s->y;i++)
   {
    r = l[3*i];
    g = l[3*i+1];
    b = l[3*i+2];
    set_color(s, idx, i, r, g, b);
   }
}

#ifndef NOSIGNALS

void catch_signal(int signo)
{
 unsigned char* line;
 rt_hlt = 1;
 line = NULL;
 printf("\nCaught signal number: %d\n", signo);
 if (signo == SIGINT)
   {
    proc_signal = 1;
    printf("Exiting...\n");
/*    printf("proc_signal = %d\n", proc_signal);*/
#ifndef NOGL
    if (use_gl) line_idx --;
#endif
    blank_line(&screen, line_idx);
    wrt_bmp(&screen, signalledbmp);
/*    rt_hlt = 0;*/
    exit(1);
   }
 else if (signo == SIGUSR1)
   {
    if (proc_signal) spanic("catch_signal: already processign signal.\n", HERE);
    printf("Writing bitmap %s...\n", demandbmp);
    proc_signal = 1;
/*    printf("proc_signal = %d\n", proc_signal);*/
    line = (unsigned char*)malloc(screen.y * 3 + 1);
#ifndef NOGL
    if (use_gl) line_idx --;
#endif
    copy_line(line, &screen, line_idx);
    blank_line(&screen, line_idx);
    wrt_bmp(&screen, demandbmp);
    unblank_line(&screen, line_idx, line);
#ifndef NOGL
    if (use_gl) line_idx ++;
#endif
    free(line);
    line = NULL;
    proc_signal = 0;
    rt_hlt = 0;
   }
 else printf("No action taken for signal: %d\n", signo);
 rt_hlt = 0;
}


void setup_signals()
{
 /*static struct sigaction act;*/
/*    printf("sighandler_inst.\n");*/
 if (no_sig) return;
 signal(SIGINT, catch_signal);
 signal(SIGUSR1, catch_signal);
 /*act.sa_handler = catch_signal;
 sigfillset(&(act.sa_mask));
 sigaction(SIGINT , &act, NULL);
 sigaction(SIGUSR1, &act, NULL);*/
}

#endif

void load_bin_preprocessed(FILE* f, Triangle* ts)
{
 int n;
 int root;
 int i,ii;
 int r,l;
 sREAL tmp;
/* printf("Binary file format detected.\n");*/
 fread(&i, 4, 1, f);			/* to skip first BINT bytes (4) */
 fread(&n, sizeof(int), 1, f);
 if (n <= 0) spanic("load_bin_preprocessed: negative or zero nElem value", HERE);
 fread(&root, sizeof(int), 1, f);
 if (root < 0) spanic("load_bin_preprocessed: bad root value", HERE);
 if ((root != n-1) || (nTriangles != n+1)) spanic("load_preprocessed: bad root index value: file appears outdated", HERE);
 btree = (BTree*)malloc(n*sizeof(BTree));
 boxes = NULL;
/* printf("n = %d, root = %d\n", n, root);*/
 for (i=0;i<n;i++)
   {
    fread(&ii, sizeof(int), 1, f);
/*    printf("ii = %d\n", ii);*/
    if (i != ii) spanic("load_bin_preprocessed: bad idx value", HERE);
    fread(&r, sizeof(int), 1, f);
    fread(&l, sizeof(int), 1, f);
    if (r < -1 || l < -1) spanic("load_bin_preprocessed: bad R/L value", HERE);
    if (r >= 0) btree[i].r = &btree[r];
    else btree[i].r = NULL;
    if (l >= 0) btree[i].l = &btree[l];
    else btree[i].l = NULL;
    btree[i].b = (Box*)malloc(sizeof(Box));
    fread(&r, sizeof(int), 1, f);
    fread(&l, sizeof(int), 1, f);
    if (r >= 0) btree[i].b->t1 = &ts[r];
    else btree[i].b->t1 = NULL;
    if (l >= 0) btree[i].b->t2 = &ts[l];
    else btree[i].b->t2 = NULL;
    fread(&tmp, sizeof(sREAL), 1, f); btree[i].b->minx = tmp;
    fread(&tmp, sizeof(sREAL), 1, f); btree[i].b->miny = tmp;
    fread(&tmp, sizeof(sREAL), 1, f); btree[i].b->minz = tmp;
    fread(&tmp, sizeof(sREAL), 1, f); btree[i].b->maxx = tmp;
    fread(&tmp, sizeof(sREAL), 1, f); btree[i].b->maxy = tmp;
    fread(&tmp, sizeof(sREAL), 1, f); btree[i].b->maxz = tmp;
   }
 btree = &btree[root];
 fclose(f);
}


void load_preprocessed(Triangle* ts)
{
 char fn[2048];
 FILE* f;
 int nr;
 int n;
 int root;
 int i,ii;
 int r,l;
 strcpy(fn, scenef);
 if (!strstr(fn, ".")) strcat(fn, ".btree");
 else
   {
    i = 0;
    while (fn[i] != '.') i++;
    fn[i] = 0;
    strcat(fn, ".btree");
   }
 printf("Loading preprocessed scene from: %s\n", fn);
 f = fopen(fn, "r");
 if (!f) spanic("load_preprocessed: cannot read from file", HERE);
 if (is_binary(f)) { load_bin_preprocessed(f, ts); return; }
 nr = fscanf(f, "nElem: %d\n", &n);
 if (nr != 1) spanic("load_preprocessed: cannot read nElem value", HERE);
 nr = fscanf(f, "Root: %d\n", &root);
 if ((root != n-1) || (nTriangles != n+1)) spanic("load_preprocessed: bad root index value: file appears outdated", HERE);
 if (nr != 1) spanic("load_preprocessed: cannot read root value", HERE);
 btree = (BTree*)malloc(n*sizeof(BTree));
 boxes = NULL;
/* printf("n = %d\n", n);*/
 for (i=0;i<n;i++)
   {
    nr = fscanf(f, "Btree: %d\n{\n", &ii);
/*    printf("nr = %d\n", nr);*/
    if (nr != 1 || i != ii) spanic("load_preprocessed: bad idx value", HERE);
    nr = fscanf(f, " R: %d\n L: %d\n", &r, &l);
    if (nr != 2) spanic("load_preprocessed: cannot read R/L value", HERE);
    if (r >= 0) btree[i].r = &btree[r];
    else btree[i].r = NULL;
    if (l >= 0) btree[i].l = &btree[l];
    else btree[i].l = NULL;
    btree[i].b = (Box*)malloc(sizeof(Box));
    nr = fscanf(f," Box:\n {\n  Triangle: %d\n  Triangle: %d\n", &r, &l);
    if (nr != 2) spanic("load_preprocessed: cannot read t1/t2", HERE);
    if (r >= 0) btree[i].b->t1 = &ts[r];
    else btree[i].b->t1 = NULL;
    if (l >= 0) btree[i].b->t2 = &ts[l];
    else btree[i].b->t2 = NULL;
    nr = fscanf(f,"  (%Lf,%Lf,%Lf)\n", &btree[i].b->minx, &btree[i].b->miny, &btree[i].b->minz);
    if (nr != 3) spanic("load_preprocessed: cannot read minimal box vaules", HERE);
    nr = fscanf(f,"  (%Lf,%Lf,%Lf)\n", &btree[i].b->maxx, &btree[i].b->maxy, &btree[i].b->maxz);
    if (nr != 3) spanic("load_preprocessed: cannot read maximal box vaules", HERE);
    fscanf(f, " }\n}\n");
   }
 btree = &btree[root];
 fclose(f);
}

void btree_height(BTree* b, int level)
{
 if (level > hlevel) hlevel = level;
 if (b->l) btree_height(b->l, level+1);
 if (b->r) btree_height(b->r, level+1);
}

void btree_info()
{
 BTree* tmp;
 hlevel = 0;
 tmp = btree;
 btree_height(tmp, 1);
 printf("Maximal btree height: %d\n", hlevel);
}


void raytrace(Screen* s, char* recfr)
{
 Ray r;
 Triangle* ts;
 int i,j;
 int idx;
 REAL load1, load2, load3;
#ifndef NOSIGNALS
 setup_signals();
#endif
 load1 = load2 = load3 = 0.;
 step = 2./256.;
/* nmalloc = nfree = 0;*/
 load_scene(&ts, scenef);
/* print_matrix(world, 4);*/
/* print_matrix(worldn, 4);*/
 btree = NULL;
 if (!want_load_pps) preprocess_scene(ts);
 else load_preprocessed(ts);
 free_mem();
 btree_info();
 g_ts = ts;
 save_preprocessed();
 loaded = 1;
 get_normal = get_normal_old;
 if (nTriangles <= 4) intersection = intersection_old;
 else intersection = intersection_new;
 r.P.x = observer.x;
 r.P.y = observer.y;
 r.P.z = observer.z;;
 r.r = 0;
 r.d.z = (s->x+s->y)*lookz;	/* FIXME: mozna zmieniac aby dokladniej zobaczyc */
 printf("\n");
 line_idx = 0;
 /*calculate_color(s, &r, ts, 400, 300);*/
 /*printf("TEST1\n");
 r.d.x = 255 - s->x/2.;
 r.d.y = 191 - s->y/2.;
 calculate_color(s, &r, ts, 255, 191);
 spanic(":", HERE);*/
 /*printf("TEST\n");
 r.d.x = 256 - s->x/2.;
 r.d.y = 256 - s->y/2.;
 calculate_color(s, &r, ts, 256, 256);
 spanic(":", HERE);*/
 idx = find_last_ok_line(recfr);
 printf("Load format: btree_intr/triang_intr/all_possible_intr;ipp_btree(ipp_tr)\n");
 printf("line: %%lines,%%btreeintr(%%intersect),Kbtint(Kint)/Kallintr;ipp_bt(ipp_tr)\n");
 for (i=idx;i<s->x;i++)
   {
/* printf("n_pl = %d\n", n_pl);*/
/*	   printf("amient = %Lf\n", ambient);*/
    r.d.x = i - s->x/2.;
    printf("%04d/%04d", i+1,s->x);
    proc_tr = proc_tr2 = proc_bt = 0.0;
    line_idx = i;
    for (j=0;j<s->y;j++)
      {
       if (!(j % 24)) { printf("."); fflush(stdout); }
       r.d.y = j - s->y/2.;
 #ifndef NOGL
       while (rt_hlt) usleep(50000);
 #endif
       calculate_color(s, &r, ts, i, j);
      }
    load1 += proc_tr;
    load2 += proc_tr2;
    load3 += proc_bt;
    printf("%3.1Lf%%,%3.1Lf%%(%3.1Lf%%),%dK(%dK)/%dK;%1.0Lf(%1.0Lf)\n", 
	    ((REAL)(i+1)*100.)/(REAL)s->x,(proc_bt*100.)/proc_tr2,(proc_tr*100.)/proc_tr2,(int)(proc_bt/1000.),
	    (int)(proc_tr/1000.), (int)(proc_tr2/1000.),
	    proc_bt/(REAL)s->y, proc_tr/(REAL)s->y
	    );
    if (i > 0 && !(i % bkup))
      {
       wrt_bmp(&screen, partialbmp);
      }
/*    printf("\n");*/
   }
 printf("\n");
 printf("Load: btree_intr/triang_intr/all_possible_intr\n");
 printf("Final Load: %3.3Lf%%,%3.3Lf%%,%13.0Lf(%13.0Lf)%13.0Lf\n", (load3*100.)/load2,(load1*100.)/load2,load3,load1,load2);
 printf("IPP: %Lf(%Lf)/%Lf\n", load3/(REAL)(s->x*s->y), load1/(REAL)(s->x*s->y), load2/(REAL)(s->x*s->y));
 wrt_bmp(&screen, screenbmp);
#ifndef NOGL
 if (use_gl)
   {
    printf("RT finished. wait for OpenGL thread. . .\n");
    if (thread) pthread_join(thread, NULL);
   }
#endif
 free_scene(&ts);
/* if (glpixels) free(glpixels);*/
 fclose(dbg);
 dbg = NULL;
}


int get_red(Screen* s, int x, int y)
{
    /* FIXME: fast version */
 if (x == s->x) x--;
 if (y == s->y) y--;
 return s->pixels[3*(s->y * x + y)];
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)];
   }
 else spanic("get_red: indices out of range", HERE);
 return -1;
}


int get_green(Screen* s, int x, int y)
{
    /* FIXME: fast version */
 if (x == s->x) x--;
 if (y == s->y) y--;
 return s->pixels[3*(s->y * x + y)+1];
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)+1];
   }
 else spanic("get_green: indices out of range", HERE);
 return -1;
}


int get_blue(Screen* s, int x, int y)
{
    /* FIXME: fast version */
 if (x == s->x) x--;
 if (y == s->y) y--;
 return s->pixels[3*(s->y * x + y)+2];
 if (x >= 0 && x < s->x && y >= 0 && y < s->y)
   {
    return s->pixels[3*(s->y * x + y)+2];
   }
 else spanic("get_blue: indices out of range", HERE);
 return -1;
}


void wrt_bmp(Screen* s, char* out_f)
{
 BMPTag bm_handle;
 FILE* plik;
#ifndef NOJPEG
 char jpeg_out[1024];
 char gjpeg_out[1024];
#endif
 char tmp[1024];
 int i,j;
 init_bmp(&bm_handle);
 plik = fopen(out_f,"wb");
 if (!plik)
   {
    printf("Error writing BMP: %s\n", out_f);
/*    if (strcmp(out_f, "panicbmp")) spanic("wrt_bmp: cannot write to file",HERE);*/
    exit(1);
   }
 fprintf(plik,"%c%c",'B', 'M');
 if (!aa)
   {
    bm_handle.bm_y = s->y;
    bm_handle.bm_x = s->x;
   }
 else
   {
    bm_handle.bm_y = s->y/2;
    bm_handle.bm_x = s->x/2;
   }
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
 if (!aa)
  {
   for (i=0;i<bm_handle.bm_y;i++)  for (j=0;j<bm_handle.bm_x;j++)
    fprintf(plik,"%c%c%c", get_blue(s, j, i), get_green(s, j, i), get_red(s, j, i));
  }
 else
  {
/*   printf("Writing AntiAliased result.\n");*/
   for (i=0;i<bm_handle.bm_y;i++)  for (j=0;j<bm_handle.bm_x;j++)
   {
    fprintf(plik,"%c%c%c",
    (get_blue(s, 2*j, 2*i)+get_blue(s,2*j+1,2*i)+
     get_blue(s, 2*j, 2*i+1)+get_blue(s, 2*j+1, 2*i+1))/4,
    (get_green(s, 2*j, 2*i)+get_green(s,2*j+1,2*i)+
     get_green(s, 2*j, 2*i+1)+get_green(s, 2*j+1, 2*i+1))/4,
    (get_red(s, 2*j, 2*i)+get_red(s,2*j+1,2*i)+
     get_red(s, 2*j, 2*i+1)+get_red(s, 2*j+1, 2*i+1))/4);
   }
  }
 fclose(plik);
 printf("Bitmap: %s written.\n", out_f);
 strcpy(tmp, out_f);
 if (aa_bkup && aa)
   {
    if (!strstr(out_f, ".")) strcat(out_f, "_ab.bmp");
    else
      {
       i =0 ;
       while (out_f[i] != '.') i++;
       out_f[i] = 0;
       strcat(out_f, "_ab.bmp");
      }
 plik = fopen(out_f,"wb");
 if (!plik)
   {
    printf("Error writing antialiased backup BMP: %s\n", out_f);
    exit(1);
   }
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
 printf("Antialiased backup bitmap: %s written.\n", out_f);
   }
strcpy(out_f, tmp);
#ifndef NOJPEG
 if (use_jpeg)
   {
    strcpy(jpeg_out, out_f);
    if (strstr(out_f, "."))
      {
       i = 0;
       while (jpeg_out[i] != '.') i++;
       jpeg_out[i] = 0;
       strcpy(gjpeg_out, jpeg_out);
       strcat(jpeg_out, ".jpeg");
       strcat(gjpeg_out, "_gs.jpeg");
      }
    else
      {
       strcpy(gjpeg_out, jpeg_out);
       strcat(jpeg_out, ".jpeg");
       strcat(gjpeg_out, "_gs.jpeg");
      }
    plik = fopen(jpeg_out, "wb");
    if (!plik) { printf("Cannot write to: %s\n", jpeg_out); exit(1); }
    save_jpeg_file(screen.pixels, screen.x, screen.y, plik);
    fclose(plik);
    printf("JPEG: %s written.\n", jpeg_out);
    if (want_gjpeg)
      {
       plik = fopen(gjpeg_out, "wb");
       if (!plik) { printf("Cannot write to: %s\n", gjpeg_out); exit(1); }
       save_gray_jpeg_file(screen.pixels, screen.x, screen.y, plik);
       fclose(plik);
       printf("Grayscale JPEG: %s written.\n", gjpeg_out);
      }
   }
#endif
}


void help(char* fn)
{
 printf("\n\nusage: %s\n", fn);
 printf("\toption; purpose [default value]\n");
 printf("Options:\n");
 printf("\t-i input.dat; scenefile to use      [scene.dat]\n");
 printf("\t-o output.bmp; write result to      [screen.bmp]\n");
 printf("\t-U debug.out; debug file            [debug.out]\n");
#ifndef NOSIGNALS
 printf("\t-D demand.bmp; write when SIGUSR1   [ondemand.bmp]\n");
 printf("\t-S signal.bmp; write interrupted to [signalled.bmp]\n");
#endif
 printf("\t-P fatalerr.bmp; write when error   [panic.bmp]\n");
 printf("\t-p part.bmp; write temporary result [partial.bmp]\n");
 printf("\t-T prefix; texture prefix search    [off,example -T '../']\n");
 printf("\t-r N; maximal recurse level         [6]\n");
 printf("\t-R recov.bmp; continue RT of file   [NULL=do not recover]\n");
 printf("\t-b N; how often write partial.bmp   [64 or read from scene.dat] \n");
 printf("\t-x N; overwrite scene width         [0=use loaded value]\n");
 printf("\t-y N; overwrite scene height        [0=use loaded value]\n");
 printf("\t-X xdistorber (internal)            [1/17]\n");
 printf("\t-Y ydistorber (internal)            [1/19]\n");
 printf("\t-s N; random seed value             [0]\n");
 printf("\t-n percent; global normal distorber [0]\n");
 printf("\t-d num; nurbs processing level      [0,read more]\n");
#ifndef NOGL
 printf("\t-t N; time period to OGL flush      [500 in milisec]\n");
#endif
 printf("\t-m minshadow; minimal shadow cast   [0.1]\n");
 printf("\t-M maxshadow; maximal shadow cast   [0.5]\n");
 printf("\t-a ambient; ambient light value     [0.3,not shadow specific]\n");
#ifndef NOJPEG
 printf("\t-q N; JPEG quality [1-100]          [90]\n");
#endif
 printf("\t-F perc; apply only %% steps of min  [off,may gen not optimal btree]\n");
 printf("\t-B ; enables binary scene gen.      [scene.dat=>scene.bin, off]\n");
 printf("\t-A ; generates antialiased output   [half resolution, off]\n");
 printf("\t-K ; enables antialias full backup  [screen_ab.bmp, off]\n");
 printf("\t-2 ; doubles resolution after all   [off]\n");
 printf("\t-l ; disable light on scene         [off=light enabled]\n");
 printf("\t-e ; disable texturing on scene     [off=textures enabled]\n");
 printf("\t-I ; disable shiness flashes        [off=shiness enabled]\n");
#ifndef NOJPEG
 printf("\t-J ; try to use JPEGs if no BMPs    [off, on=scene.jpeg]\n");
 printf("\t-g ; generate grayscale JPEG too    [off,on=scene_gs.jpeg]\n");
#endif
#ifndef NOGL
 printf("\t-G ; use OpenGL GUI to display      [time consuming, off]\n");
 printf("\t-f ; generate fast preview of scene [off,yes=not performs RT]\n");
#endif
#ifndef NOSIGNALS
 printf("\t-L ; do not use signal handlers     [default=use,off]\n");
#endif
 printf("\t-C ; save preprocessed scene file   [off,scene.btree]\n");
 printf("\t-E ; save preprocessed scene binary [off,scene.btree]\n");
 printf("\t-c ; load preprocessed scene file   [off,scene.btree]\n");
 printf("\t-N ; enable normal distorbers       [off,on=loaded from scene]\n");
 printf("\t-hH; displays help                  [terminates proram,off]\n\n");
 printf("option -d selects algorithm used to RT NURBS:\n");
 printf("0: means: triangulate nurbs and treat as triangle set\n");
 printf("1: pretriangulate and use triangles to get nurbs u,v coords\n");
 printf("the compute all other params using nurbs at this u,v\n");
 printf("normal calculation is not done, usig flat tringle's\n");
 printf("3: same as 2 plus normal calculation, options 1 & 2 are slow\n");
 printf("option -l disables light/normal angle computation\n");
 printf("but computes shadow and shiness, this option only\n");
 printf("makes image brighter, to disable shadow use -M 1 -m 0\n");
 printf("to disable all light specific computations: -l -I -m 0 -a 1 -M 0\n");
 printf("-M has higher priority than -m, so -M .1 -m .5, sets max to .1\n");
 printf("Using option -n also enables -N\n");
 printf("Transformations are applied from backward:\n");
 printf("Triangle transf, list transforms (from back), world transform\n");
 printf("Normal distorbs: if no local (tr or list) then global\n");
 printf("local priority: first own then list\n");
 printf("global priority: first world then below screen\n");
 printf("-n option overrides only global normal distorbers\n");
 printf("group options: -L(-S-D), -N(-n-s), -G(-t), -J(-q-g), -A(-K-2)\n");
 printf("Warning: optimize verbose, files *.btree shouldn't\n");
 printf("be modified by the users... they are for internal usage\n");
 printf("bad indices in BTREE files can cause to irrational program\n");
 printf("behaviour, because tree-connections are not checked\n");
 printf("to speed-up intersection computation\n");
 printf("Option -F 0 decreases btree gen algorithm from oN^3 to oN^2\n");
 printf("But for 12000 triangles it loads 48952K intr per line (IPL)\n");
 printf("In comparision to 7598K IPL with N^3  (50mins of preprocessing)\n");
 printf("Result of preprocessing can be saved using -E or -C option\n");
 printf("With 140 triangles n3 -> 209K IPL, n2 -> 292K IPL\n");
 printf("-F %%, =k%% then stop each minimalize step\n");
 printf("after k%% boxes calculated, so -F 0 means WORSE algorithm n^2\n");
 printf("-F 100 means full algorithm (n^3) as if without -F\n");
 printf("HINT: option -F 8 gives already very good results\n");
 printf("HINT: setting -F to more than 25 improver only a bit\n");
#ifndef NOGL
 printf("In OpenGL GUI:\n");
 printf("\t-/+(=): changes refresh times period\n");
#ifndef NOSIGNALS
 printf("\ts:      send SIGUSR1 to RT thread\n");
 printf("\tk:      send SIGINT to RT thread\n");
#endif
#ifndef NOJPEG
 printf("\tj/J:    manipulate JPEG quality\n");
 printf("\tg:      toggle grayscale\n");
#endif
 printf("\tp:      toggle JPEG generation\n");
 printf("\ta/A:    manipulate ambient\n");
 printf("\tm/M:    manipulate minimal shadow value\n");
 printf("\tv/V:    manipulate maximal shadow value\n");
 printf("\tb/B:    maniputale backup per lines\n");
 printf("\tr/R:    manicpulate recursion level\n");
 printf("\tuU/iI:  manipulate all triangles distorber\n");
 printf("\tt:      toggle terminate RT\n");
 printf("\td:      toggle normal distorbs\n");
 printf("\th:      display help\n");
 printf("\tn:      toggle function interpolating normal, DANGEROUS\n");
 printf("\tl:      toggle function computing intersection, DANGEROUS\n");
 printf("\t1:      toggle lighting enable/disable\n");
 printf("\t2:      toggle texturing enable/disable\n");
 printf("\t3:      toggle shiness enable/disable\n");
 printf("In OpenGL Preview:\n");
 printf("When -f option used 'fast preview'\n");
 printf("You can transform scene and if You are happy with result\n");
 printf("exit and add there transformations to the world_transform\n");
 printf("option -t is also applicable in preview mode\n");
 printf("just forces maximal framerate: -t 100 gives max 10 FPS\n");
 printf("textures in preview mode are hardware generated\n");
 printf("specially x and y sizes are truncated to power of 2\n");
 printf("so if You are using non-standard textures in RT\n");
 printf("these textures may look a bit different in preview mode\n");
 printf("preview mode uses all possible optimizations for speed\n");
 printf("so don't expect scene generated to be 100%% correctly\n");
 printf("rendered: for example triangle priority sorting is\n");
 printf("stripped to just z-middle sorting, transparency is\n");
 printf("implemented using alpha channel (just blending pixels\n");
 printf("directly under currently computed) specular is simply\n");
 printf("disabled, shiness flashes are skipped etc.\n");
 printf("for very flat objects z-middle algorithm can generate\n");
 printf("bad effects - remember: preview is only for generating\n");
 printf("some interesting scene transforms, then saving its to\n");
 printf("world_trans.dat and apply in real RT, 'w' key doing this\n");
 printf("\tZzXxCc: translate +/- X,Y,Z\n");
 printf("\taAsSdD: scale +/- X,Y,Z\n");
 printf("\t123456: rotate +/- X,Y,Z\n");
 printf("\tlLtTbB: enable/disable light,texturing,blending\n");
 printf("\tuUiIoO: light translations\n");
 printf("\tn:      invert normals\n");
 printf("\ty:      priority sort toggle\n");
 printf("\tp:      display current transforms configuration\n");
 printf("\tw:      write world transform to: world_trans.dat\n");
#endif
}


void check_options()
{
 if (skip_minimalize_algorithm)
   {
    if (apply_steps_perc < 0.) apply_steps_perc = 0;
    if (apply_steps_perc > 1.) apply_steps_perc = 1.;
   }
 if (n_pl < 0) n_pl = 0;
 if (n_pl > 2) n_pl = 2;
 if (n_pl > 0) 
     spanic("This algorithm is not supported yet\n"
	     "Must implement NURBStransform with its triangles\n"
	     "List transform handler/detector, interpolation triint\n"
	     "May be different than NURBS computed u,v so this alg\n"
	     "May not ever work at all, need implement a REAL intr\n"
	     "With a NURBS surface, but it will be very slow"
	     , HERE);
/* printf("n_pl = %d\n", n_pl);*/
 if (global_dist < 0.) global_dist = 0.;
 if (aa_bkup && !aa) aa_bkup = 0;
 if (want_gjpeg && !use_jpeg) want_gjpeg = 0;
 if (jqual < 1) jqual = 1;
 if (jqual > 100) jqual = 100;
 if (ambient < 0) ambient = 0.;
 if (ambient > 1.) ambient = 1;
 if (minshadow < 0.) minshadow = 0.;
 if (minshadow > 1.) minshadow = 1.;
 if (maxshadow < 0.) maxshadow = 0.;
 if (maxshadow > 1.) maxshadow = 1.;
 if (bkup < 0) spanic("check_options: negative bkup value", HERE);
 if (ovr_x < 0 || ovr_y < 0) spanic("check_options: ovr_x or ovr_y negative", HERE);
 if (!strcmp(scenef,"")) spanic("check_options: empty scenef", HERE);
 if (!strcmp(screenbmp,"")) spanic("check_options: empty screenbmp", HERE);
 if (!strcmp(demandbmp,"")) spanic("check_options: empty demandbmp", HERE);
 if (!strcmp(signalledbmp,"")) spanic("check_options: empty signalledbmp", HERE);
 if (!strcmp(partialbmp,"")) spanic("check_options: ewmpty partialbmp", HERE);
 if (!strcmp(panicbmp,"")) spanic("check_options: empty panicbmp", HERE);
#ifdef NOJPEG
 use_jpeg = 0;
#endif
#ifdef NOGL
 use_gl = 0;
#endif
 srand(rseed);
}

void run_fast_preview();

void parse_options(int lb, char** par, char** recfr)
{
 char u;
 int fast_preview;
 int tmout;
 fast_preview = 0;
 proc_signal = 0;
 minshadow = .9;	/* mean 1 - .8 = .2 shadow is always cast even if full transparent */
 maxshadow = .5;
 ambient = .3;
 bkup = 64;
 jqual = 90;
 timeout = 500000;
 tmout = 0;
 use_jpeg = 0;
 want_gjpeg = 0;
 use_gl = 0;
 aa = 0;
 want_bin_pps = 0;
 want_save_pps = 0;
 want_load_pps = 0;
 skip_minimalize_algorithm = 0;
 apply_steps = 0;
 apply_steps_perc = -1.;
 g_ts = NULL;
 nNURBS = 0;
 nTriNURBS = 0;
 g_nurbs = NULL;
 global_dist = 0;
 norm_dist = 0;
 aa_bkup = 0;
 rt_hlt = 0;
 double_res = 0;
 *recfr = NULL;
 max_rec = 6;
 want_bin = ovr_b = ovr_r = 0;
 no_sig = 0;
 ovr_m = ovr_a = 0;
 ovr_x = ovr_y = 0;
 ovr_ms = 0;
 ovr_no = 0;
 n_pl = 0;
 l_disabled = t_disabled = sh_disabled = 0;
#ifndef NOGL
 thread = 0;
#endif
 glpixels = NULL;
 rseed = 0;
 distorber1 = 1./17.;
 distorber2 = 1./19.;
 old_p.x = old_p.y = old_p.z = 0.;
 glob_u = glob_v = 0.;
 strcpy(scenef, "scene.dat");
 strcpy(screenbmp, "screen.bmp");
 strcpy(demandbmp, "ondemand.bmp");
 strcpy(partialbmp, "partial.bmp");
 strcpy(panicbmp, "panic.bmp");
 strcpy(signalledbmp, "signalled.bmp");
 strcpy(tprefix, "");
 strcpy(debugfile, "debug.out");
#ifndef NOGETOPT
 while ((u = getopt(lb,par,"gflehcECHNKI2BLAJGb:t:r:R:o:S:p:P:i:x:y:D:X:Y:m:M:a:q:s:n:T:F:U:d:"))!=-1)
 {
 /* printf("option '%c', optarg: '%s'\n", u, optarg);*/
  switch (u)
   {
    case 'I': sh_disabled = 1; break;
    case 'h': help(par[0]); exit(1); break;
    case 'H': help(par[0]); exit(1); break;
    case 'f': fast_preview = 1; break;
    case 'l': l_disabled = 1; break;
    case 'e': t_disabled = 1; break;
    case 'B': want_bin = 1; break;
    case 'C': want_save_pps = 1; break;
    case 'E': want_save_pps = 1; want_bin_pps = 1; break;
    case 'c': want_load_pps = 1; break;
    case 'F': skip_minimalize_algorithm = 2; apply_steps_perc = atof(optarg)/100.; break;
    case 'N': norm_dist = 1; break;
    case 'K': aa_bkup = 1; break;
    case 'L': no_sig = 1; break;
    case 'A': aa = 1; break;
    case '2': double_res = 1; break;
    case 'G': use_gl = 1; break;
    case 'J': use_jpeg = 1; break;
    case 'g': want_gjpeg = 1; break;
    case 'b': ovr_b = 1; bkup = atoi(optarg); break;
    case 'q': jqual = atoi(optarg); break;
    case 'd': n_pl = atoi(optarg); break;
    case 't': tmout = 1; timeout = 1000*atoi(optarg); break;
    case 'r': ovr_r = 1; max_rec = atoi(optarg); break;
    case 'm': ovr_m = 1;  minshadow = 1. - atof(optarg); break;
    case 'M': ovr_ms = 1; maxshadow = 1. - atof(optarg); break;
    case 'a': ovr_a = 1; ambient = atof(optarg); break;
    case 'n': norm_dist = 1; ovr_no = 1; global_dist = atof(optarg)/ 100.; break;
    case 'x': ovr_x = atoi(optarg); break;
    case 'y': ovr_y = atoi(optarg); break;
    case 'X': distorber1 = atof(optarg); break;
    case 'Y': distorber2 = atof(optarg); break;
    case 'R':  *recfr = malloc(strlen(optarg)+1); strcpy(*recfr, optarg); break;
    case 'i': strcpy(scenef, optarg); break;
    case 'U': strcpy(debugfile, optarg); break;
    case 'o': strcpy(screenbmp, optarg); break;
    case 'D': strcpy(demandbmp, optarg); break;
    case 'S': strcpy(signalledbmp, optarg); break;
    case 'p': strcpy(partialbmp, optarg); break;
    case 'P': strcpy(panicbmp, optarg); break;
    case 's': rseed = atoi(optarg); break;
    case 'T': strcpy(tprefix, optarg); break;
    default: printf("%s: Unrecognized option\n",par[0]); exit(1);
   }
 }
#else
 printf("getopt() not supported, compile in Your options here (hardcode).\n");
#endif
 check_options();
 dbg = fopen(debugfile, "wb");
 if (!dbg) spanic("raytrace: cannot open debug file", HERE);
 if (*recfr && !strcmp(*recfr, "")) spanic("check_options: empty recfr", HERE);
#ifndef NOGL
 if (fast_preview) run_fast_preview(tmout);
#endif
}

#ifndef NOGL
void help(char*);

void all_distorb(REAL di)
{
 int i;
 for (i=0;i<nTriangles;i++)
   {
    g_ts[i].n_dist += di;
    if (g_ts[i].n_dist < 0.) g_ts[i].n_dist = 0.;
    if (g_ts[i].n_dist > 1.) g_ts[i].n_dist = 1.;
   }
}


void thr_keyboard(unsigned char key, int x, int y)
{
/*    printf("key!\n");*/
 switch (key)
   {
        case 27: case 'q':
	    glutDestroyWindow(glutGetWindow());
	    if (glpixels) free(glpixels);
	    glpixels = NULL;
	    thread = 0;
	    use_gl = 0;
	    rt_hlt = 0;
	    printf("OpenGL GUI thread halted.\n");
	    pthread_exit(NULL);
	    break;
	case '1':
	    l_disabled = ! l_disabled;
	    printf("Lighting enabled: %d\n", !l_disabled);
	    break;
	case '2':
	    t_disabled = ! t_disabled;
	    printf("Texturing enabled: %d\n", !t_disabled);
	    break;
	case '3':
	    sh_disabled = ! sh_disabled;
	    printf("Shinnes enabled: %d\n", !sh_disabled);
	    break;
        case '-':
	    timeout -= 50000;
	    if (timeout<50000) timeout = 50000;
	    printf("New timeout is %d ms\n", timeout/1000);
	    break;
	case 'h':
	    help("rays-OpenGL-GUI");
	    break;
	case 'u':
	    if (!norm_dist) { printf("Enable normal distorers first!\n"); return; }
	    printf("Decreasing distorb to all triangles -0.1 %%\n");
	    all_distorb(-0.001);
	    break;
	case 'i':
	    if (!norm_dist) { printf("Enable normal distorers first!\n"); return; }
	    printf("Increasing distorb to all triangles +0.1 %%\n");
	    all_distorb(0.001);
	    break;
	case 'U':
	    if (!norm_dist) { printf("Enable normal distorers first!\n"); return; }
	    printf("Decreasing distorb to all triangles -5.0 %%\n");
	    all_distorb(-0.05);
	    break;
	case 'I':
	    if (!norm_dist) { printf("Enable normal distorers first!\n"); return; }
	    printf("Increasing distorb to all triangles +5.0 %%\n");
	    all_distorb(0.05);
	    break;
        case '=':
        case '+':
	    timeout += 50000;
	    if (timeout>120000000) timeout = 120000000;
	    printf("New timeout is %d ms\n", timeout/1000);
	    break;
#ifndef NOSIGNALS
        case 's':
	    printf("Flush demand from OpenGL GUI\n");
	    rt_hlt = 1;
	    kill(getpid(), SIGUSR1);
	    while (rt_hlt) usleep(100000);
	    break;
        case 'k':
	    printf("Terminate demand from OpenGL GUI\n");
	    rt_hlt = 1;
	    kill(getpid(), SIGINT);
	    while (1) sleep(1);
	    break;
#endif
#ifndef NOJPEG
	case 'j':
	    if (!use_jpeg) { printf("Enable JPEG generation first!\n"); break; }
	    jqual --;
	    if (jqual < 1) jqual = 1;
	    printf("Changed JPEG quality to: %d%%\n", jqual);
	    break;
	case 'J':
	    if (!use_jpeg) { printf("Enable JPEG generation first!\n"); break; }
	    jqual ++;
	    if (jqual > 100) jqual = 100;
	    printf("Changed JPEG quality to: %d%%\n", jqual);
	    break;
	case 'g':
	    if (!use_jpeg) { printf("Enable JPEG generation first!\n"); break; }
	    want_gjpeg = ! want_gjpeg;
	    printf("Grayscale JPEG generation: %d\n", want_gjpeg);
	    break;
	case 'p':
	    use_jpeg = ! use_jpeg;
	    printf("JPEG generation: %d\n", use_jpeg);
	    break;
#endif
	case 'a':
	    ambient -= 0.01;
	    if (ambient < 0.) ambient = 0.;
	    printf("New ambient is: %Lf\n", ambient);
	    break;
	case 'A':
	    ambient += 0.01;
	    if (ambient > 1.) ambient = 1.;
	    printf("New ambient is: %Lf\n", ambient);
	    break;
	case 'm':
	    minshadow += 0.01;
	    if (minshadow > 1.) minshadow = 1.;
	    printf("New minshadow is: %Lf\n", 1.-minshadow);
	    break;
	case 'M':
	    minshadow -= 0.01;
	    if (minshadow < 0.) minshadow = 0.;
	    printf("New minshadow is: %Lf\n", 1.-minshadow);
	    break;
	case 'v':
	    maxshadow += 0.01;
	    if (maxshadow > 1.) maxshadow = 1.;
	    printf("New maxshadow is: %Lf\n", 1.-maxshadow);
	    break;
	case 'V':
	    maxshadow -= 0.01;
	    if (maxshadow < 0.) maxshadow = 0.;
	    printf("New maxshadow is: %Lf\n", 1.-maxshadow);
	    break;
	case 'b':
	    if (bkup < 16) bkup --;
	    else if (bkup < 32) bkup -= 4;
	    else if (bkup < 64) bkup -= 8;
	    else if (bkup < 128) bkup -= 16;
	    else bkup -= 32;
	    if (bkup < 1) bkup = 1;
	    printf("New backup period is: %d\n", bkup);
	    break;
	case 'B':
	    if (bkup < 16) bkup ++;
	    else if (bkup < 32) bkup += 4;
	    else if (bkup < 64) bkup += 8;
	    else if (bkup < 128) bkup += 16;
	    else bkup += 32;
	    if (bkup > screen.x+1) bkup = screen.x+1;
	    printf("New backup period is: %d\n", bkup);
	    break;
	case 'r':
	    max_rec --;
	    if (max_rec < 0.) max_rec = 0.;
	    printf("New recursion level is: %d\n", max_rec);
	    break;
	case 'R':
	    max_rec ++;
	    if (max_rec > 64) max_rec = 64;
	    printf("New recursion level is: %d\n", max_rec);
	    break;
	case 't':
	    rt_hlt = ! rt_hlt;
	    printf("RayTracing Halt is: %d\n", rt_hlt);
	    break;
	case 'd':
	    norm_dist = !norm_dist;
	    printf("Normal distorbs enabled: %d\n", norm_dist);
	    break;
	case 'n':
           if (get_normal == get_normal_old)
	     {
	      get_normal = get_normal_new;
	      printf("New normal interpolating function: distance based: BROKEN, SLOW\n");
	     }
	   else
	     {
	      get_normal = get_normal_old;
	      printf("New normal interpolating function: two step intersection, THE RIGHT ONE\n");
	     }
           break;
	case 'l':
	   if (intersection == intersection_new)
	     {
	      intersection = intersection_old;
	      printf("New ntersection algorithm: SLOW, intersect all triangles, DEBUG\n");
	     }
	   else
	     {
	      intersection = intersection_new;
	      printf("New ntersection algorithm: FAST, BTree based, APPEARS OK\n");
	     }
	   break;
	default: printf("Unrecognized key: '%c'\n", key);
   }
}


void thr_resize_scene(int w, int h)
{
 pzoomx = (REAL)w/(REAL)thr_sx;
 pzoomy = (REAL)h/(REAL)thr_sy;
 glViewport(0, 0, (GLsizei)w, (GLsizei)h);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 glFrustum(-1., 1., -1., 1., 1.5, 20.);
 glMatrixMode(GL_MODELVIEW);
}

void preview_light_init();

void preview_resize_scene(int w, int h)
{
 GLfloat min;
 min = (w > h)?h:w;
/* glViewport(0, 0, (GLsizei)min, (GLsizei)min);*/
 glViewport(0, 0, (GLsizei)w, (GLsizei)h);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 glFrustum(-1., 1., -1., 1., 1.5, 20000.);
 glMatrixMode(GL_MODELVIEW);
}


void thr_compute_glpixels_all(unsigned char* glpix)
{
 int i,j,k;
 if (!aa)
   {
    for (i=0;i<thr_sy;i++)
    for (j=0;j<thr_sx;j++)
    for (k=0;k<3;k++)
       glpix[3*(thr_sx*i+j)+k] = screen.pixels[3*(thr_sy*j+i)+k];
   }
 else
   {
    for (i=0;i<thr_sy;i++)
    for (j=0;j<thr_sx;j++)
    for (k=0;k<3;k++)
       glpix[3*(thr_sx*i+j)+k] = (
	   screen.pixels[3*(2*thr_sy*(2*j)+(2*i))+k]+
	   screen.pixels[3*(2*thr_sy*(2*j+1)+(2*i))+k]+
	   screen.pixels[3*(2*thr_sy*(2*j)+(2*i+1))+k]+
	   screen.pixels[3*(2*thr_sy*(2*j+1)+(2*i+1))+k])/4;
   }
}


void thr_compute_glpixels(unsigned char* glpix)
{
 int i,j,k;
 int next_li;
 int halted;
 halted = rt_hlt;
 if (!halted) rt_hlt = 1;
 next_li = line_idx + 1;
/* printf(" [%d->%d]\n", prev_li, line_idx);*/
 if (!aa)
   {
    if (next_li > thr_sx) next_li = thr_sx;
    for (i=0;i<thr_sy;i++)
    for (j=prev_li;j<next_li;j++)
    for (k=0;k<3;k++)
       glpix[3*(thr_sx*i+j)+k] = screen.pixels[3*(thr_sy*j+i)+k];
   }
 else
   {
    prev_li /= 2;
    next_li /= 2;
    if (next_li > thr_sx) next_li = thr_sx;
    for (i=0;i<thr_sy;i++)
    for (j=prev_li;j<next_li;j++)
    for (k=0;k<3;k++)
       glpix[3*(thr_sx*i+j)+k] = (
	   screen.pixels[3*(2*thr_sy*(2*j)+(2*i))+k]+
	   screen.pixels[3*(2*thr_sy*(2*j+1)+(2*i))+k]+
	   screen.pixels[3*(2*thr_sy*(2*j)+(2*i+1))+k]+
	   screen.pixels[3*(2*thr_sy*(2*j+1)+(2*i+1))+k])/4;
   }
 prev_li = line_idx-1;
 if (!halted) rt_hlt = 0;
}


void thr_render_scene(void)
{
 thr_compute_glpixels(glpixels);
 glClear(GL_COLOR_BUFFER_BIT);
 glLoadIdentity();
 gluLookAt(0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
/* glRotated(10., 0., 0., 1.);*/
 glRasterPos3d(3.333333, -3.333333, 0.);
/* glPixelStorei(GL_UNPACK_ROW_LENGTH, 1);*/
 glPixelZoom((GLdouble)pzoomx, (GLdouble)pzoomy);
 glDrawPixels(thr_sx, thr_sy, GL_RGB, GL_UNSIGNED_BYTE, glpixels);
/* glRotated(90., 0., 0., 1.);*/
 glFlush();
 if (!rt_hlt) usleep(timeout);
 glutSwapBuffers();
}


void thr_anim(void)
{
 glutPostRedisplay();
}


void thr_visible(int vis)
{
 if (vis == GLUT_VISIBLE) glutIdleFunc(thr_anim);
 else                     glutIdleFunc(NULL);
}


void thr_Init()
{
 int i,mem;
 pzoomx = pzoomy = 1.;
 mem = (thr_sx+1)*(thr_sy+1)*3;
 glClearColor(0.,0.,0.,0.);
 glShadeModel(GL_SMOOTH);
 for (i=0;i<mem;i++) glpixels[i] = 0;
}


/*void preview_light_deinit()
{
   glDisable(GL_LIGHT0);
   glDisable(GL_LIGHTING);
}*/

void preview_light_init()
{
   int i;
   GLfloat lightPos[4];
   GLfloat aLight[4];
   GLfloat lightColor[] = { 1.0f, 1.0f, 1.0f, 1.0 };
   lightPos[0] = light.x + ltX;
   lightPos[1] = light.y + ltY;
   lightPos[2] = light.z + ltZ;
   lightPos[3] = (GLdouble)vlight;
   for (i=0;i<4;i++)
     {
      aLight[i] = ambient;
      if (aLight[i] < .334) aLight[i] = .334;
/*      printf("aLight[%d]=%f\n", i, aLight[i]);*/
     }
   /*for (i=0;i<nTex;i++)
      {
        printf("tex[%d] = %p\n", i, (void*)(texture[i].pixels));
     }*/
   glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
   glLightfv(GL_LIGHT0, GL_AMBIENT, aLight);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor);
   glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
}


void resize_texture_to_gl_format(Texture* t)
{
 int x,y;
 int i,j,k;
 unsigned char* mem;
 x = y = 1;
 while (x < t->x) x *= 2;
 if (x != t->x) x /= 2;
 while (y < t->y) y *= 2;
 if (y != t->y) y /= 2;
 mem = (unsigned char*)malloc(3*x*y*sizeof(unsigned char));
 /*for (i=0;i<x;i++)
 for (j=0;j<y;j++)
 for (k=0;k<3;k++)
   {
    mem[3*(i*y+j)+k] = t->pixels[3*(i*t->y+j)+k];
   }*/
 for (i=0;i<y;i++)
 for (j=0;j<x;j++)
 for (k=0;k<3;k++)
/*       mem[3*(x*i+((x-1)-j))+k] = t->pixels[3*(t->y*j+i)+k];*/
/*       mem[3*(x*((y-i)-1)+j)+k] = t->pixels[3*(t->y*j+i)+k];*/
       mem[3*(x*i+j)+k] = t->pixels[3*(t->y*j+i)+k];
 free(t->pixels);
/* printf("setting new size: %dx%d\n", x,y);*/
 t->x = x;
 t->y = y;
 t->pixels = mem;
}


void preview_bind_textures()
{
 int i;
 tids = (unsigned int*)malloc(nTex*sizeof(unsigned int));
 glEnable(GL_TEXTURE_2D);
 glGenTextures(nTex, tids);
 for (i=0;i<nTex;i++)
  {
   if (!texture[i].pixels)
     {
/*      tids[i] = 0; */
/*	 printf("skipping texture: %d tid -- %d\n", i, tids[i]);*/
      continue;
     }
/*   printf("binding new: %d) %dx%d --> %p\n", i, texture[i].x, texture[i].y, (void*)texture[i].pixels);*/
   resize_texture_to_gl_format(&texture[i]);
/*   printf("after resiz: %d) %dx%d --> %p\n", i, texture[i].x, texture[i].y, (void*)texture[i].pixels);*/
/*   printf("New tid -- %d\n", tids[i]);*/
   glBindTexture(GL_TEXTURE_2D, tids[i]);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture[i].x, texture[i].y, 0, GL_RGB, GL_UNSIGNED_BYTE, texture[i].pixels);
   glBindTexture(GL_TEXTURE_2D, tids[i]);
  }
}

void clear_transforms();

void preview_Init()
{
/* glClearColor(0.,0.,0.,1.);*/
 glClearColor((GLfloat)BACK_R/255.,(GLfloat)BACK_G/255., (GLfloat)BACK_B/255., 0.);
 glShadeModel(GL_SMOOTH);
 glEnable(GL_DEPTH_TEST);
 t1 = t2 = 0;
 fps = 0;
 invN = 1;
 l_enable = 1;
 preview_ps = 1;
 clear_transforms();
 preview_bind_textures();
 glEnable(GL_BLEND);
 glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
 pv_m = matrix(4);
 pv_mn = matrix(4);
 I_matrix(pv_m, 4);
 I_matrix(pv_mn, 4);
 pv_tbuf = (Triangle*)malloc(nTriangles*sizeof(Triangle));
}


void* thr_opengl_thread(void* dummy)
{
 char** pargl;
 int lbgl;
#ifndef NOSIGNALS
 setup_signals();
#endif
 while (!loaded) usleep(50000);
 printf("OpenGL thread started.\n");
 thr_sx = screen.x;
 thr_sy = screen.y;
 if (aa) { thr_sx /= 2; thr_sy /= 2; }
 glpixels = (unsigned char*)malloc(((thr_sx+1)*(thr_sy+1)*3)*sizeof(unsigned char));
 lbgl = 1;
 prev_li = 0;
 pargl = (char**)malloc(sizeof(char*));
 pargl[0] = (char*)malloc(32*sizeof(char));
 strcpy(pargl[0], "OpenGL Display Thread.");
 glutInit(&lbgl, pargl);
 glutInitDisplayMode(GLUT_DOUBLE);
 glutInitWindowSize(thr_sx, thr_sy);
 glutInitWindowPosition(10, 10);
 glutCreateWindow(pargl[0]);
 thr_Init();
 thr_compute_glpixels_all(glpixels);
 glutDisplayFunc(thr_render_scene);
 glutReshapeFunc(thr_resize_scene);
 glutKeyboardFunc(thr_keyboard);
 glutVisibilityFunc(thr_visible);
 glutMainLoop();
 free(pargl[0]);
 free(pargl);
 return NULL;
}


void opengl_display_thread()
{
 loaded = 0;
 pthread_create(&thread, NULL, thr_opengl_thread, NULL);
}


void clear_transforms()
{
 rotX = rotY = rotZ = 0.;
 scaX = scaY = scaZ = 1.;
 traX = traY = traZ = 0.;
 ltX = ltY = ltZ = 0.;
}


void write_config()
{
 FILE* f;
 f = fopen("world_trans.dat", "w");
 if (!f) spanic("write_config: cannot write to: world_trans.dat", HERE);
 fprintf(f,"ListTransform: [%d,%d]\n", 0, nTriangles-1);
 fprintf(f,"{\n");
 fprintf(f," Scale: (%f,%f,%f)\n", scaX, scaY, scaZ);
 fprintf(f," RotateX: %f\n", rotX);
 fprintf(f," RotateY: %f\n", rotY);
 fprintf(f," RotateZ: %f\n", rotZ);
 fprintf(f," Translate: (%f,%f,%f)\n", traX, traY, traZ);
 if (invN == -1.) fprintf(f, " NegateN:\n");
 fprintf(f,"}\n");
 fprintf(f,"LightTransform:\n");
 fprintf(f,"{\n");
 fprintf(f," Translate: (%f,%f,%f)\n", ltX, ltY, ltZ);
 fprintf(f,"}\n");
 fclose(f);
 printf("World configuration saved to: world_trans.dat\n");
}


void preview_current_config()
{
 printf("Current configuration:\n\trot(%f,%f,%f)\n\tsca(%f,%f,%f)\n\ttra(%f,%f,%f)\n\tlit(%f,%f,%f,%d)\n",
	 rotX,rotY,rotZ,scaX,scaY,scaZ,traX,traY,traZ,ltX,ltY,ltZ,invN);
}


void preview_check_trans()
{
 if (rotX < 0.) rotX += 360.;
 if (rotY < 0.) rotY += 360.;
 if (rotZ < 0.) rotZ += 360.;
 if (rotX >= 360.) rotX -= 360.;
 if (rotY >= 360.) rotY -= 360.;
 if (rotZ >= 360.) rotZ -= 360.;
}


void preview_keyboard(unsigned char key, int x, int y)
{
 preview_check_trans();
 switch (key)
   {
        case 27: case 'q':
            free_scene(&g_ts);
            free_screen(&screen);
	    if (tids) free(tids);
	    tids = NULL;
	    printf("Final configuration:\n");
	    preview_current_config();
            exit(0);
	    break;
	case 'y':
	    preview_ps = ! preview_ps;
	    printf("Priority sort enabled: %d\n", preview_ps);
	    break;
	case 'w':
            preview_current_config();
            write_config();
	    break;
	case 'p':
            preview_current_config();
	    break;
	case 'n':
	    invN = -invN;
	    printf("Inverting normals, current setting: %d\n", invN);
	    break;
	case 'b':
	    printf("Disabling blending.\n");
	    glDisable(GL_BLEND);
	    break;
	case 'B':
	    printf("Enabling blending.\n");
	    glEnable(GL_BLEND);
	    break;
	case 't':
	    printf("Disabling texturing.\n");
	    glDisable(GL_TEXTURE_2D);
	    break;
	case 'T':
	    printf("Enabling texturing.\n");
	    glEnable(GL_TEXTURE_2D);
	    break;
	case 'l':
	    printf("Disabling lighting.\n");
	    glDisable(GL_LIGHTING);
	    l_enable = 0;
	    break;
	case 'L':
	    printf("Enabling lighting.\n");
	    preview_light_init();
	    l_enable = 1;
	    break;
	case 'u':
	    ltX += 10.;
	    preview_light_init();
	    printf("New light position: (%Lf, %Lf, %Lf)\n", light.x+ltX, light.y+ltY, light.z+ltZ);
	    break;
	case 'U':
	    ltX -= 10.;
	    preview_light_init();
	    printf("New light position: (%Lf, %Lf, %Lf)\n", light.x+ltX, light.y+ltY, light.z+ltZ);
	    break;
	case 'i':
	    ltY += 10.;
	    preview_light_init();
	    printf("New light position: (%Lf, %Lf, %Lf)\n", light.x+ltX, light.y+ltY, light.z+ltZ);
	    break;
	case 'I':
	    ltY -= 10.;
	    preview_light_init();
	    printf("New light position: (%Lf, %Lf, %Lf)\n", light.x+ltX, light.y+ltY, light.z+ltZ);
	    break;
	case 'o':
	    ltZ += 10.;
	    preview_light_init();
	    printf("New light position: (%Lf, %Lf, %Lf)\n", light.x+ltX, light.y+ltY, light.z+ltZ);
	    break;
	case 'O':
	    ltZ -= 10.;
	    preview_light_init();
	    printf("New light position: (%Lf, %Lf, %Lf)\n", light.x+ltX, light.y+ltY, light.z+ltZ);
	    break;
	case 'h': case 'H':
	    help("OpenGL Preview");
            preview_current_config();
	    break;
	case ' ':
	    clear_transforms();
	    printf("Cleared all transformations.\n");
	    break;
	case 'A':
	    scaX *= 1.03;
	    printf("New scale X factor: %f\n", scaX);
	    break;
	case 'S':
	    scaY *= 1.03;
	    printf("New scale Y factor: %f\n", scaY);
	    break;
	case 'D':
	    scaZ *= 1.03;
	    printf("New scale Z factor: %f\n", scaZ);
	    break;
	case 'a':
	    scaX /= 1.03;
	    printf("New scale X factor: %f\n", scaX);
	    break;
	case 's':
	    scaY /= 1.03;
	    printf("New scale Y factor: %f\n", scaY);
	    break;
	case 'd':
	    scaZ /= 1.03;
	    printf("New scale Z factor: %f\n", scaZ);
	    break;
	case '1':
	    rotX += 3.;
	    printf("New rotate X angle: %f\n", rotX);
	    break;
	case '2':
	    rotX -= 3.;
	    printf("New rotate X angle: %f\n", rotX);
	    break;
	case '3':
	    rotY += 3.;
	    printf("New rotate Y angle: %f\n", rotY);
	    break;
	case '4':
	    rotY -= 3.;
	    printf("New rotate Y angle: %f\n", rotY);
	    break;
	case '5':
	    rotZ += 3.;
	    printf("New rotate Z angle: %f\n", rotZ);
	    break;
	case '6':
	    rotZ -= 3.;
	    printf("New rotate Z angle: %f\n", rotZ);
	    break;
	case 'Z':
	    traX -= 10.;
	    printf("New X translation: %f\n", traX);
	    break;
	case 'X':
	    traY -= 10.;
	    printf("New Y translation: %f\n", traY);
	    break;
	case 'C':
	    traZ -= 10.;
	    printf("New Z translation: %f\n", traZ);
	    break;
	case 'z':
	    traX += 10.;
	    printf("New X translation: %f\n", traX);
	    break;
	case 'x':
	    traY += 10.;
	    printf("New Y translation: %f\n", traY);
	    break;
	case 'c':
	    traZ += 10.;
	    printf("New Z translation: %f\n", traZ);
	    break;
	default: printf("Unrecognized key: '%c'\n", key);
   }
}


void preview_triangle(Triangle* t)
{
 GLdouble alfa;
 /*glVertex3d(0,0,0);
 glVertex3d(100,100,0);
 glVertex3d(0,100,0);
 return;*/
 GLfloat md[4];
 GLfloat ms[4];
 GLfloat mh[1];
 GLdouble i;
 mh[0] = t->ca;
 i = (GLfloat)invN;
/* printf("%p: (%Lf,%Lf,%Lf)\n", (void*)t, t->a.x, t->a.y, t->a.z);*/
/* printf("%p: (%Lf,%Lf,%Lf)\n", (void*)t, t->b.x, t->b.y, t->b.z);*/
/* printf("%p: (%Lf,%Lf,%Lf)\n", (void*)t, t->c.x, t->c.y, t->c.z);*/
 if (t->tid > 0) glBindTexture(GL_TEXTURE_2D, tids[t->tid-1]);
 else glBindTexture(GL_TEXTURE_2D, 0);
/* glBindTexture(GL_TEXTURE_2D, tids[0]);*/
 glBegin(GL_TRIANGLES);
 glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mh);
/* printf("t->tid = %d, bind %d\n", t->tid, tids[t->tid-1]);*/
 md[0] = t->mra.c;
 md[1] = t->mga.c;
 md[2] = t->mba.c;
 ms[0] = t->mra.s;
 ms[1] = t->mga.s;
 ms[2] = t->mba.s;
 alfa = 1. - (t->mra.t + t->mga.t + t->mba.t) / 3.;
 if (alfa < .2) alfa = .2;
/* printf("alfa = %f\n", alfa);*/
 ms[3] = md[3] = alfa;
 glColor4d(t->mra.c, t->mga.c, t->mba.c, alfa);
 glNormal3d(i*t->na.x, i*t->na.y, i*t->na.z);
 glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, md);
 glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, ms);
 glTexCoord2d(t->ta.x, t->ta.y);
 glVertex3d(t->a.x, t->a.y, t->a.z);
 md[0] = t->mrb.c;
 md[1] = t->mgb.c;
 md[2] = t->mbb.c;
 ms[0] = t->mrb.s;
 ms[1] = t->mgb.s;
 ms[2] = t->mbb.s;
 alfa = 1. - (t->mrb.t + t->mgb.t + t->mbb.t) / 3.;
 if (alfa < .2) alfa = .2;
 ms[3] = md[3] = alfa;
 glColor4d(t->mrb.c, t->mgb.c, t->mbb.c, alfa);
 glNormal3d(i*t->nb.x, i*t->nb.y, i*t->nb.z);
 glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, md);
 glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, ms);
 glTexCoord2d(t->tb.x, t->tb.y);
 glVertex3d(t->b.x, t->b.y, t->b.z);
 md[0] = t->mrc.c;
 md[1] = t->mgc.c;
 md[2] = t->mbc.c;
 ms[0] = t->mrc.s;
 ms[1] = t->mgc.s;
 ms[2] = t->mbc.s;
 alfa = 1. - (t->mrc.t + t->mgc.t + t->mbc.t) / 3.;
 if (alfa < .2) alfa = .2;
 ms[3] = md[3] = alfa;
 glColor4d(t->mrc.c, t->mgc.c, t->mbc.c, alfa);
 glNormal3d(i*t->nc.x, i*t->nc.y, i*t->nc.z);
 glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, md);
 glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, ms);
 glTexCoord2d(t->tc.x, t->tc.y);
 glVertex3d(t->c.x, t->c.y, t->c.z);
 glEnd();
}


void time_counter()
{
 char tstr[64];
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
    sprintf(tstr, "OpenGL preview FPS: %d, psort: %d", fps/(int)(t2-t1), preview_ps);
    t1 = t2;
    glutSetWindowTitle(tstr);
    fps = 0;
   }
}


REAL preview_medium_z(Triangle* t)
{
 return t->a.z + t->b.z + t->c.z;
}


void preview_priority_sort_triangles()
{
 int i,j,id;
 Triangle tmp;
 REAL z, mz;
 I_matrix(pv_m, 4);
 I_matrix(pv_mn, 4);
 m_scale(&pv_m, scaX, scaY, scaZ);
 m_scale(&pv_mn, scaX, scaY, scaZ);
 m_rotatex(&pv_m, rotX);
 m_rotatex(&pv_mn,rotX);
 m_rotatey(&pv_m, rotY);
 m_rotatey(&pv_mn,rotY);
 m_rotatez(&pv_m, rotZ);
 m_rotatez(&pv_mn,rotZ);
 m_translate(&pv_m, traX, traY, traZ);
/* transform_light(&light, pv_m, pv_mn, vlight);*/
 for (i=0;i<nTriangles;i++) memcpy(&pv_tbuf[i], &g_ts[i], sizeof(Triangle));
 for (i=0;i<nTriangles;i++) transform_triangle(&pv_tbuf[i], pv_m, pv_mn);
 for (i=0;i<nTriangles;i++)
   {
    id = i;
    mz = preview_medium_z(&pv_tbuf[i]);
    for (j=i+1;j<nTriangles;j++)
      {
       if ((z = preview_medium_z(&pv_tbuf[j])) > mz)
         {
          id = j;
	  mz = z;
         }
      }
    if (id != i)
      {
       memcpy(&tmp, &pv_tbuf[i], sizeof(Triangle));
       memcpy(&pv_tbuf[i], &pv_tbuf[id], sizeof(Triangle));
       memcpy(&pv_tbuf[id], &tmp, sizeof(Triangle));
      }
   }
}


void preview_render_scene(void)
{
 int i;
 if (l_enable) preview_light_init();
 time_counter();
 glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
 glLoadIdentity();
/* gluLookAt(0.0, 0.0, -50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);*/
/* printf("Look at: (%Lf,%Lf,%Lf) --> %Lf\n", observer.x, observer.y, observer.z, (screen.x+screen.y)*lookz);*/
 gluLookAt(observer.x, observer.y, observer.z,
	 0.0, 0.0, (screen.x+screen.y)*lookz, 0.0, 1.0, 0.0);
/*	 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);*/
/* glRotated(10., 0., 0., 1.);*/
 glScaled(-1., 1., 1.);
 if (!preview_ps)
   {
    glScalef(scaX, scaY, scaZ);
    glRotatef(rotX, 1., 0., 0.);
    glRotatef(rotY, 0., 1., 0.);
    glRotatef(rotZ, 0., 0., 1.);
    glTranslatef(traX, traY, traZ);
   for (i=0;i<nTriangles;i++) preview_triangle(&g_ts[i]);
   }
 else
  {
   preview_priority_sort_triangles();
   for (i=0;i<nTriangles;i++) preview_triangle(&pv_tbuf[i]);
  }
 if (pv_usetm) usleep(timeout);
 glFlush();
 glutSwapBuffers();
}


void run_fast_preview(int tmout)
{
 Triangle* ts;
 char** pargl;
 int lbgl;
 if (tmout) pv_usetm = 1;
 else pv_usetm = 0;
 load_scene(&ts, scenef);
 free_mem();
 g_ts = ts;
 printf("OpenGL preview started.\n");
 thr_sx = screen.x;
 thr_sy = screen.y;
 lbgl = 1;
 pargl = (char**)malloc(sizeof(char*));
 pargl[0] = (char*)malloc(32*sizeof(char));
 strcpy(pargl[0], "OpenGL Preview Thread.");
 glutInit(&lbgl, pargl);
 glutInitDisplayMode(GLUT_DOUBLE);
 glutInitWindowSize((thr_sx+thr_sy)/2., (thr_sx+thr_sy)/2.);
 glutInitWindowPosition(10, 10);
 glutCreateWindow(pargl[0]);
 preview_Init();
 glutDisplayFunc(preview_render_scene);
 glutReshapeFunc(preview_resize_scene);
 glutKeyboardFunc(preview_keyboard);
 glutVisibilityFunc(thr_visible);
/* preview_light_init();*/
 glutMainLoop();
 free(pargl[0]);
 free(pargl);
 free_scene(&ts);
 free_screen(&screen);
 if (tids) free(tids);
 tids = NULL;
 exit(0);
}

#endif

int main(int lb, char** par)
{
 char* recfr;
 dbg = NULL;
 printf("RAYS by morgoth dbma\nFile: %s, line: %d, compiled %s %s\n", __FILE__,__LINE__,__DATE__,__TIME__);
 if (lb == 1) { help(par[0]); return 0; }
 parse_options(lb, par, &recfr);
#ifndef NOGL
 if (use_gl) opengl_display_thread();
#endif
 raytrace(&screen, recfr);
 free_screen(&screen);
 if (boxes) blist_free(&boxes);
 return 0;
}

