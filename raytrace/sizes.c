#include <stdio.h>
#include <stdlib.h>
#define REAL long double
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

int main()
{
 printf("sizeof(%s) = %d\n", "uchar", (int)sizeof(unsigned char));
 printf("sizeof(%s) = %d\n", "short", (int)sizeof(short));
 printf("sizeof(%s) = %d\n", "int", (int)sizeof(int));
 printf("sizeof(%s) = %d\n", "long", (int)sizeof(long));
 printf("sizeof(%s) = %d\n", "float", (int)sizeof(float));
 printf("sizeof(%s) = %d\n", "double", (int)sizeof(double));
 printf("sizeof(%s) = %d\n", "long double", (int)sizeof(long double));
 printf("sizeof(%s) = %d\n", "void*", (int)sizeof(void*));
 printf("sizeof(%s) = %d\n", "Vertex", (int)sizeof(Vertex));
 printf("sizeof(%s) = %d\n", "Vector", (int)sizeof(Vector));
 printf("sizeof(%s) = %d\n", "Material", (int)sizeof(Material));
 printf("sizeof(%s) = %d\n", "TexCoord", (int)sizeof(TexCoord));
 printf("sizeof(%s) = %d\n", "Surface", (int)sizeof(Surface));
 printf("sizeof(%s) = %d\n", "Ray", (int)sizeof(Ray));
 printf("sizeof(%s) = %d\n", "Screen", (int)sizeof(Screen));
 printf("sizeof(%s) = %d\n", "Texture", (int)sizeof(Texture));
 printf("sizeof(%s) = %d\n", "Triangle", (int)sizeof(Triangle));
 printf("sizeof(%s) = %d\n", "BMPTag", (int)sizeof(BMPTag));
 printf("sizeof(%s) = %d\n", "Box", (int)sizeof(Box));
 printf("sizeof(%s) = %d\n", "BTree", (int)sizeof(BTree));
 printf("sizeof(%s) = %d\n", "VList", (int)sizeof(VList));
 printf("sizeof(%s) = %d\n", "PtrList", (int)sizeof(PtrList));
 printf("sizeof(%s) = %d\n", "ListTransform", (int)sizeof(ListTransform));
 printf("sizeof(%s) = %d\n", "BList", (int)sizeof(BList));
 return 0;
}
