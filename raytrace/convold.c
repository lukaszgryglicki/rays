#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#define HERE __FILE__,__LINE__
typedef long double REAL;
#define sREAL double

enum { R0 = 0x1000, X90, X180, X270, Y90, Y180, Y270, Z90, Z180, Z270 };

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

typedef struct _Screen
{
 unsigned char* pixels;
 int x,y;
} Screen;

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

typedef struct _BList
{
 struct _BList *next, *prev;
 Box b;
} BList;


/* variables */
int nTriangles;
int g_idx;
int rot;
Triangle* ts;
BList* boxes;
BTree* btree;
REAL gtx, gty, gtz, gsx, gsy, gsz;

void panic(char* why, ...)
{
 va_list ap;
 va_start(ap,why);
 printf("\nFatal error occured: ");
 vprintf(why,ap);
 printf("\n");
 va_end(ap);
 exit(1);
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
    (*head)->b.maxx = b->maxx;
    (*head)->b.maxy = b->maxy;
    (*head)->b.maxz = b->maxz;
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
  temp->b.maxx = b->maxx;
  temp->b.maxy = b->maxy;
  temp->b.maxz = b->maxz;
  temp->next = *head;
  temp->prev = NULL;
  (*head)->prev = temp;
  *head = temp;
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
}

void create_assoc_table(void** pi)
{
 BTree* tmp;
 g_idx = 0;
 tmp = btree;
 assoc_idx_tree(tmp, pi);
}


void save_box_to_file(FILE* f, Box* b)
{
 fprintf(f, " Box:\n {\n");
 if (b->t1) fprintf(f, "  Triangle: %d\n", b->t1->idx);
 else fprintf(f, "  Triangle: -1\n");
 if (b->t2) fprintf(f, "  Triangle: %d\n", b->t2->idx);
 else fprintf(f, "  Triangle: -1\n");
 fprintf(f, "  (%1.12Lf,%1.12Lf,%1.12Lf)\n  (%1.12Lf,%1.12Lf,%1.2Lf)\n",
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
 if (g_idx && (g_idx % 1024) == 0) { printf("."); fflush(stdout); }
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
    panic("write_btree_to_file: bad index associative table", HERE);
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
    panic("write_btree_to_file: bad index associative table", HERE);
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
 if (g_idx && (g_idx % 1024) == 0) { printf("."); fflush(stdout); }
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
    panic("write_btree_to_file: bad index associative table", HERE);
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
    panic("write_btree_to_file: bad index associative table", HERE);
   }
 else fwrite(&mi, sizeof(int), 1, f);
 l_done:
 save_bin_box_to_file(f, t->b);
 g_idx ++;
}



void make_stub_ts(int n)
{
 int i;
 nTriangles = n+1;
 ts = (Triangle*)malloc(nTriangles*sizeof(Triangle));
 if (!ts) panic("Cannot alocate ts");
 for (i=0;i<nTriangles;i++) ts[i].idx = i;
}



int load_bin_preprocessed(FILE* f)
{
 int n;
 int root;
 int i,ii;
 int r,l;
 sREAL tmp;
 fread(&i, 4, 1, f);
 fread(&n, sizeof(int), 1, f);
 if (n <= 0) panic("load_bin_preprocessed: negative or zero nElem value", HERE);
 fread(&root, sizeof(int), 1, f);
 if (root < 0) panic("load_bin_preprocessed: bad root value: n=%d, root=%d", n, root);
 if (root != n-1) panic("load_preprocessed: bad root index value: file appears outdated, n=%d, root=%d, nT=%d", n, root, nTriangles);
 make_stub_ts(n);
 btree = (BTree*)malloc(n*sizeof(BTree));
 boxes = NULL;
 for (i=0;i<n;i++)
   {
    if (i && (i % 1024) == 0) { printf("."); fflush(stdout); }
    fread(&ii, sizeof(int), 1, f);
    if (i != ii) panic("load_bin_preprocessed: bad idx value, i=%d, ii=%d", i, ii);
    fread(&r, sizeof(int), 1, f);
    fread(&l, sizeof(int), 1, f);
    if (r < -1 || l < -1) panic("load_bin_preprocessed: bad R/L value, r=%d, l=%d", r, l);
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
    blist_add_box(&boxes, btree[i].b);
   }
 printf("\n");
 btree = &btree[root];
 fclose(f);
 return 1;
}


int load_preprocessed(FILE* f)
{
 int nr;
 int n;
 int root;
 int i,ii;
 int r,l;
 if (is_binary(f)) return load_bin_preprocessed(f);
 nr = fscanf(f, "nElem: %d\n", &n);
 if (nr != 1) panic("load_preprocessed: cannot read nElem value", HERE);
 nr = fscanf(f, "Root: %d\n", &root);
 if (nr != 1) panic("load_preprocessed: cannot read root value", HERE);
 if (root != n-1) panic("load_preprocessed: bad root index value: file appears outdated, n=%d, root=%d", n, root);
 make_stub_ts(n);
 btree = (BTree*)malloc(n*sizeof(BTree));
 boxes = NULL;
 for (i=0;i<n;i++)
   {
    if (i && (i % 1024) == 0) { printf("."); fflush(stdout); }
    nr = fscanf(f, "Btree: %d\n{\n", &ii);
    if (nr != 1 || i != ii) panic("load_preprocessed: bad idx value, i=%d, ii=%d", i, ii);
    nr = fscanf(f, " R: %d\n L: %d\n", &r, &l);
    if (nr != 2) panic("load_preprocessed: cannot read R/L value, r=%d, l=%d", r, l);
    if (r >= 0) btree[i].r = &btree[r];
    else btree[i].r = NULL;
    if (l >= 0) btree[i].l = &btree[l];
    else btree[i].l = NULL;
    btree[i].b = (Box*)malloc(sizeof(Box));
    nr = fscanf(f," Box:\n {\n  Triangle: %d\n  Triangle: %d\n", &r, &l);
    if (nr != 2) panic("load_preprocessed: cannot read t1/t2", HERE);
    if (r >= 0) btree[i].b->t1 = &ts[r];
    else btree[i].b->t1 = NULL;
    if (l >= 0) btree[i].b->t2 = &ts[l];
    else btree[i].b->t2 = NULL;
	/*printf("%d) %p %p\n", i, btree[i].b->t1, btree[i].b->t2);*/
    nr = fscanf(f,"  (%Lf,%Lf,%Lf)\n", &btree[i].b->minx, &btree[i].b->miny, &btree[i].b->minz);
    if (nr != 3) panic("load_preprocessed: cannot read minimal box vaules", HERE);
    nr = fscanf(f,"  (%Lf,%Lf,%Lf)\n", &btree[i].b->maxx, &btree[i].b->maxy, &btree[i].b->maxz);
    if (nr != 3) panic("load_preprocessed: cannot read maximal box vaules", HERE);
    fscanf(f, " }\n}\n");
    blist_add_box(&boxes, btree[i].b);
   }
 printf("\n");
 btree = &btree[root];
 fclose(f);
 return 1;
}


void save_preprocessed(FILE* f, int bin)
{
 void **p_idx;
 int n,i;
 BTree* tmp;
 g_idx = 0;
 count_indices();
 n = g_idx;
 p_idx = (void**)malloc(n*sizeof(void*));
 create_assoc_table(p_idx);
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

void initialize()
{
 ts = NULL;
 btree = NULL;
 boxes = NULL;
 nTriangles = 0;
 g_idx = 0;
}

/*void fixup(REAL* v)
{
 REAL Feps = 100;
 if ((fabs(*v) < Feps) && *v >= 0.) *v = -Feps;
 if ((fabs(*v) < Feps) && *v <= 0.) *v = Feps;
}*/

void make_order(Box* b)
{
 REAL temp;
 if (b->minx >= b->maxx) { temp = b->minx; b->minx = b->maxx; b->maxx = temp; }
 if (b->miny >= b->maxy) { temp = b->miny; b->miny = b->maxy; b->maxy = temp; }
 if (b->minz >= b->maxz) { temp = b->minz; b->minz = b->maxz; b->maxz = temp; }

 /*
 fixup(&(b->minx));
 fixup(&(b->miny));
 fixup(&(b->minz));
 fixup(&(b->maxx));
 fixup(&(b->maxy));
 fixup(&(b->maxz));
 */
}

void rotate_box(Box* b)
{
 REAL temp;
 /*printf("rot = %p\n", rot);*/
 switch(rot)
   {
	case X270:		
		temp = b->miny;
		b->miny = b->minz;
		b->minz = -temp;

		temp = b->maxy;
		b->maxy = b->maxz;
		b->maxz = -temp;

		make_order(b);
		break;

	case X180:
		b->miny = -b->miny;
		b->minz = -b->minz;

		b->maxy = -b->maxy;
		b->maxz = -b->maxz;

		make_order(b);
		break;

	case X90:
		temp = b->minz;
		b->minz = b->miny;
		b->miny = -temp;

		temp = b->maxz;
		b->maxz = b->maxy;
		b->maxy = -temp;

		make_order(b);
		break;

	case Y270:		
		temp = b->minz;
		b->minz = b->minx;
		b->minx = -temp;

		temp = b->maxz;
		b->maxz = b->maxx;
		b->maxx = -temp;

		make_order(b);
		break;

	case Y180:
		b->minz = -b->minz;
		b->minx = -b->minx;

		b->maxz = -b->maxz;
		b->maxx = -b->maxx;

		make_order(b);
		break;

	case Y90:
		temp = b->minx;
		b->minx = b->minz;
		b->minz = -temp;

		temp = b->maxx;
		b->maxx = b->maxz;
		b->maxz = -temp;

		make_order(b);
		break;

	case Z270:		
		temp = b->minx;
		b->minx = b->miny;
		b->miny = -temp;

		temp = b->maxx;
		b->maxx = b->maxy;
		b->maxy = -temp;

		make_order(b);
		break;

	case Z180:
		b->minx = -b->minx;
		b->miny = -b->miny;

		b->maxx = -b->maxx;
		b->maxy = -b->maxy;

		make_order(b);
		break;

	case Z90:
		temp = b->miny;
		b->miny = b->minx;
		b->minx = -temp;

		temp = b->maxy;
		b->maxy = b->maxx;
		b->maxx = -temp;

		make_order(b);
		break;

	default: return;
   }
}

void transform_box(Box* b)
{
 b->minx *= gsx;
 b->miny *= gsy;
 b->minz *= gsz;
 b->maxx *= gsx;
 b->maxy *= gsy;
 b->maxz *= gsz;

 rotate_box(b);

 b->minx += gtx;
 b->miny += gty;
 b->minz += gtz;
 b->maxx += gtx;
 b->maxy += gty;
 b->maxz += gtz;
}

void transform_tree(BTree* root)
{
 g_idx ++;
 if (g_idx && (g_idx % 1024) == 0) { printf("."); fflush(stdout); }
 if (root->r) transform_tree(root->r);
 if (root->l) transform_tree(root->l);
 transform_box(root->b);
}

void calculate_rotation_type(char* sprot)
{
 int i;
 rot = R0;
 if (!sprot) return;
 for (i=0;i<strlen(sprot);i++) if (sprot[i] >= 'A' && sprot[i] <= 'Z') sprot[i] += 0x20;
 if (!strcmp(sprot, "x90"))       rot = X90;
 else if (!strcmp(sprot, "x180")) rot = X180;
 else if (!strcmp(sprot, "x270")) rot = X270;
 else if (!strcmp(sprot, "y90"))  rot = Y90;
 else if (!strcmp(sprot, "y180")) rot = Y180;
 else if (!strcmp(sprot, "y270")) rot = Y270;
 else if (!strcmp(sprot, "z90"))  rot = Z90;
 else if (!strcmp(sprot, "z180")) rot = Z180;
 else if (!strcmp(sprot, "z270")) rot = Z270;
}

/*void recalculate_box_indices(Box* b, int offset)
{
}

void recalculate_indices(BTree* root, int offset)
{
 if (root->l) recalculate_indices(root->l, offset);
 if (root->r) recalculate_indices(root->r, offset);
 recalculate_box_indices(root->b, offset);
}*/

void merge_trees(BTree* btree1, BTree* btree2, BList* boxes1, BList* boxes2, 
	Triangle* ts1, Triangle* ts2, int nT1, int nT2)
{
 int i;
 ts = (Triangle*)malloc((nT1+nT2)*sizeof(Triangle));
 for (i=0;i<nT1;i++) 
   {
    memcpy((void*)(&ts[i]), (void*)(&ts1[i]), sizeof(Triangle));
    ts[i].idx = i;
   }
 for (i=nT1;i<nT1+nT2;i++) 
   {
    memcpy((void*)(&ts[i]), (void*)(&ts2[i-nT1]), sizeof(Triangle));
    ts[i].idx = i;
   }
 for (i=0;i<nT2;i++)
   {
    ts2[i].idx = i+nT1;
   }
 nTriangles = nT1+nT2;
/* recalculate_indices(btree2, nT1);*/
 btree = (BTree*)malloc(sizeof(BTree));
 btree->l = btree1;
 btree->r = btree2;
 btree->b = (Box*)malloc(sizeof(Box*));
 btree->b->t1 = NULL;
 btree->b->t2 = NULL;
 btree->b->minx = (btree1->b->minx < btree2->b->minx) ? btree1->b->minx : btree2->b->minx;
 btree->b->miny = (btree1->b->miny < btree2->b->miny) ? btree1->b->miny : btree2->b->miny;
 btree->b->minz = (btree1->b->minz < btree2->b->minz) ? btree1->b->minz : btree2->b->minz;
 btree->b->maxx = (btree1->b->maxx > btree2->b->maxx) ? btree1->b->maxx : btree2->b->maxx;
 btree->b->maxy = (btree1->b->maxy > btree2->b->maxy) ? btree1->b->maxy : btree2->b->maxy;
 btree->b->maxz = (btree1->b->maxz > btree2->b->maxz) ? btree1->b->maxz : btree2->b->maxz;
}

void mergetree(char fmt, char* inf1, char* inf2, char* outf)
{
 FILE  *in1, *in2, *out;
 BTree *btree1, *btree2;
 BList *boxes1, *boxes2;
 Triangle *ts1, *ts2;
 int nT1, nT2;
 int iret,binary;
 in1 = fopen(inf1, "rb");
 if (!in1) panic("Cannot read from file: %s\n", inf1);
 in2 = fopen(inf2, "rb");
 if (!in2) { fclose(in1); panic("Cannot read from file: %s\n", inf2); }
 out = fopen(outf, "wb");
 if (!outf) { fclose(in1); fclose(in2); panic("Cannot write to file: %s\n", outf); }
 
 binary = (fmt == 'b') ? 1 : 0;

 initialize();
 
 printf("Reading BTree from: %s\n", inf1);
 iret = load_preprocessed(in1);
 if (!iret) panic("Reading BTree failed.");
 btree1 = btree;
 boxes1 = boxes;
 ts1    = ts;
 nT1    = nTriangles;

 initialize();
 
 printf("Reading BTree from: %s\n", inf2);
 iret = load_preprocessed(in2);
 if (!iret) panic("Reading BTree failed.");
 btree2 = btree;
 boxes2 = boxes;
 ts2    = ts;
 nT2    = nTriangles;

 printf("Merging BTrees...\n");
 merge_trees(btree1, btree2, boxes1, boxes2, ts1, ts2, nT1, nT2);
 printf("Done.\n");
 
 printf("Writing BTree to: %s\n", outf);
 save_preprocessed(out, binary);
 printf("\n");

}


void btreeconv(char fmt, char* inf, char* outf, REAL tx, REAL ty, REAL tz, REAL sx, REAL sy, REAL sz, char* sprot)
{
 FILE *in, *out;
 int binary;
 int iret;
 in = fopen(inf, "rb");
 if (!in) panic("Cannot read from file: %s\n", inf);
 out = fopen(outf, "wb");
 if (!outf) { fclose(in); panic("Cannot write to file: %s\n", outf); }

 binary = (fmt == 'b') ? 1 : 0;
 initialize();

 printf("Reading BTree from: %s\n", inf);
 iret = load_preprocessed(in);
 if (!iret) panic("Reading BTree failed.");

 printf("Transforming BTree: T: (%3.2Lf,%3.2Lf,%3.2Lf), S: (%2.3Lf,%2.3Lf,%2.3Lf), R: %s\n", tx, ty, tz, sx, sy, sz, sprot);
 calculate_rotation_type(sprot);
 gtx = tx;
 gty = ty;
 gtz = tz;
 gsx = sx;
 gsy = sy;
 gsz = sz;
 g_idx = 0;
 transform_tree(btree);
 printf("\n");

 printf("Writing BTree to: %s\n", outf);
 save_preprocessed(out, binary);
 printf("\n");


 fclose(in);
 fclose(out);
}

void help()
{
 printf("\n");
 printf("arguments: [b|t] infile.btree outfile.btree tx ty tz sx sy sz rot\n");
 printf("b means binary output, t means text output.\n");
 printf("although rotations of AABB tree are imposible\nsome special cases can be implemented\n");
 printf("Set as rot parameter: X90, X180, X270, Y90, Y180, Y270, Z90, Z180, Z270 or R0\n");
 printf("R0 means no rotation, other (X,Y,Z)(90,180,270) are some special cases\n");
 printf("\nprogram can be used to merge trees:\n");
 printf("arguments: mt|mb infile1.btree infile2.btree outfile.btree\n");
 printf("\nIT SHOULD BE NOTICED that merging two trees don't gives You\n");
 printf("optimal BTree unless merged trees don't penetrate each other\n");
 printf("if they do so, generated tree is just a root gateway to choose\n");
 printf("which way to go and in this situation it will probably go both\n");
 printf("ALSO NOTE that merging many separate object this way is a very\n");
 printf("fast way to generate optimal trees for big scenes, but You should\n");
 printf("use some special order while generating structures, for example\n");
 printf("such usage is not good, because it creates unbalanced tree\n");
 printf("merge(c=a+b) merge (e=c+d) merge(g=e+f)\n");
 printf("INSTEAD You should create tree-like object structure, example:\n");
 printf("merge(c=a+b) merge(f=d+e) merge(g=c+f)\n");
 printf("first has 3 height while second has 2 height\n\n");
}

int main(int lb, char** par)
{
 if (lb >= 5 && par[1][0] == 'm' && (strlen(par[1]) >= 2))
   {
    mergetree(par[1][1], par[2], par[3], par[4]);
    return 0;
   }
 if (lb < 11) 
   { 
    help(); 
    return 1; 
   }
 btreeconv(par[1][0], par[2], par[3], atof(par[4]), atof(par[5]), atof(par[6]), atof(par[7]), atof(par[8]), atof(par[9]), par[10]);
 return 0;
}
