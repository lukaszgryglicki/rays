#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define HERE __FILE__,__LINE__
#define REAL double

typedef struct _V
{
 REAL x,y,z;
 REAL nx,ny,nz;
 REAL tx,ty;
 REAL idx;
} Vertex;

typedef struct _T
{
 Vertex *a, *b, *c;
} Triangle;
int idx; 
REAL tr; 
REAL tg; 
REAL tb; 
REAL sr; 
REAL sg; 
REAL sb;
REAL cr; 
REAL cg; 
REAL cb; 
REAL fr; 
REAL fg; 
REAL fb; 
REAL sf;
int tid;
int fac; 
int ran;
int hdr;

void err(char* f, int l)
{
 printf("Error at %s: %d\n", f, l);
 exit(1);
}

REAL d_rand()
{
 return (REAL)rand() / (REAL)RAND_MAX;
}

int i_rand(int i)
{
 return rand() % i;
}

void write_triangles(Triangle* t, int n)
{
 int i;
 int ntex;
 ntex = 18;
 if (hdr)
  {
   printf("Screen: (800,600)\n");
/*printf(" NormalDistorber: %f%%\n",d_rand());*/ 
   printf("MinShadow: %f\n", d_rand()/2.5);
   printf("Ambient: %f\n", 0.2+d_rand()/2.5);
   printf("MaxRecurse: %d\n", 5 + i_rand(5));
   printf("Backup: 200\n");
   printf("Observer: Vertex: (%f,%f,%f)\n", 0., 0., -200.-150*d_rand());
   printf("LookZ: %f%%\n", 30.+20.*d_rand());	
   printf("Light: Vertex: (%f,%f,%f)\n", 200.*d_rand()-100., 200.*d_rand()-100., -400.*d_rand());	
/*printf(" WorldTransform:\n");
printf(" }\n");
printf(" {\n");*/
   printf("TexDirectory: textures\n");
   printf("NumTextures: %d\n",ntex);
   printf("nTriangles: %d\n",n);
  }
 printf("ListTransform: [%d,%d]\n",idx, idx+n-1);
 printf("{\n");
 printf("\n");
 printf("} \n");
 for (i=0;i<n;i++)
   {
    printf("Triangle: %d\n",i+idx);
    printf("{\n");
    printf(" a: Vertex: (%f,%f,%f)\n", t[i].a->x, t[i].a->y, t[i].a->z);
    printf(" b: Vertex: (%f,%f,%f)\n", t[i].b->x, t[i].b->y, t[i].b->z);
    printf(" c: Vertex: (%f,%f,%f)\n", t[i].c->x, t[i].c->y, t[i].c->z);
    printf(" texA: TexCoord: (%f,%f)\n", t[i].a->tx, t[i].a->ty);
    printf(" texB: TexCoord: (%f,%f)\n", t[i].b->tx, t[i].b->ty);
    printf(" texC: TexCoord: (%f,%f)\n", t[i].c->tx, t[i].c->ty);
    printf(" na: Vector: (%f,%f,%f)\n",t[i].a->nx, t[i].a->ny,t[i].a->nz);
    printf(" nb: Vector: (%f,%f,%f)\n",t[i].b->nx, t[i].b->ny,t[i].b->nz);
    printf(" nc: Vector: (%f,%f,%f)\n",t[i].c->nx, t[i].c->ny,t[i].c->nz);
    if (ran)
      {
       printf(" transparency: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
       printf(" specular: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
       printf(" diffuse: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
       printf(" transparencyFactR: (1,%f)\n", 1.+d_rand());
       printf(" transparencyFactG: (1,%f)\n", 1.+d_rand());
       printf(" transparencyFactB: (1,%f)\n", 1.+d_rand());
       printf(" normalDist: %f%%\n",d_rand());	
       printf(" specularFact: %f\n", d_rand()*200.);
       printf(" faces: %d\n",fac);		
       printf(" texture: %d\n", i_rand(ntex+1));
      }
    else
      {
       printf(" transparency: RGB: (%f,%f,%f)\n",tr,tg,tb);
       printf(" specular: RGB: (%f,%f,%f)\n",sr,sg,sb);
       printf(" diffuse: RGB: (%f,%f,%f)\n",cr,cg,cb);
       printf(" transparencyFactR: (1,%f)\n", fr);
       printf(" transparencyFactG: (1,%f)\n", fg);
       printf(" transparencyFactB: (1,%f)\n", fb);
       /*printf("   normalDist: %f%%\n",d_rand());*/
       printf(" specularFact: %f\n", sf);
       printf(" faces: %d\n",fac);		
       printf(" texture: %d\n", tid);
      }
 /*printf("   Transform:\n");
 printf("    {\n");
 printf("    }\n");*/
    printf("  }\n");
   }
}

void uli2dat()
{
 int nr;
 int n,i;
 int nt;
 int ai,bi,ci;
 Vertex* v;
 Triangle* t;
 nr = scanf("%d\n", &n);
 if (nr != 1) err(HERE);
 if (n < 0) err(HERE);
 v = (Vertex*)malloc(n*sizeof(Vertex));
 for (i=0;i<n;i++)
   {
    nr = scanf("%lf %lf %lf %lf %lf %lf %lf %lf\n", &v[i].x, &v[i].y, &v[i].z, &v[i].nx, &v[i].ny, &v[i].nz, &v[i].tx, &v[i].ty);
    if (nr != 8) err(HERE);
   }
 nr = scanf("%d\n", &nt);
 if (nr != 1) err(HERE);
 if (nt < 0) err(HERE);
 t = (Triangle*)malloc(nt*sizeof(Triangle));
 for (i=0;i<nt;i++)
   {
    nr = scanf("%d %d %d\n", &ai, &bi, &ci);
    if (nr != 3) err(HERE);
    if (ai < 1 || ai > n) err(HERE);
    if (bi < 1 || bi > n) err(HERE);
    if (ci < 1 || ci > n) err(HERE);
    t[i].a = &v[ai-1];
    t[i].b = &v[bi-1];
    t[i].c = &v[ci-1];
   }
 write_triangles(t, nt);
}

void help(char* par)
{
 printf("usage: %s params\n", par);
 printf("params: idx,t(rgb),s(rgb),d(rgb),tf(rgb)\n");
 printf("sf,tid,fac,rand,hdr\n");
}

int main(int lb, char** par)
{
 time_t tm;
 time(&tm);
 srand((int)tm);
 if (lb != 19) { help(par[0]); err(HERE); }
 idx =atoi(par[1]);
 tr = atof(par[2]); 
 tg = atof(par[3]); 
 tb = atof(par[4]); 
 sr = atof(par[5]);
 sg = atof(par[6]); 
 sb = atof(par[7]); 
 cr = atof(par[8]);
 cg = atof(par[9]); 
 cb = atof(par[10]); 
 fr = atof(par[11]);
 fg = atof(par[12]); 
 fb = atof(par[13]); 
 sf = atof(par[14]);
 tid =atoi(par[15]);
 fac =atoi(par[16]); 
 ran =atoi(par[17]); 
 hdr =atoi(par[18]); 
 uli2dat();
 return 0;
}
