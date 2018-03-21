#include <stdio.h>
#include <stdlib.h>
void faces(int idx, double tr, double tg, double tb, double sr, double sg, double sb, 
	double cr, double cg, double cb, double fr, double fg, double fb, double z_raise,
	double sf, int tid, int fac, int lt)
{
 int i;
 /*double tr,tg,tb;
 double sr,sg,sb;
 double cr,cg,cb;
 double fr,fg,fb;
 double z_raise;
 double sf;
 int tid,fac;
 tr = tg = tb = 0.;
 sr = sg = sb = 0;
 cr = cg = cb = 1.;
 fr = fg = fb = 1.5;
 z_raise = 0.;
 sf = 150.;
 tid = 2;
 fac = 2;*/
 if (z_raise == 100.)
 {
  z_raise = 0.587785;
 }
 if (lt == 1 || lt == 2)
 {
 printf("ListTransform: [%d,%d]\n", idx, idx+59);
 printf("{\n");
 printf("#apply your own transforms here.\n");
 printf("#Scale: (50,50,50)\n");
 printf("#NegateN: #(to see from outside)\n");
 printf("}\n");
 printf("ListTransform: [%d,%d]\n",idx+30,idx+59);
 printf("{\n");
 printf(" Translate: (0,0,-2.227033)\n");
 printf(" RotateX: 180\n");
 printf("}\n");
 printf("ListTransform: [%d,%d]\n",idx,idx+29);
 printf("{\n");
 printf(" Translate: (0,0,2.227033)\n");
 printf("}\n");
 for (i=1;i<=5;i++)
   {
    printf("ListTransform: [%d,%d]\n",idx+i*5,idx+i*5+4);
    printf("{\n");
    printf(" RotateZ: %f\n", (i-1) * 72.);
    printf(" Translate: (0,-1.376382,0)\n");
    printf(" RotateX: 63.434949\n");
    printf(" Translate: (0,-1.376382,0)\n");
    printf(" RotateZ: 180\n");
    printf("}\n");
   }
 for (i=1;i<=5;i++)
   {
    printf("ListTransform: [%d,%d]\n",idx+30+i*5,idx+30+i*5+4);
    printf("{\n");
    printf(" RotateZ: %f\n", (i-1) * 72.);
    printf(" Translate: (0,-1.376382,0)\n");
    printf(" RotateX: 63.434949\n");
    printf(" Translate: (0,-1.376382,0)\n");
    printf(" RotateZ: 180\n");
    printf("}\n");
   }
 }
 if (lt == 0 || lt == 2)
 {
 for (i=idx;i<idx+5;i++)
   {
     printf("Triangle: %i\n", i);
     printf("{\n");
/*     printf(" a:  Vertex: (0,0,0.587785)\n");*/
     printf(" a: Vertex: (0,0,%f)\n", z_raise);
     printf(" b: Vertex: (1,-1.376382,0)\n");
     printf(" c: Vertex: (-1,-1.376382,0)\n");
     if ((i - idx) == 0)
       {
        printf(" texA: TexCoord: (0.5,0.428571)\n");
        printf(" texB: TexCoord: (0.785714,0.035319)\n");
        printf(" texC: TexCoord: (0.214286,0.035319)\n");
      }
     if ((i - idx) == 1)
       {
        printf("  texA: TexCoord: (0.5,0.428571)\n");
        printf("  texB: TexCoord: (0.962295,0.57878)\n");
        printf("  texC: TexCoord: (0.785714,0.035319)\n");
      }
     if ((i - idx) == 2)
       {
        printf("  texA: TexCoord: (0.5,0.428571)\n");
        printf("  texB: TexCoord: (0.5,0.914658)\n");
        printf("  texC: TexCoord: (0.962295,0.57878)\n");
      }
     if ((i - idx) == 3)
       {
        printf("  texA: TexCoord: (0.5,0.428571)\n");
        printf("  texB: TexCoord: (0.037705,0.57878)\n");
        printf("  texC: TexCoord: (0.5,0.914658)\n");
      }
     if ((i - idx) == 4)
       {
        printf("  texA: TexCoord: (0.5,0.428571)\n");
        printf("  texB: TexCoord: (0.214286,0.035319)\n");
        printf("  texC: TexCoord: (0.037705,0.57878)\n");
      }
     printf(" na: Vector: (0,0,0)\n");
     printf(" nb: Vector: (0,0,0)\n");
     printf(" nc: Vector: (0,0,0)\n");
     printf(" transparency: RGB: (%f,%f,%f)\n",tr,tg,tb);
     printf(" specular: RGB: (%f,%f,%f)\n",sr,sg,sb);
     printf(" diffuse: RGB: (%f,%f,%f)\n",cr,cg,cb);
     printf(" transparencyFactR: (1,%f)\n",fr);
     printf(" transparencyFactG: (1,%f)\n",fg);
     printf(" transparencyFactB: (1,%f)\n",fb);
     printf(" specularFact: %f\n",sf);
     printf(" faces: %d\n",fac);
     printf(" texture: %d\n",tid);
     printf(" Transform:\n");
     printf(" {\n");
     printf(" RotateZ: %f\n", (i - idx) * 72.);
     printf(" }\n");
     printf("}\n");
   }
for (i=5;i<60;i++)
  {
   printf("CopyTriangle: %d<-%d\n", i + idx, (i % 5) + idx);
  }
 }
}

void help()
{
 printf("params: idx,t(rgb),s(rgb),d(rgb),tf(rgb)\n");
 printf("zraise: use 100 for special 60faces\n");
 printf("sf,tid,fac,lt(=1=>lt, =0=>tr, =2=>lt,tr\n");
}

int main(int lb, char** par)
{
 if (lb != 19) { help(); return 1; }
 faces(atoi(par[1]), atof(par[2]), atof(par[3]), atof(par[4]),	/* idx,tr,tg,tb */
	 atof(par[5]), atof(par[6]), atof(par[7]),	/* sr,sg,sb */
	 atof(par[8]), atof(par[9]), atof(par[10]),	/* cr,cg,cb */
	 atof(par[11]), atof(par[12]), atof(par[13]),	/* fr,fg,fb */
	 atof(par[14]), atof(par[15]),			/* zraise, sf */
	 atoi(par[16]), atoi(par[17]), atoi(par[18])    /* tid,fac,lt */
	 );
 return 0;
}
