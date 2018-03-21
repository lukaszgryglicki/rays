/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI 3.14159265F	/* math constant */


void cylindergen(double r, double h, int n, int idx, double cr, double cg, double cb,
	double sr, double sg, double sb, double tr, double tg, double tb, int tid, int faces, 
	double sf, double tfr, double tfg, double tfb)
{
 int i,I;
 double x,z;
 double x1,z1;
 double nx,nz;
 double nx1,nz1;
 for (i=0;i<n;i++)
   {
    I = i+1;
    
    x  = r*cos(2.*PI*(double)i/(double)n);
    z  = r*sin(2.*PI*(double)i/(double)n);
    
    x1 = r*cos(2.*PI*(double)I/(double)n);
    z1 = r*sin(2.*PI*(double)I/(double)n);

    nx = cos(2.*PI*(double)i/(double)n);
    nz = sin(2.*PI*(double)i/(double)n);
    
    nx1 = cos(2.*PI*(double)I/(double)n);
    nz1 = sin(2.*PI*(double)I/(double)n);
    
       printf("Triangle: %d\n", idx);
       printf("{\n");
       printf("a: Vertex: (%f,%f,%f)\n",0.,0.,0.);
       printf("b: Vertex: (%f,%f,%f)\n",x,0.,z);
       printf("c: Vertex: (%f,%f,%f)\n",x1,0.,z1);
       printf("texA: TexCoord: (%f,%f)\n", 0., 0.);
       printf("texB: TexCoord: (%f,%f)\n", 1., (double)i/(double)n);
       printf("texC: TexCoord: (%f,%f)\n", 1., (double)I/(double)n);
       printf("na: Vector: (%f,%f,%f)\n", 0., -1., 0.);
       printf("nb: Vector: (%f,%f,%f)\n", 0., -1., 0.);
       printf("nc: Vector: (%f,%f,%f)\n", 0., -1., 0.);
       printf("transparency: RGB: (%f,%f,%f)\n", tr, tg, tb);
       printf("specular: RGB: (%f,%f,%f)\n", sr, sg, sb);
       printf("diffuse: RGB: (%f,%f,%f)\n", cr, cg, cb);
       printf("surface: ABCD: (0,0,0,0)\n");
       printf("transparencyFactR: (1,%f)\n", tfr);
       printf("transparencyFactG: (1,%f)\n", tfg);
       printf("transparencyFactB: (1,%f)\n", tfb);
       printf("specularFact: %f\n", sf);
       printf("faces: %d\n",faces);
       printf("texture: %d\n", tid);
       printf("}\n");
       idx ++;
       printf("Triangle: %d\n", idx);
       printf("{\n");
       printf("a: Vertex: (%f,%f,%f)\n",x,h,z);
       printf("b: Vertex: (%f,%f,%f)\n",x1,0.,z1);
       printf("c: Vertex: (%f,%f,%f)\n",x,0.,z);
       printf("texA: TexCoord: (%f,%f)\n", 1., (double)i/(double)n);
       printf("texB: TexCoord: (%f,%f)\n", 0., (double)I/(double)n);
       printf("texC: TexCoord: (%f,%f)\n", 0., (double)i/(double)n);
       printf("na: Vector: (%f,%f,%f)\n",nx,0.,nz);
       printf("nb: Vector: (%f,%f,%f)\n",nx1,0.,nz1);
       printf("nc: Vector: (%f,%f,%f)\n",nx,0.,nz);
       printf("transparency: RGB: (%f,%f,%f)\n", tr, tg, tb);
       printf("specular: RGB: (%f,%f,%f)\n", sr, sg, sb);
       printf("diffuse: RGB: (%f,%f,%f)\n", cr, cg, cb);
       printf("surface: ABCD: (0,0,0,0)\n");
       printf("transparencyFactR: (1,%f)\n", tfr);
       printf("transparencyFactG: (1,%f)\n", tfg);
       printf("transparencyFactB: (1,%f)\n", tfb);
       printf("specularFact: %f\n", sf);
       printf("faces: %d\n",faces);
       printf("texture: %d\n",tid);
       printf("}\n");
       idx++;
       printf("Triangle: %d\n", idx);
       printf("{\n");
       printf("a: Vertex: (%f,%f,%f)\n",0.,h,0.);
       printf("b: Vertex: (%f,%f,%f)\n",x1,h,z1);
       printf("c: Vertex: (%f,%f,%f)\n",x,h,z);
       printf("texA: TexCoord: (%f,%f)\n", 1., 1.);
       printf("texB: TexCoord: (%f,%f)\n", 0., (double)I/(double)n);
       printf("texC: TexCoord: (%f,%f)\n", 0., (double)i/(double)n);
       printf("na: Vector: (%f,%f,%f)\n", 0., 1., 0.);
       printf("nb: Vector: (%f,%f,%f)\n", 0., 1., 0.);
       printf("nc: Vector: (%f,%f,%f)\n", 0., 1., 0.);
       printf("transparency: RGB: (%f,%f,%f)\n", tr, tg, tb);
       printf("specular: RGB: (%f,%f,%f)\n", sr, sg, sb);
       printf("diffuse: RGB: (%f,%f,%f)\n", cr, cg, cb);
       printf("surface: ABCD: (0,0,0,0)\n");
       printf("transparencyFactR: (1,%f)\n", tfr);
       printf("transparencyFactG: (1,%f)\n", tfg);
       printf("transparencyFactB: (1,%f)\n", tfb);
       printf("specularFact: %f\n", sf);
       printf("faces: %d\n",faces);
       printf("texture: %d\n", tid);
       printf("}\n");
       idx ++;
       printf("Triangle: %d\n", idx);
       printf("{\n");
       printf("a: Vertex: (%f,%f,%f)\n",x1,h,z1);
       printf("b: Vertex: (%f,%f,%f)\n",x1,0.,z1);
       printf("c: Vertex: (%f,%f,%f)\n",x,h,z);
       printf("texA: TexCoord: (%f,%f)\n", 1., (double)I/(double)n);
       printf("texB: TexCoord: (%f,%f)\n", 0., (double)I/(double)n);
       printf("texC: TexCoord: (%f,%f)\n", 1., (double)i/(double)n);
       printf("na: Vector: (%f,%f,%f)\n",nx1,0.,nz1);
       printf("nb: Vector: (%f,%f,%f)\n",nx1,0.,nz1);
       printf("nc: Vector: (%f,%f,%f)\n",nx,0.,nz);
       printf("transparency: RGB: (%f,%f,%f)\n", tr, tg, tb);
       printf("specular: RGB: (%f,%f,%f)\n", sr, sg, sb);
       printf("diffuse: RGB: (%f,%f,%f)\n", cr, cg, cb);
       printf("surface: ABCD: (0,0,0,0)\n");
       printf("transparencyFactR: (1,%f)\n", tfr);
       printf("transparencyFactG: (1,%f)\n", tfg);
       printf("transparencyFactB: (1,%f)\n", tfb);
       printf("specularFact: %f\n", sf);
       printf("faces: %d\n",faces);
       printf("texture: %d\n",tid);
       printf("}\n");
       idx++;
      }
    printf("\n");
}

int main(int lb, char** par)
{
 if (lb != 20)
   {
    printf("usage: cone r h n idx cr cg cb sr sg sb tr tg tb tid faces spec_fac tfr tfg tfb >> scene.dat\n");
    return 1;
   }
 else cylindergen(atof(par[1]), atof(par[2]), atoi(par[3]), atoi(par[4]),
	 atof(par[5]), atof(par[6]), atof(par[7]), 
	 atof(par[8]), atof(par[9]), atof(par[10]), 
	 atof(par[11]), atof(par[12]), atof(par[13]), atoi(par[14]), 
	 atoi(par[15]), atof(par[16]), atof(par[17]), atof(par[18]), 
	 atof(par[19]));
 return 0;
}

/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
