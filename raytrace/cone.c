/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI 3.1415926F	/* math constant */


void normalize(double* x, double* y, double* z)
{
 double len;
 double xx,yy,zz;
 xx = *x;
 yy = *y;
 zz = *z;
 len = sqrt(xx*xx+yy*yy+zz*zz);
 xx /= len;
 yy /= len;
 zz /= len;
 *x = xx;
 *y = yy;
 *z = zz;
}
	/* generates torrus grid with: */
	/* big radius R, small radius r */
	/* and big circle points# N */
	/* small circle points# n */
	/* finally points generated: N*n */
	/* generated output should be used */
	/* as input file for torrusrend */

void conegen(double r, double h, int n, int idx, double cr, double cg, double cb,
	double sr, double sg, double sb, double tr, double tg, double tb, int tid, int faces, 
	double sf, double tfr, double tfg, double tfb)
{
 int i,I;
 double x,y,z;
 double x1,y1,z1;
 double nx,ny,nz;
 double nx1,ny1,nz1;
 double nxh,nyh,nzh;
 double alfa;
 alfa = atan(h/r);
 for (i=0;i<n;i++)
   {
    I = i+1;
/*    if (I == n) I = 0;*/
    
    x  = r*cos(2.*PI*(double)i/(double)n);
    y  = 0.;
    z  = r*sin(2.*PI*(double)i/(double)n);
    
    x1 = r*cos(2.*PI*(double)I/(double)n);
    y1 = 0.;
    z1 = r*sin(2.*PI*(double)I/(double)n);

    nx = sin(alfa) * cos(2.*PI*(double)i/(double)n);
    ny = cos(alfa);
    nz = sin(alfa) * sin(2.*PI*(double)i/(double)n);
    
    nx1 = sin(alfa) * cos(2.*PI*(double)I/(double)n);
    ny1 = cos(alfa);
    nz1 = sin(alfa) * sin(2.*PI*(double)I/(double)n);
    
    nxh = sin(alfa) * cos(PI*(double)(i+I)/(double)n);
    nyh = cos(alfa);
    nzh = sin(alfa) * sin(PI*(double)(i+I)/(double)n);
    
    normalize(&nx, &ny, &nz);
    normalize(&nx1, &ny1, &nz1);
    normalize(&nxh, &nyh, &nzh);
       
       printf("Triangle: %d\n", idx);
       printf("{\n");
       printf("a: Vertex: (%f,%f,%f)\n",0.,0.,0.);
       printf("b: Vertex: (%f,%f,%f)\n",x,y,z);
       printf("c: Vertex: (%f,%f,%f)\n",x1,y1,z1);
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
       printf("a: Vertex: (%f,%f,%f)\n",0.,h,0.);
       printf("b: Vertex: (%f,%f,%f)\n",x1,y1,z1);
       printf("c: Vertex: (%f,%f,%f)\n",x,y,z);
       printf("texA: TexCoord: (%f,%f)\n", 0., 0.);
       printf("texB: TexCoord: (%f,%f)\n", (double)I/(double)n, 1.);
       printf("texC: TexCoord: (%f,%f)\n", (double)i/(double)n, 1.);
       printf("na: Vector: (%f,%f,%f)\n",nxh,nyh,nzh);
       printf("nb: Vector: (%f,%f,%f)\n",nx1,ny1,nz1);
       printf("nc: Vector: (%f,%f,%f)\n",nx,ny,nz);
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
 else conegen(atof(par[1]), atof(par[2]), atoi(par[3]), atoi(par[4]),
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
