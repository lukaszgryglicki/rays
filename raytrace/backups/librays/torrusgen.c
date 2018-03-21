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

void torrusgen(double R, double r, int N, int n, int idx, double cr, double cg, double cb,
	double sr, double sg, double sb, double tr, double tg, double tb, int tid, int faces, 
	double sf, double tfr, double tfg, double tfb)
{
 int i,j;
 int I,J;
 double x,y,z;
 double x01,x10,x11;
 double y01,y10,y11;
 double z01,z10,z11;
 double nx,ny,nz;
 double nx01,nx10,nx11;
 double ny01,ny10,ny11;
 double nz01,nz10,nz11;
/* printf("%d,%d\n", N,n);*/
 for (i=0;i<N;i++)	/* using torrus properties */
   {
    I = i+1;
    if (I == N) I = 0;
    for (j=0;j<n;j++)
      {
       J = j+1;
       if (J == n) J = 0;
       x = R*cos((2.0*PI*(double)i)/(double)N)
	   +r*cos((2.0*PI*(double)i)/(double)N)
	   *cos((2.0*PI*(double)j)/(double)n);
       y = R*sin((2.0*PI*(double)i)/(double)N)
	   +r*sin((2.0*PI*(double)i)/(double)N)
	   *cos((2.0*PI*(double)j)/(double)n);
       z=  r*sin((2.0*PI*(double)j)/(double)n);
       
       x10 = R*cos((2.0*PI*(double)I)/(double)N)
	   +r*cos((2.0*PI*(double)I)/(double)N)
	   *cos((2.0*PI*(double)j)/(double)n);
       
       y10 = R*sin((2.0*PI*(double)I)/(double)N)
	   +r*sin((2.0*PI*(double)I)/(double)N)
	   *cos((2.0*PI*(double)j)/(double)n);
       
       z10=  r*sin((2.0*PI*(double)j)/(double)n);
       
       x01 = R*cos((2.0*PI*(double)i)/(double)N)
	   +r*cos((2.0*PI*(double)i)/(double)N)
	   *cos((2.0*PI*(double)J)/(double)n);
       
       y01 = R*sin((2.0*PI*(double)i)/(double)N)
	   +r*sin((2.0*PI*(double)i)/(double)N)
	   *cos((2.0*PI*(double)J)/(double)n);
       
       z01=  r*sin((2.0*PI*(double)J)/(double)n);
       
       x11 = R*cos((2.0*PI*(double)I)/(double)N)
	   +r*cos((2.0*PI*(double)I)/(double)N)
	   *cos((2.0*PI*(double)J)/(double)n);
       
       y11 = R*sin((2.0*PI*(double)I)/(double)N)
	   +r*sin((2.0*PI*(double)I)/(double)N)
	   *cos((2.0*PI*(double)J)/(double)n);
       
       z11=  r*sin((2.0*PI*(double)J)/(double)n);


       nx = cos(2.0*PI*(double)i/(double)N)*cos(2.0*PI*(double)j/(double)n);
       ny = sin(2.0*PI*(double)i/(double)N)*cos(2.0*PI*(double)j/(double)n);
       nz = sin(2.0*PI*(double)j/(double)n);
       
       nx10 = cos(2.0*PI*(double)I/(double)N)*cos(2.0*PI*(double)j/(double)n);
       ny10 = sin(2.0*PI*(double)I/(double)N)*cos(2.0*PI*(double)j/(double)n);
       nz10 = sin(2.0*PI*(double)j/(double)n);
       
       nx01 = cos(2.0*PI*(double)i/(double)N)*cos(2.0*PI*(double)J/(double)n);
       ny01 = sin(2.0*PI*(double)i/(double)N)*cos(2.0*PI*(double)J/(double)n);
       nz01 = sin(2.0*PI*(double)J/(double)n);
       
       nx11 = cos(2.0*PI*(double)I/(double)N)*cos(2.0*PI*(double)J/(double)n);
       ny11 = sin(2.0*PI*(double)I/(double)N)*cos(2.0*PI*(double)J/(double)n);
       nz11 = sin(2.0*PI*(double)J/(double)n);

/*       printf("normal: (%f,%f,%f)\n", nx, ny, nz);*/

       normalize(&nx, &ny, &nz);
       normalize(&nx10, &ny10, &nz10);
       normalize(&nx01, &ny01, &nz01);
       normalize(&nx11, &ny11, &nz11);
       
       printf("Triangle: %d\n", idx);
       printf("{\n");
       printf("a: Vertex: (%f,%f,%f)\n",x,y,z);
       printf("b: Vertex: (%f,%f,%f)\n",x01,y01,z01);
       printf("c: Vertex: (%f,%f,%f)\n",x11,y11,z11);
       printf("texA: TexCoord: (%f,%f)\n", (double)i/(double)N, (double)j/(double)n);
       printf("texB: TexCoord: (%f,%f)\n", (double)i/(double)N, (double)(j+1)/(double)n);
       printf("texC: TexCoord: (%f,%f)\n", (double)(i+1)/(double)N, (double)(j+1)/(double)n);
/*       printf("(i,j) = (%d,%d)\n", i,j);*/
       printf("na: Vector: (%f,%f,%f)\n",nx,ny,nz);
       printf("nb: Vector: (%f,%f,%f)\n",nx01,ny01,nz01);
       printf("nc: Vector: (%f,%f,%f)\n",nx11,ny11,nz11);
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
       printf("a: Vertex: (%f,%f,%f)\n",x,y,z);
       printf("b: Vertex: (%f,%f,%f)\n",x11,y11,z11);
       printf("c: Vertex: (%f,%f,%f)\n",x10,y10,z10);
       printf("texA: TexCoord: (%f,%f)\n", (double)i/(double)N, (double)j/(double)n);
       printf("texB: TexCoord: (%f,%f)\n", (double)(i+1)/(double)N, (double)(j+1)/(double)n);
       printf("texC: TexCoord: (%f,%f)\n", (double)(i+1)/(double)N, (double)(j)/(double)n);
       printf("na: Vector: (%f,%f,%f)\n",nx,ny,nz);
       printf("nb: Vector: (%f,%f,%f)\n",nx11,ny11,nz11);
       printf("nc: Vector: (%f,%f,%f)\n",nx10,ny10,nz10);
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
}

	/* You shoul pass 4 arguments to program: */
	/* R a double value - big radius */
	/* r a double value - small radius */
	/* N a integer value - how many points in big circle */
	/* n a integer value - how many points in small circle */
	/* the will be N*n generated points */

int main(int lb, char** par)
{
 if (lb != 21)
   {
    printf("usage: torrusgen R r N n idx cr cg cb sr sg sb tr tg tb tid faces spec_faci tfr tfg tfb >> scene.dat\n");
    printf("then:  torrusrend file.txt\n");
    printf("suggested values: R: 6-12, r: 1-4, N: 5-100, n: 5-100\n");
    return 1;
   }
 else torrusgen(atof(par[1]), atof(par[2]), atoi(par[3]), atoi(par[4]),
	 atoi(par[5]), atof(par[6]), atof(par[7]), atof(par[8]), 
	 atof(par[9]), atof(par[10]), atof(par[11]), 
	 atof(par[12]), atof(par[13]), atof(par[14]), atoi(par[15]), 
	 atoi(par[16]), atof(par[17]), atof(par[18]), atof(par[19]), 
	 atof(par[20]));
 return 0;
}

/* Written by MorgothDBMA, morgothdbma@o2.pl, tel: +48693582014 */
/* Lukasz Gryglicki MiNI M1 CC */
/* License BSD */
