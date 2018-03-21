#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double d_rand()
{
 return (double)rand() / (double)RAND_MAX;
}

int i_rand(int i)
{
 return rand() % i;
}

void rtriangle(int num)
{
 double dx,dy,dz;
 int i;
 int ntex;
 double bd,sd;
 ntex = 18;
printf(" Screen: (800,600)\n");
/*printf(" NormalDistorber: %f%%\n",d_rand());*/
printf(" MinShadow: %f\n", d_rand()/2.5);
printf(" Ambient: %f\n", 0.2+d_rand()/2.5);
printf(" MaxRecurse: %d\n", 5 + i_rand(5));
printf(" Backup: 200\n");
printf(" Observer: Vertex: (%f,%f,%f)\n", 0., 0., -600.*d_rand());
printf(" LookZ: %f%%\n", 30.+20.*d_rand());	
printf(" Light: Vertex: (%f,%f,%f)\n", 200.*d_rand()-100., 200.*d_rand()-100., -400.*d_rand());	
/*printf(" WorldTransform:\n");
printf(" }\n");
printf(" {\n");*/
printf(" TexDirectory: textures\n");
printf(" NumTextures: %d\n",ntex);
printf(" nTriangles: %d\n",num);
printf(" ListTransform: [0,%d]\n",num-1);
printf(" {\n");
printf("\n");
printf(" } \n");
bd = 400.;
sd = 100.;
for (i=0;i<num;i++)
{
 printf(" Triangle: %d\n",i);
 printf("  {\n");
 dx = d_rand()* bd - bd/2.;
 dy = d_rand()* bd - bd/2.;
 dz = d_rand()* bd - bd/2.;
 printf("   a: Vertex: (%f,%f,%f)\n", dx+d_rand()*sd-sd/2.,dy+d_rand()*sd-sd/2.,dz+d_rand()*sd-sd/2.);
 printf("   b: Vertex: (%f,%f,%f)\n", dx+d_rand()*sd-sd/2.,dy+d_rand()*sd-sd/2.,dz+d_rand()*sd-sd/2.);
 printf("   c: Vertex: (%f,%f,%f)\n", dx+d_rand()*sd-sd/2.,dy+d_rand()*sd-sd/2.,dz+d_rand()*sd-sd/2.);
 printf("   texA: TexCoord: (0,0)\n");
 printf("   texB: TexCoord: (0,1)\n");
 printf("   texC: TexCoord: (1,1)\n");
 printf("   na: Vector: (0,0,0)\n");
 printf("   nb: Vector: (0,0,0)\n");
 printf("   nc: Vector: (0,0,0)\n");
 printf("   transparencyA: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   specularA: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   diffuseA: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   transparencyB: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   specularB: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   diffuseB: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   transparencyC: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   specularC: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   diffuseC: RGB: (%f,%f,%f)\n",d_rand(),d_rand(),d_rand());
 printf("   transparencyFactR: (1,%f)\n", 1.+d_rand());
 printf("   transparencyFactG: (1,%f)\n", 1.+d_rand());
 printf("   transparencyFactB: (1,%f)\n", 1.+d_rand());
 printf("   normalDist: %f%%\n",d_rand());	
 printf("   specularFact: %f\n", d_rand()*200.);
 printf("   faces: 1\n");		
 printf("   texture: %d\n", i_rand(ntex+1));
 /*printf("   Transform:\n");
 printf("    {\n");
 printf("    }\n");*/
 printf("  }\n");
}
}

int main(int lb, char** par)
{
 time_t tm;
 time(&tm);
 srand((int)tm);
 if (lb != 2) { printf("%s numTriangles\n", par[0]); return 1; }
 else rtriangle(atoi(par[1]));
 return 0;
}
