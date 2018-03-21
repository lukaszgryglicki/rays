#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int tetragen()
{
 double ax, ay, az;
 double bx, by, bz;
 double cx, cy, cz;
 double dx, dy, dz;

 ax = -0.5;
 ay = -sqrt(6.0) / 9.0;
 az = -sqrt(3.0) / 6.0; 

 bx = 0.5;
 by = -sqrt(6.0) / 9.0;
 bz = -sqrt(3.0) / 6.0;

 cx = 0.0;
 cy = -sqrt(6.0) / 9.0;
 cz = sqrt(3.0) / 3.0;

 dx = 0.0;
 dy = sqrt(6.0) / 4.5;
 dz = 0.0;

 printf(
"Triangle: 0\n"
"{\n"
" a: Vertex: (%f,%f,%f)\n"
" b: Vertex: (%f,%f,%f)\n"
" c: Vertex: (%f,%f,%f)\n"
" texA: TexCoord: (0,0)\n"
" texB: TexCoord: (0,1)\n"
" texC: TexCoord: (1,0.5)\n"
" na: Vector: (0,0,0)\n"
" nb: Vector: (0,0,0)\n"
" nc: Vector: (0,0,0)\n"
" transparency: RGB: (0,0,0)\n"
" specular: RGB: (0,0,0)\n"
" diffuse: RGB: (1,1,1)\n"
" transparencyFactR: (1,1.220000)\n"
" transparencyFactG: (1,1.230000)\n"
" transparencyFactB: (1,1.250000)\n"
" specularFact: 25.000000\n"
" faces: 1\n"
" texture: 0\n"
"}\n"
,ax, ay, az, bx, by, bz, cx, cy, cz);


 printf(
"Triangle: 1\n"
"{\n"
" a: Vertex: (%f,%f,%f)\n"
" b: Vertex: (%f,%f,%f)\n"
" c: Vertex: (%f,%f,%f)\n"
" texA: TexCoord: (0,0)\n"
" texB: TexCoord: (0,1)\n"
" texC: TexCoord: (1,0.5)\n"
" na: Vector: (0,0,0)\n"
" nb: Vector: (0,0,0)\n"
" nc: Vector: (0,0,0)\n"
" transparency: RGB: (0,0,0)\n"
" specular: RGB: (0,0,0)\n"
" diffuse: RGB: (1,1,1)\n"
" transparencyFactR: (1,1.220000)\n"
" transparencyFactG: (1,1.230000)\n"
" transparencyFactB: (1,1.250000)\n"
" specularFact: 25.000000\n"
" faces: 1\n"
" texture: 0\n"
"}\n"
,ax, ay, az, bx, by, bz, dx, dy, dz);


 printf(
"Triangle: 2\n"
"{\n"
" a: Vertex: (%f,%f,%f)\n"
" b: Vertex: (%f,%f,%f)\n"
" c: Vertex: (%f,%f,%f)\n"
" texA: TexCoord: (0,0)\n"
" texB: TexCoord: (0,1)\n"
" texC: TexCoord: (1,0.5)\n"
" na: Vector: (0,0,0)\n"
" nb: Vector: (0,0,0)\n"
" nc: Vector: (0,0,0)\n"
" transparency: RGB: (0,0,0)\n"
" specular: RGB: (0,0,0)\n"
" diffuse: RGB: (1,1,1)\n"
" transparencyFactR: (1,1.220000)\n"
" transparencyFactG: (1,1.230000)\n"
" transparencyFactB: (1,1.250000)\n"
" specularFact: 25.000000\n"
" faces: 1\n"
" texture: 0\n"
"}\n"
,bx, by, bz, cx, cy, cz, dx, dy, dz);


 printf(
"Triangle: 3\n"
"{\n"
" a: Vertex: (%f,%f,%f)\n"
" b: Vertex: (%f,%f,%f)\n"
" c: Vertex: (%f,%f,%f)\n"
" texA: TexCoord: (0,0)\n"
" texB: TexCoord: (0,1)\n"
" texC: TexCoord: (1,0.5)\n"
" na: Vector: (0,0,0)\n"
" nb: Vector: (0,0,0)\n"
" nc: Vector: (0,0,0)\n"
" transparency: RGB: (0,0,0)\n"
" specular: RGB: (0,0,0)\n"
" diffuse: RGB: (1,1,1)\n"
" transparencyFactR: (1,1.220000)\n"
" transparencyFactG: (1,1.230000)\n"
" transparencyFactB: (1,1.250000)\n"
" specularFact: 25.000000\n"
" faces: 1\n"
" texture: 0\n"
"}\n"
,cx, cy, cz, ax, ay, az, dx, dy, dz);

	return 0;
}

int main()
{
	return tetragen();	
}

