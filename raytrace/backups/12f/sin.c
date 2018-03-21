#include <math.h>
#include <stdio.h>

#define PI 3.14159265

int main(int lb, char** par)
{
 double d,r;
 if (lb < 2) { printf("%s degrees\n", par[0]); return 1; }
 sscanf(par[1], "%lf", &d);
 printf("degrees: %f\n", d);
 r = (PI*d)/180.;
 printf("radians: %f\n", r);
 printf("sine = %1.8f\n"  , sin(r));
 printf("cosine = %1.8f\n", cos(r));
 printf("tang = %1.8f\n"  , tan(r));
 printf("cotang = %1.8f\n", 1./tan(r));
 return 0;
}
