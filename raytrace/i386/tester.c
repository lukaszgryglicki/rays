#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define REAL long double
typedef struct _V
{
 REAL x,y,z;
} Vector;

int nearly_equal(REAL, REAL, REAL);
REAL compute_angle(Vector*, Vector*, Vector*);

REAL ccompute_angle(Vector* a, Vector* b, Vector* c)
{
 register REAL x1,x2,y1,y2,z1,z2,len;
 x1 = b->x - a->x;
 x2 = c->x - a->x;
 y1 = b->y - a->y;
 y2 = c->y - a->y;
 z1 = b->z - a->z;
 z2 = c->z - a->z;
 /*printf("x = %Lf, %Lf\n", x1, x2);
 printf("y = %Lf, %Lf\n", y1, y2);
 printf("z = %Lf, %Lf\n", z1, z2);*/
 len = sqrt(x1*x1+y1*y1+z1*z1);
 if (len < 1e-7) return 4.;
/* printf("len1 = %Lf\n", len);*/
 x1 /= len;
 y1 /= len;
 z1 /= len;
 len = sqrt(x2*x2+y2*y2+z2*z2);
 if (len < 1e-7) return 4.;
/* printf("len2 = %Lf\n", len);*/
 x2 /= len;
 y2 /= len;
 z2 /= len;
 /*printf("x = %Lf, %Lf\n", x1, x2);
 printf("y = %Lf, %Lf\n", y1, y2);
 printf("z = %Lf, %Lf\n", z1, z2);*/
 return (x1*x2+y1*y2+z1*z2);
}

int main()
{
 REAL rs;
 Vector a,b,c;
 long long times;
 time_t t_start, t;
 a.x = 1.;
 a.y = 2.;
 a.z = 3.;
 b.x = 4.;
 b.y = 6.;
 b.z = 8.;
 c.x = 9.;
 c.y = 12.;
 c.z = 15.;
 /*nearly_equal(0.0001, 0.0003, 0.1);
 nearly_equal(0.0001, 0.0003, 0.1);*/
 times = 0;
 time(&t_start);
 while (1)
   {
    rs = compute_angle(&a, &b, &c);
    times ++;
    if (!(times % 100000)) 
      {
       time(&t);
       if (t >= t_start + 10) goto fine;
      }
   }
fine:
 printf("made %lld (%lld/s)\n", times, times / (t-t_start));
 /*printf("Ars = %Lf\n", rs);
 rs = compute_angle(&a, &b, &c);
 printf("Crs = %Lf\n", rs);*/
 return 0;
}

