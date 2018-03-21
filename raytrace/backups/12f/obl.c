#include <math.h>
#include <stdio.h>

#define PI 3.14159265

int main()
{
 double a;
 double kat;
 double b;
 double hb;
 double hp;
 double hc;
 double h;
 double sin_ahalf;
 double alfa;
 double h5;
 double hh;
 double kat2;
 a = 1.;
 printf("a = %f\n",a);
 kat = (36. * PI) / 180.;
 kat2 = (54. * PI) / 180.;
 b = 2. * a * cos(kat);
 printf("b = %f\n", b);
 h5 = sin(kat2) * a;
 printf("h5 = %f\n", h5);
 hh = sqrt(a*a - h5*h5);
 printf("hh = %f\n", hh);
 hb = sqrt(a*a - 0.25*b*b);
 printf("hb = %f\n", hb);
 hp = 0.5 * b * sqrt(3.);
 printf("hp = %f\n", hp);
 h = sqrt(hb*hb - (hp*hp)/9.);
 printf("h = %f\n", h);
 h = sqrt(a*a - (hp*hp)*(4./9.));
 printf("h = %f\n", h);
 hc = sin(kat) * b;
 printf("hc = %f\n", hc);
 sin_ahalf = (0.5*b) / hc;
 printf("sin(a/2) = %f\n", sin_ahalf);
 alfa = 2. *  asin(sin_ahalf);
 alfa = 180. - (alfa * 180.) / PI;
 printf("alfa = %f\n", alfa);
 return 0;
}

