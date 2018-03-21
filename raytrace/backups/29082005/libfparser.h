#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

double fpar_f(double);
double fpar2_f(double, double);
double fpar3_f(double, double, double);
int fpar_function(char* tab);
void fpar_info();
void fpar_free();
int fpar_ok();
int fpar2_ok();
int fpar3_ok();

/* supported funnctions: sin,cos,tan,tg,ctan,ctg,exp,ln,log,sqrt,abs
 * supported functions: asin,acos,asinh,acosh,sinh,cosh,cbrt,ceil,sgn
 * supported functions: tanh,atan,atanh,tgh,atg,atgh
 * supported operators: +,-,*,/,^, unary -
 * blank characters are skipped */
