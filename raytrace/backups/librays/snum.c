#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void display_num(double fp)
{
 char ndouble[sizeof(double)+1];
 char nldouble[sizeof(long double)+1];
 char nfloat[sizeof(float)+1];
 int lf, ld, lld;
 int i;
 unsigned char it;
 float fnum;
 double dnum;
 long double ldnum;
 lf  = sizeof(float);
 ld  = sizeof(double);
 lld = sizeof(long double);
 fnum = (float)fp;
 dnum = fp;
 ldnum = (long double)fp;
 memcpy((void*)nfloat, (void*)&fnum, lf);  
 memcpy((void*)ndouble, (void*)&dnum, ld);  
 memcpy((void*)nldouble, (void*)&ldnum, lld);  
 nfloat[lf] = 0;
 ndouble[ld] = 0;
 nldouble[lld] = 0;
 printf("floatH:\t\t ");
 for (i=0;i<lf;i++) 
    {
     if (i > 0 && !(i % 2)) printf(" ");
     it = nfloat[i];
     it %= 0x100;
     printf("%02x", it);
    }
 printf("\n");
 printf("doubleH:\t ");
 for (i=0;i<ld;i++) 
    {
     if (i > 0 && !(i % 2)) printf(" ");
     it = ndouble[i];
     it %= 0x100;
     printf("%02x", it);
    }
 printf("\n");
 printf("long doubleH:\t ");
 for (i=0;i<lld;i++) 
    {
     if (i > 0 && !(i % 2)) printf(" ");
     it %= 0x100;
     it = nldouble[i];
     printf("%02x", it);
    }
 printf("\n");
 for (i=0;i<lf;i++) if (nfloat[i] == 0) nfloat[i] = 1;
 for (i=0;i<ld;i++) if (ndouble[i] == 0) ndouble[i] = 1;
 for (i=0;i<lld;i++) if (nldouble[i] == 0) nldouble[i] = 1;
 printf("float:\t\t '%s'\n", nfloat);
 printf("double:\t\t '%s'\n", ndouble);
 printf("long double:\t '%s'\n", nldouble);
}

int main(int lb, char** par)
{
 double fp;
 if (lb < 2) fp = 0.;
 else fp = atof(par[1]);
 display_num(fp);
 return 0;
}

