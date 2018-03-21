#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void open_files(FILE** iF, FILE** oF, char* i, char* o)
{
 *iF = *oF = NULL;
 *iF = fopen(i, "rb");
 if (!(*iF)) { fprintf(stderr, "Cannot read from %s\n", i); exit(1); }
 *oF = fopen(o, "wb");
 if (!(*oF)) { fprintf(stderr, "Cannot write to %s\n", o); exit(1); }
}

void d2u(char* f1, char* f2)
{
 int zn;
 FILE *in, *out;
 open_files(&in, &out, f1, f2);
 while ((zn = fgetc(in)) != EOF)
  {
   if (zn != '\r') fputc(zn, out);
  }
 fclose(in);
 fclose(out);
}

void u2d(char* f1, char* f2)
{
 int zn;
 FILE *in, *out;
 open_files(&in, &out, f1, f2);
 while ((zn = fgetc(in)) != EOF)
  {
   if (zn == '\n') fputc('\r', out);
   fputc(zn, out);
  }
 fclose(in);
 fclose(out);
}
int main(int lb, char** par)
{
 if (lb != 4) 
   { 
    fprintf(stderr,"need parameter: (d2u or u2d)\n"
		   "example: ./compress d2u f1.txt f2.txt\n");
    exit(1);
   }
 if (!strcmp(par[1], "d2u")) d2u(par[2], par[3]);
 else if (!strcmp(par[1], "u2d")) u2d(par[2], par[3]);
 else { fprintf(stderr,"bad parameter\n"); exit(2); }
 return 0;
}
