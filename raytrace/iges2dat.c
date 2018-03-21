#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#define PI 3.1415926
#define ENTIRE_ARC 23040
#define eps 1e-3

int dbg;
int interval;
int interval_stub;
int ttid;
double tr,tg,tb,sr,sg,sb,cr,cg,cb, tfr, tfg, tfb, sf, nd;
double scx, scy, scz, rtx, rty, rtz;

struct Point
{
 double x,y,z;
};

typedef unsigned long ulong;

struct IGSObject
{
 int Did;		/* K1,K2 indeksy gorne sum */
 int K1, K2;		/* N1 = K1-M1+1, N2 = K2-M2+1; */
 int M1, M2;		/* stopnie krzywych bazowych */
 int A,B,C;		/* A = N1+2*M1, B = N2+2*M2 */
 			/* C = (K1+1)*(K2+1) */
 int N1,N2;		/* computed; */
 double* S;		/* S(-M1) ... S(M1+N1), S[A+1] */
 double* T;		/* T(-M2) ... T(M2+N2), T[B+1] */
 			/* FIXME: A dotyczy K2,M2,N2, B dotyczy K1,M1,N1 */
 double** W;		/* Wagi C elementow, W[K1+1][K2+1] */
 struct Point** P;	/* Punkty kontr. C elementow P[K1+1][K2+1].{x,y,z} */
 double U0, U1;		/* Przedz parametryzacji w/g X */
 double V0, V1;		/* Przedz parametryzacji w/g Y */
 			/* wyliczone trojkaty: */
 double* triangles;	/* jeden trojkat to 9 kolejnych double */
 int   ntriangles;	/* ilosc trojkatow */
};


static struct IGSObject* igs_array = NULL;
static int igs_num = 0;



void error(char* fmt, ...)
{
 va_list lst;
 va_start(lst,fmt);
 printf("CRITICAL ERROR: \t");
 vprintf(fmt,lst);
 printf("\n\tXCAM HALTED\n");
 fflush(stdout);
 va_end(lst);
 exit(1);
}


void debug(char* fmt, ...)
{
 va_list lst;
 if (dbg>=1)
   {
    va_start(lst,fmt);
    vprintf(fmt,lst);
    printf("\n");
    fflush(stdout);
    va_end(lst);
   }
}

void vdebug(char* fmt, ...)
{
 va_list lst;
 if (dbg>=2)
   {
    va_start(lst,fmt);
    vprintf(fmt,lst);
    printf("\n");
    fflush(stdout);
    va_end(lst);
   }
}


int readline(FILE* f, char* str)
{
 int zn;
 int i;
 vdebug("readline");
 i=0;
 while (1)
   {
    zn = fgetc(f);
    if (zn==EOF) return EOF;
    if (zn=='\n') { str[i] = 0; return 0; }
    str[i] = zn;
    i++;
   }
 return 0;
}


int param_id(char* line)
{
 char scan[82];
 int idx;
 vdebug("param_id"); 
 strcpy(scan, line);
 scan[72] = 0;
 idx = 71;
 while (scan[idx] != ' ' && idx>=0) idx--;
 return atoi(scan+idx);
}


int get_word_count(char* line)
{
 int i;
 int words = 0;
 vdebug("get_word_count"); 
 for (i=0;i<72;i++)
    if (line[i]==',') words++;
 vdebug("words = %d", words);
 return words;
}


void get_line(double* array, char* line, int idx, int words)
{
 double tmp;
 int nsc,i;
 int offset = 0;
 vdebug("get_line\n");
 line[72] = 0;
 line[73] = 0;
 for (i=0;i<words;i++)
   {
    nsc = sscanf(line+offset, "%lf", &tmp);
    vdebug("tmp = %f, offset=%d\n", tmp, offset);
/*    printf("%s", line+offset);*/
    if (nsc==1)
      {
       while (line[offset]!=',') offset++;
       offset++;
       array[idx+i] = tmp;
/*       printf("idx.\n");*/
      }
   }
}


void fill_args(double* array, char** par, int len, int maxidx)
{
 int i;
 int idx;
 int wc;
 idx = 0;
 debug("fill_args"); 
 for (i=0;i<len;i++)
   {
    wc = get_word_count(par[i]);
    get_line(array, par[i], idx, wc);
    idx += wc;
    if (idx>maxidx) error("max_index exceeded.");
   }
 if (idx!=maxidx) error("bad index count");
 for (i=0;i<idx;i++) vdebug("%f\n", array[i]);
}


double basis_bspline(double* T, double val, int dim, int off, int I, int J)
{
 double up;
 double down;
 double factor1;
 double factor2;
 factor1 = factor2 = 0.0;
 if (I>J) error("I more than J in BSPLINE basis");
 vdebug("ibasis_bspline: val=%f, [%f-%f], I=%d, J=%d, dim=%d\n", val,T[off+I], T[off+J], I, J, dim);
 if (dim == 0)
   {
    if (val >= T[I+off] && val <= T[J+off]) { vdebug("ret 1"); return 1.0; }
    else { vdebug("ret 0"); return 0.0; }
   }
 vdebug("recurse");
 up = (val-T[I+off]) * basis_bspline(T, val, dim-1, off, I, J-1);
 down = T[J+off-1] - T[I+off];
 if (down==0.0) factor1 = 0.0;
 else factor1 = up/down;
 up = (T[J+off]-val) * basis_bspline(T, val, dim-1, off, I+1, J);
 down = T[J+off] - T[I+off+1];
 if (down==0.0) factor2 += 0.0;
 else factor2 = up/down;
 return factor1+factor2;
}


double bspline_Integral_Z(double s, double t, int K1, int K2, int idx)
{
 int i,j;
 double up,down;
 double factor;
 int M1,M2;
 double b1,b2;
/* return 0.0;*/
 vdebug("bspline_Integral_Z\n");
 M1 = igs_array[idx].M1;
 M2 = igs_array[idx].M2;
 up = 0.0;
 down = 0.0;
 for (i=0;i<=K1;i++)
 {
 for (j=0;j<=K2;j++)
   {
    b1 = basis_bspline(igs_array[idx].S, s, M1, M1, i-M1, i+1);
    b2 = basis_bspline(igs_array[idx].T, t, M2, M2, j-M2, j+1);
    vdebug("Z> %f %f %f %f| %f %f %d %d %d\n", igs_array[idx].W[i][j], igs_array[idx].P[i][j].z, b1, b2,s,t,K1,K2,idx);
    factor = igs_array[idx].W[i][j]*igs_array[idx].P[i][j].z*b1*b2;
    up += factor;
    factor = igs_array[idx].W[i][j]*b1*b2;
    down += factor;
   }
 }
 if (down==0.0) return 0.0;
 vdebug("Z> %f\n", up/down);
 return up/down;
}


double bspline_Integral_Y(double s, double t, int K1, int K2, int idx)
{
 int i,j;
 double up,down;
 double factor;
 int M1,M2;
 double b1,b2;
 vdebug("bspline_Integral_Y\n");
/* return 0.0;*/
 M1 = igs_array[idx].M1;
 M2 = igs_array[idx].M2;
 up = 0.0;
 down = 0.0;
 for (i=0;i<=K1;i++)
 {
 for (j=0;j<=K2;j++)
   {
    b1 = basis_bspline(igs_array[idx].S, s, M1, M1, i-M1, i+1);
    b2 = basis_bspline(igs_array[idx].T, t, M2, M2, j-M2, j+1);
    vdebug("Y> %f %f %f %f| %f %f %d %d %d\n", igs_array[idx].W[i][j], igs_array[idx].P[i][j].z, b1, b2,s,t,K1,K2,idx);
    factor = igs_array[idx].W[i][j]*igs_array[idx].P[i][j].y*b1*b2;
    up += factor;
    factor = igs_array[idx].W[i][j]*b1*b2;
    down += factor;
   }
 }
 if (down==0.0) return 0.0;
 vdebug("Y> %f\n", up/down);
 return up/down;
}


double bspline_Integral_X(double s, double t, int K1, int K2, int idx)
{
 int i,j;
 double up,down;
 double factor;
 int M1,M2;
 double b1,b2;
 vdebug("bspline_Integral_X\n");
/* return 0.0;*/
 M1 = igs_array[idx].M1;
 M2 = igs_array[idx].M2;
 up = 0.0;
 down = 0.0;
 for (i=0;i<=K1;i++)
 {
 for (j=0;j<=K2;j++)
   {
    b1 = basis_bspline(igs_array[idx].S, s, M1, M1, i-M1, i+1);
    b2 = basis_bspline(igs_array[idx].T, t, M2, M2, j-M2, j+1);
    vdebug("X> %f %f %f %f| %f %f %d %d %d\n", igs_array[idx].W[i][j], igs_array[idx].P[i][j].z, b1, b2,s,t,K1,K2,idx);
    factor = igs_array[idx].W[i][j]*igs_array[idx].P[i][j].x*b1*b2;
    up += factor;
    factor = igs_array[idx].W[i][j]*b1*b2;
    down += factor;
   }
 }
 if (down==0.0) return 0.0;
 vdebug("X> %f\n", up/down);
 return up/down;
}


void generate_structures(double* args, int nargs, int idx)
{
 int M1,M2,K1,K2,N1,N2,A,B,C,i,j,x;
 double U0,U1,V0,V1,xm,ym,stepx,stepy;
 double t1,t2,t3;
 K1 = (int)args[2];
 K2 = (int)args[1];			/* FIXME points are inverted in my structures :-( */
 M1 = (int)args[3];			
 M2 = (int)args[4];
 N1 = (int)K1-M1+1;
 N2 = (int)K2-M2+1;
 A = N2+2*M1;				/* FIXME, SO A<->B also, this causes problems only... */
 B = N1+2*M2;
 C = (K1+1)*(K2+1);
 if (args[5]!=0 || args[6]!=0 || args[7]!=1 || args[8]!=0 || args[9]!=0) printf("PROP invalid for load_128\n");
 debug("generate_structures> K1=%d, K2=%d, M1=%d, M2=%d, N1=%d, N2=%d, A=%d, B=%d, C=%d\n",K1,K2,M1,M2,N1,N2,A,B,C);
 igs_array[idx].M1 = M1;
 igs_array[idx].M2 = M2;
 igs_array[idx].K1 = K1;
 igs_array[idx].K2 = K2;
 igs_array[idx].N1 = N1;
 igs_array[idx].N2 = N2;
 igs_array[idx].A  = A;
 igs_array[idx].B  = B;
 igs_array[idx].C  = C;
 igs_array[idx].S = (double*)malloc((A+2)*sizeof(double));
 if (!igs_array[idx].S) error("malloc S");
 for (i=0;i<=A;i++) igs_array[idx].S[i] = args[10+i];
 igs_array[idx].S[A+1] = igs_array[idx].S[A];
 igs_array[idx].T = (double*)malloc((B+2)*sizeof(double));
 if (!igs_array[idx].T) error("malloc T");
 for (i=0;i<=B;i++) igs_array[idx].T[i] = args[11+A+i];
 igs_array[idx].T[B+1] = igs_array[idx].T[B];
 igs_array[idx].W = (double**)malloc((K1+1)*sizeof(double*));
 if (!igs_array[idx].W) error("malloc W");
 for (i=0;i<=K1;i++)
   {
    igs_array[idx].W[i] = (double*)malloc((K2+1)*sizeof(double));
    if (!igs_array[idx].W[i]) error("malloc W[]");
   }
 for (i=0;i<=K1;i++)
    for (j=0;j<=K2;j++)
       igs_array[idx].W[i][j] = args[12+A+B+((K2+1)*i)+j];
 igs_array[idx].P = (struct Point**)malloc((K1+1)*sizeof(struct Point*));
 if (!igs_array[idx].P) error("malloc P");
 for (i=0;i<=K1;i++)
   {
    igs_array[idx].P[i] = (struct Point*)malloc((K2+1)*sizeof(struct Point));
    if (!igs_array[idx].P[i]) error("malloc P[]");
   }
 for (i=0;i<=K1;i++)
    for (j=0;j<=K2;j++)
      {
       igs_array[idx].P[i][j].x = args[12+A+B+C+3*(((K2+1)*i)+j)];
       igs_array[idx].P[i][j].y = args[12+A+B+C+3*(((K2+1)*i)+j)+1];
       igs_array[idx].P[i][j].z = args[12+A+B+C+3*(((K2+1)*i)+j)+2];
      }
 debug("book %d - prog %d\n", 12+A+B+4*C, 12+A+B+C+3*((K1+1)*(K2+1)));
 U0 = args[12+A+B+4*C];
 U1 = args[13+A+B+4*C];
 V0 = args[14+A+B+4*C];
 V1 = args[15+A+B+4*C];
 igs_array[idx].U0 = U0;
 igs_array[idx].U1 = U1;
 igs_array[idx].V0 = V0;
 igs_array[idx].V1 = V1;
 debug("S\n");
 for (i=0;i<=A;i++) debug("%f\n", igs_array[idx].S[i]);
 debug("T\n");
 for (i=0;i<=B;i++) debug("%f\n", igs_array[idx].T[i]);
 debug("W\n");
 for (i=0;i<=K1;i++)
 {
  for (j=0;j<=K2;j++) debug("%2.3f ", igs_array[idx].W[i][j]);
  debug("\n");
 }
 debug("P\n");
 for (i=0;i<=K1;i++)
   {
    for (j=0;j<=K2;j++)
      debug("(%2.3f,%3.2f,%3.2f) ",
    igs_array[idx].P[i][j].x,igs_array[idx].P[i][j].y,igs_array[idx].P[i][j].z);
    debug("\n");
   }
 debug("U=[%f,%f], V=[%f,%f]\n", U0,U1,V0,V1);
 igs_array[idx].ntriangles = 18*interval*interval;
 igs_array[idx].triangles = (double*)malloc(igs_array[idx].ntriangles*sizeof(double));
 if (!igs_array[idx].triangles) error("malloc failed while generating triangles");
 x = 0;
 stepx = (V1-V0)/(double)interval;
 stepy = (U1-U0)/(double)interval;
 for (xm=U0;xm<U1;xm+=stepx)
 for (ym=V0;ym<V1;ym+=stepy)
    {
     if (xm>U1-stepy-eps || ym>V1-stepy-eps) continue;
     if (xm<U0+eps || ym<V0+eps) continue;
     t1 = bspline_Integral_X(xm, ym, K1, K2, idx);
     t2 = bspline_Integral_Y(xm, ym, K1, K2, idx);
     t3 = bspline_Integral_Z(xm, ym, K1, K2, idx);
     debug("XYZ> %f,%f,%f",t1,t2,t3);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: %f,%f\n", xm,ym);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_Integral_X(xm+stepx, ym+stepy, K1, K2, idx);
     t2 = bspline_Integral_Y(xm+stepx, ym+stepy, K1, K2, idx);
     t3 = bspline_Integral_Z(xm+stepx, ym+stepy, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: +%f,+%f\n", xm+stepx,ym+stepy);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_Integral_X(xm+stepx, ym, K1, K2, idx);
     t2 = bspline_Integral_Y(xm+stepx, ym, K1, K2, idx);
     t3 = bspline_Integral_Z(xm+stepx, ym, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: +%f,%f\n", xm+stepx,ym);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_Integral_X(xm, ym, K1, K2, idx);
     t2 = bspline_Integral_Y(xm, ym, K1, K2, idx);
     t3 = bspline_Integral_Z(xm, ym, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: %f,%f\n", xm,ym);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_Integral_X(xm, ym+stepy, K1, K2, idx);
     t2 = bspline_Integral_Y(xm, ym+stepy, K1, K2, idx);
     t3 = bspline_Integral_Z(xm, ym+stepy, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: %f,+%f\n", xm,ym+stepy);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_Integral_X(xm+stepx, ym+stepy, K1, K2, idx);
     t2 = bspline_Integral_Y(xm+stepx, ym+stepy, K1, K2, idx);
     t3 = bspline_Integral_Z(xm+stepx, ym+stepy, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: +%f,+%f\n", xm+stepx,ym+stepy);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
    }
/* for (i=0;i<igs_array[idx].ntriangles;i++) if ((i%3)==2) igs_array[idx].triangles[i] = 0.0;*/
 debug("X = %d, ntr = %d\n", x, igs_array[idx].ntriangles);
 igs_array[idx].ntriangles = x;
}


void load_128(FILE* f, int n)
{
 char line[82];
 int llen,entry,i,j;
 int np=0;
 int nr=0;
 int pmid;
 int ping=1;
 int* lines;
 int* nlines;
 int* nargs;
 double** args;
 char*** param_data;
 debug("load_128");
 lines = (int*)malloc(igs_num*sizeof(int));
 if (!lines) error("malloc failed");
 for (i=0;i<igs_num;i++) lines[i] = 0;
 fseek(f, 0, SEEK_SET);
 while (readline(f, line)!=EOF)
   {
    llen = strlen(line);
    if (llen!=81 && llen!=80) error("BAD IGS file, linelen = %d", strlen(line));
    if (line[72] == 'D')
      {
       if (sscanf(line, "%d", &entry)<1) error("scanning entry error");
       if (entry==128)
         {
	  if (ping) igs_array[np].Did = atoi(line+73);
	  if (!ping) np++;
	  ping = !ping;
         }
      }
    if (line[72] == 'P')
    	{
	 pmid = param_id(line);
	 if (igs_array[nr].Did==pmid ||
	    (igs_array[nr].Did!=pmid && nr<igs_num-1 && igs_array[nr+1].Did==pmid))
	      {
	       if (igs_array[nr].Did!=pmid) nr++;
	       lines[nr]++;
	      }
	 else if (nr==igs_num) { debug("NOT intrested in rest DATA\n"); goto readprm; }
	}
   }
 readprm:
 for (i=0;i<igs_num;i++) debug("lines[%d] = %d\n", i, lines[i]);
 nlines = (int*)malloc(igs_num*sizeof(int));
 if (!nlines) error("malloc failed");
 for (i=0;i<igs_num;i++) nlines[i] = 0;
 param_data = (char***)malloc(igs_num*sizeof(char**));
 if (!param_data) error("malloc char***");
 for (i=0;i<igs_num;i++)
   {
    param_data[i] = (char**)malloc(lines[i]*sizeof(char*));
    if (!param_data[i]) error("malloc char**");
    for (j=0;j<lines[i];j++)
      {
       param_data[i][j] = (char*)malloc(82*sizeof(char));
       if (!param_data[i][j]) error("malloc char*");
      }
   }
 fseek(f, 0, SEEK_SET);
 nr = 0;
 while (readline(f, line)!=EOF)
   {
    if (line[72] == 'P')
    	{
	 pmid = param_id(line);
	 if (igs_array[nr].Did==pmid ||
	    (igs_array[nr].Did!=pmid && nr<igs_num-1 && igs_array[nr+1].Did==pmid))
	      {
	       if (igs_array[nr].Did!=pmid) nr++;
	       strcpy(param_data[nr][nlines[nr]], line);
	       nlines[nr]++;
	      }
	 else if (nr==igs_num) { debug("IGS READ successfully.\n"); goto rsdis; }
	}
   }
 rsdis:
 nargs = (int*)malloc(igs_num*sizeof(int));
 if (!nargs) error("malloc int*");
 for (i=0;i<igs_num;i++) nargs[i] = 0;
 for (i=0;i<igs_num;i++)
    for (j=0;j<lines[i];j++) nargs[i] += get_word_count(param_data[i][j]);
 for (i=0;i<igs_num;i++) debug("nargs[%d] = %d\n", i, nargs[i]);
 args = (double**)malloc(igs_num*sizeof(double*));
 if (!args) error("malloc double**");
 for (i=0;i<igs_num;i++)
   {
    args[i] = (double*)malloc(nargs[i]*sizeof(double));
    if (!args[i]) error("malloc int*");
   }
 for (i=0;i<igs_num;i++)
   {
    fill_args(args[i], param_data[i], lines[i], nargs[i]);
    generate_structures(args[i], nargs[i], i);
   }
 for (i=0;i<igs_num;i++)
  {
   for (j=0;j<lines[i];j++) free(param_data[i][j]);
   free(args[i]);
   free(param_data[i]);
  }
 free(args); args = NULL;
 free(nargs); nargs = NULL;
 free(lines); lines = NULL;
 free(nlines); nlines = NULL;
 free(param_data); param_data = NULL;
}


void load_igsfile(char* fpath)
{
 FILE* f;
 char line[256];	/* >= 82 bytes MIN */
 int n;
 int entry;
 int llen;
 debug("load_IGSfile: %s", fpath);
 f = fopen(fpath,"r");
 if (!f) error("cannot open file: %s", fpath);
 n=0;
 strcpy(line,"");
 debug("determining how much lines read");
 while (readline(f, line)!=EOF)
   {
    llen = strlen(line);
    if (llen!=81 && llen!=80) error("BAD IGS file, linelen = %d", strlen(line));
    if (line[72] == 'D')
      {
       if (sscanf(line, "%d", &entry)<1) error("scanning entry error");
       if (entry==128) n++;
      }
    if (line[72] == 'P') goto parameter;
   }
 error("No PARAMETER SECTION found.");
 parameter:
 debug("N128=%d\n", n);
 if (n<=0) error("no 128 entires found");
 if (n)
   {
    if (n%2) error("bad lines count");
    n /= 2;
    fseek(f, 0, SEEK_SET);
    igs_array = (struct IGSObject*)malloc(n*sizeof(struct IGSObject));
    if (!igs_array) error("malloc igs_array");
    igs_num = n;
    debug("SIZEOF IGSObject is %d\n", sizeof(struct IGSObject));
    load_128(f,n);
   }
 else { igs_array = 0; igs_num = 0; }
 debug("IGS file read.");
}

void free_igsobj(struct IGSObject* obj)
{
 int i;
 debug("free_IGS_obj");
 if (!obj) return;
 free(obj->triangles);
 free(obj->T);
 free(obj->S);
 for (i=0;i<=obj->K1;i++) { free(obj->P[i]); free(obj->W[i]); }
 free(obj->P);
 free(obj->W);
 obj = NULL;
}



void free_temporary_structures()
{
 int i;
 debug("free_temporary_structures");
 if (igs_num>0)
   {
    if (!igs_array) return;
    for (i=0;i<igs_num;i++) free_igsobj(&igs_array[i]);
    free(igs_array);
    igs_array = NULL;
   }
}

void write_header(FILE* f)
{
fprintf(f, "Screen: (1024,768)\n");	
fprintf(f, "Backup: 256\n");
fprintf(f, "Observer: Vertex: (0,0,-400)\n");
fprintf(f, "Light: Vertex: (-50,75,-100)\n");
fprintf(f, "TexDirectory: textures\n");
fprintf(f, "NumTextures: 33\n");
fprintf(f, "nTriangles: 0\n");
fprintf(f, "nNURBS: %d\n", igs_num);
}

void write_nurbs(struct IGSObject* obj, FILE* f, int num)
{
 int i,j;
 fprintf(f, "NURBS: %d\n {\n", num);
 fprintf(f, " DimU: %d\n", obj->M1);
 fprintf(f, " DimV: %d\n", obj->M2);
 fprintf(f, " nptsU: %d\n", obj->K1+1);
 fprintf(f, " nptsV: %d\n", obj->K2+1);
 fprintf(f, " divU: %d\n", interval_stub);
 fprintf(f, " divV: %d\n", interval_stub);
 fprintf(f, " ITex: %d\n", 0);
 fprintf(f, " Faces: %d\n", 1);
 fprintf(f, " Tid: %d\n", ttid);
 fprintf(f, " userNodes: %d\n", 0);
 fprintf(f, " userKnots: %d\n", 1);
 fprintf(f, " UKnots:");
 for (i=0;i<=obj->B;i++)
   {
    fprintf(f, " %f", obj->T[i]);
   }
 fprintf(f, "\n");
 fprintf(f, " VKnots:");
 for (i=0;i<=obj->A;i++)
   {
    fprintf(f, " %f", obj->S[i]);
   }
 fprintf(f, "\n");
 for (i=0;i<=obj->K1;i++)
   {     
    fprintf(f, "  Line=%d\n", i);
    for (j=0;j<=obj->K2;j++)
       fprintf(f, "   {v=(%f,%f,%f),t=(%f,%f,%f),s=(%f,%f,%f),c=(%f,%f,%f),f=(%f,%f,%f),sf=%f,nd=%f%%,tc=(0.000000,0.000000),w=%f}\n"
	       ,obj->P[i][j].x, obj->P[i][j].y, obj->P[i][j].z, tr, tg, tb, sr, sg, sb, cr, cg, cb, tfr, tfg, tfb, sf, nd, obj->W[i][j]);
   }
 fprintf(f, " }\n");
}

void write_transform(FILE* f)
{
 fprintf(f, "NURBSTransform: [0,%d]\n {\n Scale: (%f,%f,%f)\n Rotate: (%f,%f,%f)\n }\n\n"
	 , igs_num-1, scx, scy, scz, rtx, rty, rtz);
}

void save_datfile(char* datfile)
{
 FILE* f;
 int i;
 if (igs_num <0) { printf("NO NURBS read from input file.\n"); return; }
 f = fopen(datfile, "w");
 if (!f) { printf("cannot open: %s, exiting...\n", datfile); return; }
 write_header(f);
 for (i=0;i<igs_num;i++) write_nurbs(&igs_array[i], f, i);
 write_transform(f);
 fclose(f);
}


void process_IGES(char* igsfile, char* datfile)
{
 load_igsfile(igsfile);
 save_datfile(datfile);
 free_temporary_structures();
 free(igs_array);
 igs_array = NULL;
}


void help()
{
 printf("switches: -i input, -o output, -h (help), -d (debug), -t num (triang_fact)");
 printf(", -C \"n1 n2 n3\" (color 3val), -T (trans 3val), -S (spec 3val), -F (trans_fact 3val)");
 printf(", -N (ndist), -I (shiness) -X (scale 3val), -R (rotate 3val) -m (texID)\n");
}

void set_3val(double* v1, double* v2, double* v3, char* arg)
{
 int i;
/* printf("Setting 3 values from: \"%s\"\n", arg);*/
 i = sscanf(arg, "%lf %lf %lf", v1, v2, v3);
 if (i != 3) printf("set_3val: WARNING: cannot properly set 3 values.\n");
}


int main(int lb, char** par)
{
 char u;
 char igsfile[1024];
 char datfile[1024];
 dbg = 0;
 interval = 2;
 interval_stub = 12;			/* FIXME: stub for leaving triangulation */
 					/* to external application (rayslib) */
 strcpy(igsfile, "Iges/BSPLINE.igs");
 strcpy(datfile, "NurbsFromIges.dat");
 tr = tg = tb = 0.;
 sr = sg = sb = 0.;
 cr = cg = cb = 1.;
 scx = scy = scz = 1.;
 rtx = rty = rtz = 0.;
 tfr = tfb = tfg = 1.1;
 ttid = 2;
 sf = 20.;
 nd = 0.;
 while ((u = getopt(lb, par, "m:R:X:T:C:S:F:N:I:i:o:t:dh")) != -1)
 {
  switch (u)
   {
    case 'h': help(); return 1; break;
    case 'T': set_3val(&tr, &tg, &tb, optarg); break;
    case 'C': set_3val(&cr, &cg, &cb, optarg); break;
    case 'S': set_3val(&sr, &sg, &sb, optarg); break;
    case 'F': set_3val(&tfr, &tfg, &tfb, optarg); break;    
    case 'd': dbg = 1; break;
    case 'N': nd = atof(optarg); break;
    case 'I': sf = atof(optarg); break;
    case 'm': ttid = atoi(optarg); break;
    case 't': interval_stub = atoi(optarg); break;
    case 'i': strcpy(igsfile, optarg); break;
    case 'o': strcpy(datfile, optarg); break;
    case 'X': set_3val(&scx, &scy, &scz, optarg); break;
    case 'R': set_3val(&rtx, &rty, &rtz, optarg); break;
    default: printf("%s: Unrecognized option\n",par[0]); return 1;
   }
 }
 process_IGES(igsfile, datfile);
 return 0;
}

