#include <unistd.h>			//FIXME MAX DEBUG ALL
#include <stdio.h>			//SET MALLOC_PROTECT BELOW & ABOVE
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <math.h>
#include <pthread.h>
#include "xkeys.h"
#define MAX_LINE 0x400
#define MAXPATH  1024
#define ZMAGIC 0X6FFFFFFF
#define FAR_AWAY 1e7
#define FREE_MUTEX  pthread_mutex_unlock(&mutex);
#define WAIT_MUTEX  pthread_mutex_lock(&mutex);
#define PI 3.1415926
#define MILL_POINTS 222
#define MAX_LPOS 10.0
#define ENTIRE_ARC 23040
#define MIN_TIMER 80000
#define MAX_INTERVAL 0x200
#define eps 1e-3
#define __asm extern
static int   interval = 50;
static double OINT=0.2;
static double LX;
static double LY;
static double LZ;
static int START_AX = 0;
static int START_AY = 0;
static int START_AZ = 0;
static unsigned int line_r,line_g,line_b;
static unsigned int mill_r,mill_g,mill_b;
static unsigned int igs6_r,igs6_g,igs6_b;
static unsigned int igs_r,igs_g,igs_b;
static double ZOFFSET = 0.0;
static double YOFFSET = 0.0;
static double XOFFSET = 0.0;
static int AX_STEP = 15;
static int AY_STEP = 15;
static int AZ_STEP = 15;
static int TX_STEP = 270;
static int TY_STEP = 0;
static int TZ_STEP = 0;
static int MAXLINES = 1024;
static int TIMERINT = 90000;
static int LINESINC = 10;
static double SCALEX  = 1.41;
static double SCALEY = 1.41;
static double zoff = 12.0;
static int angle_x;
static int angle_y;
static int angle_z;
static int rotate_tool_x = 0;
static int rotate_tool_y = 0;
static int rotate_tool_z = 0;
static double sines[360];
static double cosines[360];
static int time_async;
static int Xup;
static int done;
static int stopped;
static int dbg=0;
static double nx;
static double ny;
static double nz;
static double lx;
static double ly;
static double lz;
static double lxb;
static double lyb;
static double lzb;
static double light_1;
static double light_2;
static double light_3;
static int want_triangles=0;
static int t1x,t2x,t3x,t1y,t2y,t3y;
static double ox=0.0;
static double oy=0.0;
static double oz=0.0;
static int ldiff = 10;
static int want_igs  = 1;
static int want_cam  = 1;
pthread_mutex_t mutex;
pthread_t timer_thr;
Display* dsp;
GC gc;
Window win;

__asm unsigned long RGB16(int,int,int);
__asm unsigned long RGB24(int,int,int);
__asm int RR16(unsigned long);
__asm int RG16(unsigned long);
__asm int RB16(unsigned long);
__asm int RR24(unsigned long);
__asm int RG24(unsigned long);
__asm int RB24(unsigned long);

struct Point
{
 double x,y,z;
};
struct Point* path;
struct Point* path_buff;
int npoints;
typedef unsigned long ulong;

struct Zbuff
{
 ulong col;
 double Z;
};
struct Zbuff** zbuff = NULL;
static double* mill = NULL;
static double* mill_buff = NULL;

struct IGSObject
{
 int Did;		// K1,K2 indeksy gorne sum
 int K1, K2;		// N1 = K1-M1+1, N2 = K2-M2+1;
 int M1, M2;		//stopnie krzywych bazowych
 int A,B,C;		// A = N1+2*M1, B = N2+2*M2
 			// C = (K1+1)*(K2+1)
 int N1,N2;		//computed;
 double* S;		// S(-M1) ... S(M1+N1), S[A]
 double* T;		// T(-M2) ... T(M2+N2), T[B]
 double** W;		//Wagi C elementow, W[K1+1][K2+1]
 struct Point** P;	//Punkty kontr. C elementow P[K1+1][K2+1].{x,y,z}
 double U0, U1;		//Przedz parametryzacji w/g X
 double V0, V1;		//Przedz parametryzacji w/g Y
 			//wyliczone trojkaty:
 double* triangles;	//jeden trojkat to 9 kolejnych double
 int   ntriangles;	//ilosc trojkatow
};

struct IGS6Object
{
 int Did;		// K1 indeks gorny sumy
 int K;			// N = K-M+1
 int M;			//stopien krzywej bazowej
 int A;			// A = N+2*M
 int N;			//computed;
 double* T;		// T(-M) ... T(M+N), T[A+1]
 double* W;		//Wagi W[K+1]
 struct Point* P;	//Punkty kontr  P[K+1].{x,y,z}
 double U0, U1;		//Przedz parametryzacji w/g X
 			//wyliczone trojkaty:
 double* points;		//jeden trojkat to 9 kolejnych double
 int   npoints;		//ilosc trojkatow
};

static int WX = 0x200;
static int WY = 0x200;
static struct IGSObject* igs_array = NULL;
static int igs_num = 0;
static struct IGS6Object* igs6_array = NULL;
static int igs6_num = 0;
static int* ntr = 0;
static double** tr = 0;
static double** tr_buff = 0;
static int* npt = 0;
static double** pt = 0;
static double** pt_buff = 0;
void debug(char*,...);
void vdebug(char*,...);
#define BITS16
//#define BITS24
#ifdef BITS24

#define RGB RGB24
#define ReturnRed   RR24
#define ReturnGreen RG24
#define ReturnBlue  RB24

#endif
#ifdef BITS16

#define RGB RGB16
#define ReturnRed   RR16
#define ReturnGreen RG16
#define ReturnBlue  RB16

#endif


void flush_zbuff()
{
 int i,j;
 debug("flush_zbuff");
 if (want_triangles) return;
 for (i=0;i<WX;i++)
	 for (j=0;j<WY;j++)
	 {
	   if (zbuff[i][j].col != ZMAGIC)
	   {
	    XSetForeground(dsp,gc,zbuff[i][j].col);
	    XDrawPoint(dsp,win,gc,i,j);
	   }
	  zbuff[i][j].col = ZMAGIC;
	  zbuff[i][j].Z   = FAR_AWAY;
	 }
}


void clear_zbuff()
{
 int i,j;
 debug("clear_zbuff");
 if (want_triangles) return;
 for (i=0;i<WX;i++)
	 for (j=0;j<WY;j++) { zbuff[i][j].col = ZMAGIC; zbuff[i][j].Z = FAR_AWAY; }
}


void create_zbuff()
{
 int i;
 debug("create_zbuff");
 if (want_triangles) return;
 zbuff = (struct Zbuff**)malloc(sizeof(struct ZBuff*)*WX);
 for (i=0;i<WX;i++) zbuff[i] = (struct Zbuff*)malloc(sizeof(struct Zbuff)*WY);
}


void delete_zbuff()
{
 int i;
 debug("delete_zbuff");
 if (!zbuff) return;
 if (want_triangles) return;
 for (i=0;i<WX;i++) free(zbuff[i]);
 free(zbuff);
 zbuff = NULL;
}

void error(char*,...);

void info_out()
{
 debug("info_out");
 printf("START_AX=%d, START_AY=%d, START_AZ=%d\n",START_AX,START_AY,START_AZ);
 printf("STEP_AX=%d, STEP_AY=%d, STEP_AZ=%d\n",AX_STEP,AY_STEP,AZ_STEP);
 printf("TX_STEP=%d, TY_STEP=%d, TZ_STEP=%d\n",TX_STEP,TY_STEP,TZ_STEP);
 printf("XOFFSET=%f, YOFFSET=%f, ZOFFSET=%f\n",XOFFSET,YOFFSET,ZOFFSET);
 printf("MAXLINES=%d, TIMERINT=%d, LINESINC=%d\n",MAXLINES,TIMERINT,LINESINC);
 printf("SCALEX=%f, SCALEY=%f, LINECOL: %02x:%02x:%02x\n",SCALEX,SCALEY,line_r,line_g,line_b);
 printf("MILLCOL %02x:%02x:%02x,CURVCOL %02x:%02x:%02x,SURFCOL %02x:%02x:%02x\n",mill_r,mill_g,mill_b,igs_r,igs_g,igs_b,igs6_r,igs6_g,igs6_b);
 printf("zoff=%f, ldiff=%d, raw_triangs=%d, interval=%d\n",zoff, ldiff, want_triangles, interval);
 printf("LX=%f, LY=%f, LZ=%f, OINT=%f\n", LX,LY,LZ,OINT);
 printf("num_igs=%d, num_igs6=%d\n",igs_num,igs6_num);
 printf("WANTS:, CAM->%d, IGS->%d\n", want_cam, want_igs);
 printf("WX-WY: %dx%d\n",WX,WY);
 printf("\n");
 if (START_AX<=-360 || START_AX>=360) error("bad value: START_AX=%d", START_AX);
 if (START_AY<=-360 || START_AY>=360) error("bad value: START_AY=%d", START_AY);
 if (START_AZ<=-360 || START_AZ>=360) error("bad value: START_AZ=%d", START_AZ);
 if (AX_STEP<0  || AX_STEP>=360)  error("bad value: AX_STEP=%d", AX_STEP);
 if (AY_STEP<0  || AY_STEP>=360)  error("bad value: AY_STEP=%d", AY_STEP);
 if (AZ_STEP<0  || AZ_STEP>=360)  error("bad value: AZ_STEP=%d", AZ_STEP);
 if (TX_STEP<0  || TX_STEP>=360)  error("bad value: TX_STEP=%d", TX_STEP);
 if (TY_STEP<0  || TY_STEP>=360)  error("bad value: TY_STEP=%d", TY_STEP);
 if (TZ_STEP<0  || TZ_STEP>=360)  error("bad value: TZ_STEP=%d", TZ_STEP);
 if (XOFFSET>1e5) error("bad value: XOFFSET=%f", XOFFSET);
 if (YOFFSET>1e5) error("bad value: YOFFSET=%f", YOFFSET);
 if (ZOFFSET>1e5) error("bad value: ZOFFSET=%f", ZOFFSET);
 if (MAXLINES<1 || MAXLINES>2048) error("bad value: MAXLINES=%d", MAXLINES);
 if (TIMERINT<30000 || TIMERINT>60000000) error("bad value: TIMER_INT=%d", TIMERINT);
 if (LINESINC<1     || LINESINC>1024) error("bad value: LINESINC=%d", LINESINC);
 if (SCALEX==0 && SCALEY==0) error("bad value: SCALEX && SCALEY cannot be both 0");
 if (SCALEX>1e4) error("bad value: SCALEX too high: %f", SCALEX);
 if (SCALEY>1e4) error("bad value: SCALEY too high: %f", SCALEY);
 if (line_r<0 || line_r>=0x100) error("bad value: line_r=%d", line_r);
 if (line_g<0 || line_g>=0x100) error("bad value: line_g=%d", line_g);
 if (line_b<0 || line_b>=0x100) error("bad value: line_b=%d", line_b);
 if (zoff<0.0) printf("\n\nWARNING:\tzoff<0\n\n");
 if (ldiff<=1) printf("\n\nWARNING:\tdiff awfully small: %d\n\n", ldiff);
 if (fabs(LX)>1e5) printf("\n\nWARNING:\tlight far away from home: LX=%f\n\n", LX);
 if (fabs(LY)>1e5) printf("\n\nWARNING:\tlight far away from home: LY=%f\n\n", LY);
 if (fabs(LZ)>1e5) printf("\n\nWARNING:\tlight far away from home: LZ=%f\n\n", LZ);
 if (OINT<=0.0) error("bad value: OINT=%f\n", OINT);
 if (interval<=0 || interval>512) error("bad value: interval: %d", interval);
 if (WX<32 || WX>1024) error("bad value: WX: %d", WX);
 if (WY<32 || WY>768) error("bad value: WY: %d", WY);
}


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



void create_f_table()
{
 int i;
 debug("creating function table");
 for (i=0;i<360;i++)
   { sines[i] = sin(((double)i*3.1415926F)/180.0); cosines[i] = cos(((double)i*3.1415926F)/180.0); }
 angle_x = START_AX;
 angle_y = START_AY;
 angle_z = START_AZ;
}


void copy_to_buffers()
{
 int i,j;
 debug("copying points to buffers");
 if (want_cam)
   {
    for (i=0;i<npoints;i++)  path_buff[i] = path[i];
    for (i=0;i<MILL_POINTS;i++)      mill_buff[i] = mill[i];
   }
 if (want_igs)
   {
    for (i=0;i<igs_num;i++)
    for (j=0;j<ntr[i];j++) tr_buff[i][j] = tr[i][j];
    for (i=0;i<igs6_num;i++)
    for (j=0;j<npt[i];j++) pt_buff[i][j] = pt[i][j];
   }
 lxb=lx;
 lyb=ly;
 lzb=lz;
}


void copy_from_buffers()
{
 int i,j;
 debug("copying points from buffers");
 if (want_cam)
   {
    for (i=0;i<npoints;i++)  path[i] = path_buff[i];
    for (i=0;i<MILL_POINTS;i++)      mill[i] = mill_buff[i];
   }
 if (want_igs)
   {
    for (i=0;i<igs_num;i++)
    for (j=0;j<ntr[i];j++) tr[i][j] = tr_buff[i][j];
    for (i=0;i<igs6_num;i++)
    for (j=0;j<npt[i];j++) pt[i][j] = pt_buff[i][j];
   }
 lx=lxb;
 ly=lyb;
 lz=lzb;
}


void check_angles()
{
 debug("checking constraints");
 if (angle_x<0)     angle_x +=  360;
 if (angle_x>=360)  angle_x -=  360;
 if (angle_y<0)     angle_y +=  360;
 if (angle_y>=360)  angle_y -=  360;
 if (angle_z<0)     angle_z +=  360;
 if (angle_z>=360)  angle_z -=  360;
 if (lx> MAX_LPOS)      lx=MAX_LPOS;
 if (ly> MAX_LPOS)      ly=MAX_LPOS;
 if (lz> MAX_LPOS)      ly=MAX_LPOS;
 if (lx<-MAX_LPOS)      lx=-MAX_LPOS;
 if (ly<-MAX_LPOS)      ly=-MAX_LPOS;
 if (lz<-MAX_LPOS)      lz=-MAX_LPOS;
}


void world_transforms()
{
 int i,j;
 double y,x;
 debug("world transformations");
 check_angles();
 if (want_cam)
   {
    for (i=0;i<npoints;i++)
      {
       path[i].x += ox;
       path[i].y += oy;
       path[i].z += oz;
      }
    for (i=0;i<74;i++)
      {
       mill[3*i]   += ox;
       mill[3*i+1] += oy;
       mill[3*i+2] += oz;
      }
    for (i=0;i<npoints;i++)
      {
       y         = path[i].y*cosines[angle_x] - path[i].z*  sines[angle_x];
       path[i].z = path[i].y*  sines[angle_x] + path[i].z*cosines[angle_x];
       path[i].y = y;
      }
    for (i=0;i<npoints;i++)
      {
       x         = path[i].x*cosines[angle_y] - path[i].z*  sines[angle_y];
       path[i].z = path[i].x*  sines[angle_y] + path[i].z*cosines[angle_y];
       path[i].x = x;
      }
    for (i=0;i<npoints;i++)
      {
       x         = path[i].x*cosines[angle_z] - path[i].y*  sines[angle_z];
       path[i].y = path[i].x*  sines[angle_z] + path[i].y*cosines[angle_z];
       path[i].x = x;
      }
    for (i=0;i<74;i++)
      {
       y             = mill[3*i+1]*cosines[angle_x] - mill[3*i+2]*sines[angle_x];
       mill[3*i+2] = mill[3*i+1]*sines[angle_x] + mill[3*i+2]*cosines[angle_x];
       mill[3*i+1] = y;
      }
    for (i=0;i<74;i++)
      {
       x             = mill[3*i]*cosines[angle_y] - mill[3*i+2]*sines[angle_y];
       mill[3*i+2] = mill[3*i]*sines[angle_y]   + mill[3*i+2]*cosines[angle_y];
       mill[3*i]   = x;
      }
    for (i=0;i<74;i++)
      {
       x             = mill[3*i]*cosines[angle_z] - mill[3*i+1]*sines[angle_z];
       mill[3*i+1] = mill[3*i]*sines[angle_z]   + mill[3*i+1]*cosines[angle_z];
       mill[3*i] = x;
      }
   }
 if (want_igs)
   {
    for (j=0;j<igs_num;j++)
    for (i=0;i<ntr[j]/3;i++)
      {
       tr[j][3*i]   += ox;
       tr[j][3*i+1] += oy;
       tr[j][3*i+2] += oz;
      }
    for (j=0;j<igs_num;j++)
    for (i=0;i<ntr[j]/3;i++)
      {
       y            = tr[j][3*i+1]*cosines[angle_x] - tr[j][3*i+2]*sines[angle_x];
       tr[j][3*i+2] = tr[j][3*i+1]*sines[angle_x] + tr[j][3*i+2]*cosines[angle_x];
       tr[j][3*i+1] = y;
      }
    for (j=0;j<igs_num;j++)
    for (i=0;i<ntr[j]/3;i++)
      {
       x            = tr[j][3*i]*cosines[angle_y] - tr[j][3*i+2]*sines[angle_y];
       tr[j][3*i+2] = tr[j][3*i]*sines[angle_y]   + tr[j][3*i+2]*cosines[angle_y];
       tr[j][3*i]   = x;
      }
    for (j=0;j<igs_num;j++)
    for (i=0;i<ntr[j]/3;i++)
      {
       x            = tr[j][3*i]*cosines[angle_z] - tr[j][3*i+1]*sines[angle_z];
       tr[j][3*i+1] = tr[j][3*i]*sines[angle_z]   + tr[j][3*i+1]*cosines[angle_z];
       tr[j][3*i] = x;
      }
    for (j=0;j<igs6_num;j++)
    for (i=0;i<npt[j]/3;i++)
      {
       pt[j][3*i]   += ox;
       pt[j][3*i+1] += oy;
       pt[j][3*i+2] += oz;
      }
    for (j=0;j<igs6_num;j++)
    for (i=0;i<npt[j]/3;i++)
      {
       y            = pt[j][3*i+1]*cosines[angle_x] - pt[j][3*i+2]*sines[angle_x];
       pt[j][3*i+2] = pt[j][3*i+1]*sines[angle_x] + pt[j][3*i+2]*cosines[angle_x];
       pt[j][3*i+1] = y;
      }
    for (j=0;j<igs6_num;j++)
    for (i=0;i<npt[j]/3;i++)
      {
       x            = pt[j][3*i]*cosines[angle_y] - pt[j][3*i+2]*sines[angle_y];
       pt[j][3*i+2] = pt[j][3*i]*sines[angle_y]   + pt[j][3*i+2]*cosines[angle_y];
       pt[j][3*i]   = x;
      }
    for (j=0;j<igs6_num;j++)
    for (i=0;i<npt[j]/3;i++)
      {
       x            = pt[j][3*i]*cosines[angle_z] - pt[j][3*i+1]*sines[angle_z];
       pt[j][3*i+1] = pt[j][3*i]*sines[angle_z]   + pt[j][3*i+1]*cosines[angle_z];
       pt[j][3*i] = x;
      }
   }
}


void ZDrawLine(int sx, int sy, int ex, int ey, int sz, int ez, ulong lcol)
{
 int i,y,z,prvy,y1,y2,j,x,x1,x2;
 if (want_triangles)
   {
    XSetForeground(dsp,gc,RGB(line_r, line_g, line_b));
    XDrawLine(dsp,win,gc,sx,sy,ex,ey);
    return;
   }
 if (sx<0 || sy<0 || ex<0 || ey<0) return;
 if (sx>=WX || sy>=WY || ex>=WX || ey>=WY) return;
 if (sx==ex)
   {
    y1 = (sy<ey)?sy:ey;
    y2 = (sy<ey)?ey:sy;
    for (j=y1;j<y2;j++)
      {
       z = sz + ((j-sy)*(ez-sz))/(ey-sy);
       if (zbuff[sx][j].Z>z)
         {
          zbuff[sx][j].Z = z;
	  zbuff[sx][j].col = lcol;
         }
      }
    return;
   }
 y = sy;
 x1 = (sx<ex)?sx:ex;
 x2 = (sx<ex)?ex:sx;
 for (i=x1;i<=x2;i++)
 {
  prvy = y;
  y = sy + ((i-sx)*(ey-sy))/(ex-sx);
  z = sz + ((i-sx)*(ez-sz))/(ex-sx);
  if (zbuff[i][y].Z>z)
   {
    zbuff[i][y].Z = z;
    zbuff[i][y].col = lcol;
   }
  if (fabs(y-prvy)>1)
   {
    y1 = (y>prvy)?prvy:y;
    y2 = (y<prvy)?prvy:y;
    for (j=y1;j<=y2;j++)
     {
      x = sx + ((j-sy)*(ex-sx))/(ey-sy);
      z = sz + ((j-sy)*(ez-sz))/(ey-sy);
      if (zbuff[x][j].Z>z)
        {
         zbuff[x][j].Z = z;
	 zbuff[x][j].col = lcol;
        }
     }
   }
 }
}


void check_lights()
{
 vdebug("check_lights\n");
 if (want_triangles) return;
 if (light_1<0.0) light_1=0.0;
 if (light_2<0.0) light_2=0.0;
 if (light_3<0.0) light_3=0.0;
 if (light_1>0.99) light_1=0.99;
 if (light_2>0.99) light_2=0.99;
 if (light_3>0.99) light_3=0.99;
}


void XNormalF(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3)
{
 double Ax,Ay,Az,Bx,By,Bz,l;
 double N1x,N1y,N1z,N2x,N2y,N2z,N3x,N3y,N3z;
 if (want_triangles) return;
 Ax = x2-x1;
 Bx = x3-x1;
 Ay = y2-y1;
 By = y3-y1;
 Az = z2-z1;
 Bz = z3-z1;
 nx = Ay*Bz-Az*By;
 ny = Az*Bx-Ax*Bz;
 nz = Ax*By-Ay*Bx;
 l = sqrt(nx*nx+ny*ny+nz*nz);
 nx /= l;
 ny /= l;
 nz /= l;
 N1x = lx - x1;
 N1y = ly - y1;
 N1z = lz - z1;
 N2x = lx - x2;
 N2y = ly - y2;
 N2z = lz - z2;
 N3x = lx - x3;
 N3y = ly - y3;
 N3z = lz - z3;
 l = sqrt(N1x*N1x+N1y*N1y+N1z*N1z);
 N1x /= l;
 N1y /= l;
 N1z /= l;
 l = sqrt(N2x*N2x+N2y*N2y+N2z*N2z);
 N2x /= l;
 N2y /= l;
 N2z /= l;
 l = sqrt(N3x*N3x+N3y*N3y+N3z*N3z);
 N3x /= l;
 N3y /= l;
 N3z /= l;
 light_1 = nx*N1x+ny*N1y+nz*N1z;
 light_2 = nx*N2x+ny*N2y+nz*N2z;
 light_3 = nx*N3x+ny*N3y+nz*N3z;
 check_lights();
}


void draw_light()
{
 debug("draw_light");
 if (want_triangles) return;
 XSetForeground(dsp,gc,RGB(0x40,0x40,0xFF));
 XFillArc(dsp,win,gc,t1x,t1y,10,10,0,ENTIRE_ARC);
}


void draw_triangle(ulong col, double z1, double z2, double z3)
{
 int wmp,wst,wen;
 double fmp,cz;
 int w1x,w2x,w1y,w2y,w3x,w3y,i;
 int minx,maxx,minxy,maxxy,j;
 double start,end,step_st,step_en,step_zst,step_zen;
 double ist,ien;
 double base_r,base_g,base_b;
 double delta_r,delta_g,delta_b;
 double startr,startg,startb,startz;
 double endr,endg,endb,endz;
 double minxl,maxxl,minxz,maxxz;
 double step_lrst,step_lren;
 double step_lgst,step_lgen;
 double step_lbst,step_lben;
 int   solid;
 int cx,cy;
 cx = WX;
 cy = WY;
 ulong solid_col;
 if (z1 < -1.5*zoff || z2 < -1.5*zoff || z3 < -1.5*zoff) return;
 if (want_triangles)
 {
  XSetForeground(dsp,gc,RGB(line_b, line_g, line_r));
  XDrawLine(dsp,win,gc,t1x,t1y,t2x,t2y);
  XDrawLine(dsp,win,gc,t1x,t1y,t3x,t3y);
  XDrawLine(dsp,win,gc,t2x,t2y,t3x,t3y);
  return ;
 }
 solid = 0;
 solid_col = 0;
 if (light_1<0.0)
   {
    solid = 1;
    solid_col = RGB(ReturnRed(col)/ldiff,ReturnGreen(col)/ldiff,ReturnBlue(col)/ldiff);
   }
 base_r = (double)ReturnRed(col)/(double)ldiff;
 base_g = (double)ReturnGreen(col)/(double)ldiff;
 base_b = (double)ReturnBlue(col)/(double)ldiff;		/*gdy brak swiatla*/
 delta_r = base_r*((double)ldiff-1.0);
 delta_g = base_g*((double)ldiff-1.0);
 delta_b = base_b*((double)ldiff-1.0);
 w1x=t1x;
 w1y=t1y;
 w2x=t2x;
 w2y=t2y;
 w3x=t3x;
 w3y=t3y;
 if (w2y>w3y)
   {
    wmp = w2y;
    w2y = w3y;
    w3y = wmp;
    wmp = w2x;
    w2x = w3x;
    w3x = wmp;
    fmp = light_2;
    light_2 = light_3;
    light_3 = fmp;
    fmp = z2;
    z2 = z3;
    z3 = fmp;
   }
 if (w1y>w2y)
   {
    wmp = w1y;
    w1y = w2y;
    w2y = wmp;
    wmp = w1x;
    w1x = w2x;
    w2x = wmp;
    fmp = light_1;
    light_1 = light_2;
    light_2 = fmp;
    fmp = z1;
    z1 = z2;
    z2 = fmp;
   }
 if (w2y>w3y)
   {
    wmp = w2y;
    w2y = w3y;
    w3y = wmp;
    wmp = w2x;
    w2x = w3x;
    w3x = wmp;
    fmp = light_2;
    light_2 = light_3;
    light_3 = fmp;
    fmp = z2;
    z2 = z3;
    z3 = fmp;
   }			/* szybki sort po Y-kach*/
 minx  = w2x;
 maxx  = w3x;
 minxy = w2y;
 maxxy = w3y;
 minxl = light_2;
 maxxl = light_3;
 minxz = z2;
 maxxz = z3;
 if (minx>maxx)		/* ustalenie startowych X-ow, swiatel, Ztow natezen barw itp. */
  {
    wmp   = minx;
    minx  = maxx;
    maxx  = wmp;
    wmp   = minxy;
    minxy = maxxy;
    maxxy = wmp;
    fmp   = minxl;
    minxl = maxxl;
    maxxl = fmp;
    fmp   = minxz;
    minxz = maxxz;
    maxxz = fmp;
  }
 start = (double)w1x;
 end   = (double)w1x;
 step_st = ((double)minx-(double)w1x)/((double)minxy-(double)w1y);
 step_en = ((double)maxx-(double)w1x)/((double)maxxy-(double)w1y);
 startz = z1;
 endz   = z1;
 step_zst = (minxz-z1)/((double)minxy-(double)w1y);
 step_zen = (maxxz-z1)/((double)maxxy-(double)w1y);
 startr = base_r + light_1*delta_r;
 startg = base_g + light_1*delta_g;
 startb = base_b + light_1*delta_b;
 endr=startr;
 endg=startg;
 endb=startb;
 step_lrst = ((minxl-light_1)/((double)minxy-(double)w1y))*delta_r;
 step_lren = ((maxxl-light_1)/((double)maxxy-(double)w1y))*delta_r;
 step_lgst = ((minxl-light_1)/((double)minxy-(double)w1y))*delta_g;
 step_lgen = ((maxxl-light_1)/((double)maxxy-(double)w1y))*delta_g;
 step_lbst = ((minxl-light_1)/((double)minxy-(double)w1y))*delta_b;
 step_lben = ((maxxl-light_1)/((double)maxxy-(double)w1y))*delta_b;
 ist = w1y;
 ien = w2y;
 if (ist<0) ist=0;
 if (ien>=cy) ien = cy-1;
 for (i=ist;i<ien;i++)		/* pierwsza czesc trojkata */
 {
  start+=step_st;
  end+=step_en;
  startz+=step_zst;
  endz+=step_zen;
  startr+=step_lrst;
  endr+=step_lren;
  startg+=step_lgst;
  endg+=step_lgen;
  startb+=step_lbst;
  endb+=step_lben;
  wst=(int)start+1;
  wen=(int)end;
  if (wen<wst) { wmp=wst; wst=wen; wen=wmp; }
  if (wst<0) wst=0;
  if (wen>=cx) wen = cx-1;
  for (j=wst;j<=wen;j++)
    {
     cz = startz+(((double)j-start)/(end-start)*(endz-startz));
     if (zbuff[j][i].Z>cz)
       {
	zbuff[j][i].Z=cz;
	if (!solid) zbuff[j][i].col = RGB((int)(startr+(((double)j-start)/(end-start))*(endr-startr)),
			     (int)(startg+(((double)j-start)/(end-start))*(endg-startg)),
			     (int)(startb+(((double)j-start)/(end-start))*(endb-startb)));
	else zbuff[j][i].col = solid_col;
       }
    }
 }
 if (minx==w2x)
  {
   start = (double)minx;
   step_st = ((double)maxx-(double)minx)/((double)maxxy-(double)minxy);
   startr = base_r + minxl*delta_r;
   startg = base_g + minxl*delta_g;
   startb = base_b + minxl*delta_b;
   startz = minxz;
   step_zst = (z3-minxz)/((double)maxxy-(double)minxy);
   step_lrst = (light_3-minxl)/((double)maxxy-(double)minxy)*delta_r;
   step_lgst = (light_3-minxl)/((double)maxxy-(double)minxy)*delta_g;
   step_lbst = (light_3-minxl)/((double)maxxy-(double)minxy)*delta_b;
  }
 else
  {
   end = (double)maxx;
   step_en = ((double)minx-(double)maxx)/((double)minxy-(double)maxxy);
   endr = base_r + maxxl*delta_r;
   endg = base_g + maxxl*delta_g;
   endb = base_b + maxxl*delta_b;
   endz = maxxz;
   step_zen = (maxxz-z3)/((double)maxxy-(double)minxy);
   step_lren = (maxxl-light_3)/((double)maxxy-(double)minxy)*delta_r;
   step_lgen = (maxxl-light_3)/((double)maxxy-(double)minxy)*delta_g;
   step_lben = (maxxl-light_3)/((double)maxxy-(double)minxy)*delta_b;
  }
 if (end<start)						/* gdy poczatek dalej niz koniec */
  {
   fmp       = end;
   end       = start;
   start     = fmp;
   fmp       = step_st;
   step_st   = step_en;
   step_en   = fmp;
   fmp       = endz;   		/* not needed?*/
   endz      = startz;
   startz    = fmp;
   fmp       = step_zst;
   step_zst  = step_zen;
   step_zen  = fmp;      	/* to here? */
   fmp       = startr;
   startr    = endr;
   endr      = fmp;
   fmp       = startg;
   startg    = endg;
   endg      = fmp;
   fmp       = startb;
   startb    = endb;
   endb      = fmp;
   fmp       = step_lrst;
   step_lrst = step_lren;
   step_lren = fmp;
   fmp       = step_lgst;
   step_lgst = step_lgen;
   step_lgen = fmp;
   fmp       = step_lbst;
   step_lbst = step_lben;
   step_lben = fmp;
  }
 ist = w2y;
 ien = w3y;
 if (ist<0) ist=0;
 if (ien>=cy) ien = cy-1;
 for (i=ist;i<ien;i++)		/* druga czesc trojkata */
 {
  start+=step_st;
  end+=step_en;
  startz+=step_zst;
  endz+=step_zen;
  startr+=step_lrst;
  endr+=step_lren;
  startg+=step_lgst;
  endg+=step_lgen;
  startb+=step_lbst;
  endb+=step_lben;
  wst = (int)start+1;
  wen = (int)end;
  if (wst<0) wst=0;
  if (wen>=cx) wen = cx-1;
  for (j=wst;j<=wen;j++)
    {
     cz = startz+(((double)j-start)/(end-start)*(endz-startz));
     if (zbuff[j][i].Z>cz)
       {
	zbuff[j][i].Z=cz;
	if (!solid) zbuff[j][i].col = RGB((int)(startr+(((double)j-start)/(end-start))*(endr-startr)),
			     (int)(startg+(((double)j-start)/(end-start))*(endg-startg)),
			     (int)(startb+(((double)j-start)/(end-start))*(endb-startb)));
	else zbuff[j][i].col = solid_col;
       }
    }
 }
}


void draw_bspline_2D()
{
 debug("draw_bspline_2D");
 double perspectivex;
 double perspectivey;
 int wx,wy;
 wx=WX/2;
 wy=WY/2;
 perspectivex = (double)WX/2.0;
 perspectivey = (double)WY/2.0;
 int T1x,T2x,T1y,T2y;
 ulong curr_col;
 int i,j;
 if (!want_igs || !igs6_num) return;
 curr_col = RGB(igs6_r,igs6_g,igs6_b);
 XSetForeground(dsp,gc,curr_col);
 for (j=0;j<igs6_num;j++)
 for (i=3;i<npt[j]/3;i++)
      {
       T1x = wx+(int)(SCALEX*(perspectivex*pt[j][3*i-3]/(zoff+pt[j][3*i-1])));
       T1y = wy-(int)(SCALEY*(perspectivey*pt[j][3*i-2]/(zoff+pt[j][3*i-1])));
       T2x = wx+(int)(SCALEX*(perspectivex*pt[j][3*i  ]/(zoff+pt[j][3*i+2])));
       T2y = wy-(int)(SCALEY*(perspectivey*pt[j][3*i+1]/(zoff+pt[j][3*i+2])));
       if (!want_triangles) ZDrawLine(T1x,T1y,T2x,T2y,pt[j][3*i-1],pt[j][3*i+2],curr_col);
       else XDrawLine(dsp,win,gc,T1x,T1y,T2x,T2y);
      }
}


void draw_bspline()
{
 debug("draw_bspline");
 ulong curr_col;
 int i,j;
 if (!want_igs || !igs_num) return;
/* printf("DRAW BSPLINE!\n");*/
 double perspectivex;
 double perspectivey;
 int wx,wy;
 wx=WX/2;
 wy=WY/2;
 perspectivex = (double)WX/2.0;
 perspectivey = (double)WY/2.0;
 curr_col = RGB(igs_r,igs_g,igs_b);
 for (i=0;i<igs_num;i++)
 for (j=0;j<ntr[i];j+=9)
   {
/*    printf("%f\n", tr[i][3*j]);*/
    t1x = wx + (int)(SCALEX*(perspectivex*tr[i][j  ]/(zoff+tr[i][j+2])));
    t1y = wy - (int)(SCALEY*(perspectivey*tr[i][j+1]/(zoff+tr[i][j+2])));
    t2x = wx + (int)(SCALEX*(perspectivex*tr[i][j+3]/(zoff+tr[i][j+5])));
    t2y = wy - (int)(SCALEY*(perspectivey*tr[i][j+4]/(zoff+tr[i][j+5])));
    t3x = wx + (int)(SCALEX*(perspectivex*tr[i][j+6]/(zoff+tr[i][j+8])));
    t3y = wy - (int)(SCALEY*(perspectivey*tr[i][j+7]/(zoff+tr[i][j+8])));
/*    printf("%d-%d-%d-%d-%d-%d\n",t1x,t1y,t2x,t2y,t3x,t3y);*/
    XNormalF(tr[i][j  ], tr[i][j+1], tr[i][j+2],
             tr[i][j+3], tr[i][j+4], tr[i][j+5],
             tr[i][j+6], tr[i][j+7], tr[i][j+8]);
    draw_triangle(curr_col, tr[i][j+2],tr[i][j+5],tr[i][j+8]);
   }
}


void draw_tool()
{
 debug("draw_tool");
 ulong curr_col;
 int i,ni,j,nj;
 if (!want_cam) return;
 double perspectivex;
 double perspectivey;
 int wx,wy;
 wx=WX/2;
 wy=WY/2;
 perspectivex = (double)WX/2.0;
 perspectivey = (double)WY/2.0;
 debug("draw_tool");
 curr_col = RGB(mill_r,mill_g,mill_b);
 t1x = wx + (int)(SCALEX*(perspectivex*mill[0 ]/(zoff+mill[2 ])));
 t1y = wy - (int)(SCALEY*(perspectivey*mill[1 ]/(zoff+mill[2 ])));
 for (i=2;i<37;i++)
   {
    t2x = wx + (int)(SCALEX*(perspectivex*mill[3*i  ]/(zoff+mill[3*i+2])));
    t2y = wy - (int)(SCALEY*(perspectivey*mill[3*i+1]/(zoff+mill[3*i+2])));
    t3x = wx + (int)(SCALEX*(perspectivex*mill[3*i+3]/(zoff+mill[3*i+5])));
    t3y = wy - (int)(SCALEY*(perspectivey*mill[3*i+4]/(zoff+mill[3*i+5])));
    XNormalF(mill[0],mill[1],mill[2],
            mill[3*i],mill[3*i+1],mill[3*i+2],
	    mill[3*i+3],mill[3*i+4],mill[3*i+5]);
    draw_triangle(curr_col, mill[2],mill[3*i+2],mill[3*i+5]);
   }
 t2x = wx + (int)(SCALEX*(perspectivex*mill[111]/(zoff+mill[113])));
 t2y = wy - (int)(SCALEY*(perspectivey*mill[112]/(zoff+mill[113])));
 t3x = wx + (int)(SCALEX*(perspectivex*mill[6  ]/(zoff+mill[8  ])));
 t3y = wy - (int)(SCALEY*(perspectivey*mill[7  ]/(zoff+mill[8  ])));
 XNormalF(mill[0],mill[1],mill[2],
            mill[111],mill[112],mill[113],
	    mill[6  ],mill[7  ],mill[8  ]);
 draw_triangle(curr_col,mill[2],mill[113],mill[8]);
 t1x = wx + (int)(SCALEX*(perspectivex*mill[3 ]/(zoff+mill[5 ])));
 t1y = wy - (int)(SCALEY*(perspectivey*mill[4 ]/(zoff+mill[5 ])));
 for (i=38;i<73;i++)
   {
    t3x = wx + (int)(SCALEX*(perspectivex*mill[3*i  ]/(zoff+mill[3*i+2])));
    t3y = wy - (int)(SCALEY*(perspectivey*mill[3*i+1]/(zoff+mill[3*i+2])));
    t2x = wx + (int)(SCALEX*(perspectivex*mill[3*i+3]/(zoff+mill[3*i+5])));
    t2y = wy - (int)(SCALEY*(perspectivey*mill[3*i+4]/(zoff+mill[3*i+5])));
    XNormalF(mill[3],mill[4],mill[5],
	    mill[3*i+3],mill[3*i+4],mill[3*i+5],
            mill[3*i],mill[3*i+1],mill[3*i+2]);
    draw_triangle(curr_col, mill[5],mill[3*i+5],mill[3*i+2]);
   }
 t3x = wx + (int)(SCALEX*(perspectivex*mill[219]/(zoff+mill[221])));
 t3y = wy - (int)(SCALEY*(perspectivey*mill[220]/(zoff+mill[221])));
 t2x = wx + (int)(SCALEX*(perspectivex*mill[117]/(zoff+mill[119])));
 t2y = wy - (int)(SCALEY*(perspectivey*mill[118]/(zoff+mill[119])));
 XNormalF(mill[3],mill[4],mill[5],
	    mill[117],mill[118],mill[119],
            mill[219],mill[220],mill[221]);
 draw_triangle(curr_col,mill[5],mill[119],mill[221]);
 for (i=2;i<37;i++)
   {
    ni = i+1;
    j  = i+36;
    nj = j+1;
    if (ni==37) ni = 2;
    if (nj==73) nj = 38;
    t1x = wx + (int)(SCALEX*(perspectivex*mill[3*i   ]/(zoff+mill[3*i+ 2])));
    t1y = wy - (int)(SCALEY*(perspectivey*mill[3*i+1 ]/(zoff+mill[3*i+ 2])));
    t3x = wx + (int)(SCALEX*(perspectivex*mill[3*ni  ]/(zoff+mill[3*ni+2])));
    t3y = wy - (int)(SCALEY*(perspectivey*mill[3*ni+1]/(zoff+mill[3*ni+2])));
    t2x = wx + (int)(SCALEX*(perspectivex*mill[3*j   ]/(zoff+mill[3*j+ 2])));
    t2y = wy - (int)(SCALEY*(perspectivey*mill[3*j+1 ]/(zoff+mill[3*j+ 2])));
    XNormalF(
	    mill[3*i],mill[3*i+ 1],mill[3*i +2],
            mill[3*j], mill[3*j+ 1],mill[3*j +2],
	    mill[3*ni],mill[3*ni+1],mill[3*ni+2]);
    draw_triangle(curr_col, mill[3*i+2],mill[3*j+2],mill[3*ni+2]);
    t1x = wx + (int)(SCALEX*(perspectivex*mill[3*j   ]/(zoff+mill[3*j+ 2])));
    t1y = wy - (int)(SCALEY*(perspectivey*mill[3*j+1 ]/(zoff+mill[3*j+ 2])));
    t3x = wx + (int)(SCALEX*(perspectivex*mill[3*ni  ]/(zoff+mill[3*ni+2])));
    t3y = wy - (int)(SCALEY*(perspectivey*mill[3*ni+1]/(zoff+mill[3*ni+2])));
    t2x = wx + (int)(SCALEX*(perspectivex*mill[3*nj  ]/(zoff+mill[3*nj+2])));
    t2y = wy - (int)(SCALEY*(perspectivey*mill[3*nj+1]/(zoff+mill[3*nj+2])));
    XNormalF(
	    mill[3*j],mill[3*j+ 1],mill[3*j +2],
            mill[3*nj],mill[3*nj+1],mill[3*nj+2],
	    mill[3*ni],mill[3*ni+1],mill[3*ni+2]);
    draw_triangle(curr_col, mill[3*j+2],mill[3*nj+2],mill[3*ni+2]);
   }
}


void move_tool(int idx)
{
 debug("move_tool");
 int i;
 int to = MILL_POINTS/3;
 if (!want_cam) return;
 debug("move_tool");
 if (idx<0 || idx>=npoints) return;
 for (i=0;i<to;i++)
  {
   mill[3*i]   += path[idx].x;
   mill[3*i+1] += path[idx].y;
   mill[3*i+2] += path[idx].z;
  }
}


void refresh_window(int isBlack)
{
 int i,T1x,T2x,T1y,T2y,from,to;
 int mid;
 double perspectivex;
 double perspectivey;
 int wx,wy;
 wx=WX/2;
 wy=WY/2;
 perspectivex = (double)WX/2.0;
 perspectivey = (double)WY/2.0;
 from = to = 0;
 ulong solid_col;
 debug("refreshing window, arg=%d",isBlack);
 if (isBlack) { clear_zbuff(); XClearWindow(dsp,win); return; }
 copy_from_buffers();
 if (want_cam)
   {
    from = ((time_async-MAXLINES)>=1)?(time_async-MAXLINES):1;
    to = (time_async<npoints)?time_async:npoints;
    mid = (from+to)/2;
    move_tool(to);
   }
 world_transforms();
 if (want_cam) debug("drawing lines, from:%d, to:%d",from,to);
 solid_col = RGB(line_r, line_g, line_b);
 if (want_igs) { draw_bspline(); draw_bspline_2D(); }
 if (want_cam)
   {
    XSetForeground(dsp,gc,solid_col);
    for (i=from;i<to;i++)
      {
       T1x = wx+(int)(SCALEX*(perspectivex*path[i-1].x/(zoff+path[i-1].z)));
       T1y = wy-(int)(SCALEY*(perspectivey*path[i-1].y/(zoff+path[i-1].z)));
       T2x = wx+(int)(SCALEX*(perspectivex*path[i].x/(zoff+path[i].z)));
       T2y = wy-(int)(SCALEY*(perspectivey*path[i].y/(zoff+path[i].z)));
       if (!want_triangles)  ZDrawLine(T1x,T1y,T2x,T2y,path[i-1].z,path[i].z,solid_col);
       else XDrawLine(dsp,win,gc,T1x,T1y,T2x,T2y);
      }
    draw_tool();
   }
 debug("window refreshed");
 if (!want_triangles) flush_zbuff();
 t1x = wx + (int)(SCALEX*(perspectivex*lx/(zoff+lz)));
 t1y = wy - (int)(SCALEY*(perspectivey*ly/(zoff+lz)));
 draw_light();
}


int readline(FILE* f, char* str)
{
 vdebug("readline");
 int zn;
 int i;
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
 vdebug("param_id"); 
 char scan[82];
 int idx;
 strcpy(scan, line);
 scan[72] = 0;
 idx = 71;
 while (scan[idx] != ' ' && idx>=0) idx--;
 return atoi(scan+idx);
}


int get_word_count(char* line)
{
 vdebug("get_word_count"); 
 int i;
 int words = 0;
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
 debug("fill_args"); 
 int i;
 int idx;
 int wc;
 idx = 0;
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


double bspline_CALKA_Z(double s, double t, int K1, int K2, int idx)
{
 int i,j;
 double up,down;
 double factor;
 int M1,M2;
 double b1,b2;
/* return 0.0;*/
 vdebug("bspline_CALKA_Z\n");
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


double bspline_CALKA_Y(double s, double t, int K1, int K2, int idx)
{
 int i,j;
 double up,down;
 double factor;
 int M1,M2;
 double b1,b2;
 vdebug("bspline_CALKA_Y\n");
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


double bspline_CALKA_X(double s, double t, int K1, int K2, int idx)
{
 int i,j;
 double up,down;
 double factor;
 int M1,M2;
 double b1,b2;
 vdebug("bspline_CALKA_X\n");
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


double bspline_CALKA_126Y(double t, int K, int idx)
{
 int i;
 double up,down;
 double factor;
 int N,M;
 double b;
 vdebug("bspline_CALKA_126Y\n");
/* return 0.0;*/
 N = igs6_array[idx].N;
 M = igs6_array[idx].M;
 up = 0.0;
 down = 0.0;
 for (i=0;i<=K;i++)
   {
    b = basis_bspline(igs6_array[idx].T, t, M, M, i-M, i+1);
    vdebug("6Y> b = %f\n", b);
    factor = igs6_array[idx].W[i]*igs6_array[idx].P[i].y*b;
    up += factor;
    factor = igs6_array[idx].W[i]*b;
    down += factor;
   }
 if (down==0.0) return 0.0;
 vdebug("6Y> %f\n", up/down);
 return up/down;
}


double bspline_CALKA_126X(double t, int K, int idx)
{
 int i;
 double up,down;
 double factor;
 int N,M;
 double b;
 vdebug("bspline_CALKA_126X\n");
/* return 0.0;*/
 N = igs6_array[idx].N;
 M = igs6_array[idx].M;
 up = 0.0;
 down = 0.0;
 for (i=0;i<=K;i++)
   {
    b = basis_bspline(igs6_array[idx].T, t, M, M, i-M, i+1);
    vdebug("6X> b = %f\n", b);
    factor = igs6_array[idx].W[i]*igs6_array[idx].P[i].x*b;
    up += factor;
    factor = igs6_array[idx].W[i]*b;
    down += factor;
   }
 if (down==0.0) return 0.0;
 vdebug("6X> %f\n", up/down);
 return up/down;
}


void generate_structures_126(double* args, int nargs, int idx)
{
 int M,K,N,A,i,x;
 double U0,U1,xm,step;
 K = (int)args[1];
 M = (int)args[2];
 N = K-M+1;
 A = N+2*M;
 debug("generate_structures_126> K=%d, M=%d, N=%d, A=%d\n",K,M,N,A);
 igs6_array[idx].M  = M;
 igs6_array[idx].K  = K;
 igs6_array[idx].N  = N;
 igs6_array[idx].A  = A;
 if (args[3]!=0 || args[4]!=0 || args[5]!=1 || args[6]!=0) printf("PROP invalid for load_126\n");
 igs6_array[idx].T = (double*)malloc((A+1)*sizeof(double));
 if (!igs6_array[idx].T) error("malloc T");
 for (i=0;i<=A;i++) igs6_array[idx].T[i] = args[7+i];
 igs6_array[idx].W = (double*)malloc((K+1)*sizeof(double));
 if (!igs6_array[idx].W) error("malloc W");
 for (i=0;i<=K;i++) igs6_array[idx].W[i] = args[8+A+i];
 igs6_array[idx].P = (struct Point*)malloc((K+1)*sizeof(struct Point));
 if (!igs6_array[idx].P) error("malloc P");
 for (i=0;i<=K;i++)
      {
       igs6_array[idx].P[i].x = args[9+A+K+3*i];
       igs6_array[idx].P[i].y = args[9+A+K+3*i+1];
       igs6_array[idx].P[i].z = args[9+A+K+3*i+2];
      }
 U0 = args[12+A+4*K];
 U1 = args[13+A+4*K];
 igs6_array[idx].U0 = U0;
 igs6_array[idx].U1 = U1;
 debug("T\n");
 for (i=0;i<=A;i++) debug("%f\n", igs6_array[idx].T[i]);
 debug("W\n");
 for (i=0;i<=K;i++)
    debug("%2.3f\n", igs6_array[idx].W[i]);
 debug("P\n");
 for (i=0;i<=K;i++)
      debug("(%2.3f,%3.2f,%3.2f)\n",
    igs6_array[idx].P[i].x,igs6_array[idx].P[i].y,igs6_array[idx].P[i].z);
 debug("U=[%f,%f]\n", U0,U1);
 igs6_array[idx].npoints = 3*(interval+1);
 igs6_array[idx].points = (double*)malloc(igs6_array[idx].npoints*sizeof(double));
 if (!igs6_array[idx].points) error("malloc failed while generating points");
 x = 0;
 step = (U1-U0)/(double)interval;
 debug("step = %f\n", step);
 double t1,t2;
 for (xm=U0;xm<U1;xm+=step)
   {
     if (xm>U1-eps) continue;
     if (xm<U0+eps) continue;			//FIXME is this needed??
     t1 = bspline_CALKA_126X(xm, K, idx);
     t2 = bspline_CALKA_126Y(xm, K, idx);
     debug("XY> %f,%f",t1,t2);
     igs6_array[idx].points[x++] = t1;
     igs6_array[idx].points[x++] = t2;
     igs6_array[idx].points[x++] = 0.0;
   }
 debug("X = %d, npt = %d\n", x, igs6_array[idx].npoints);
 igs6_array[idx].npoints = x;
}


void generate_structures(double* args, int nargs, int idx)
{
 int M1,M2,K1,K2,N1,N2,A,B,C,i,j,x;
 double U0,U1,V0,V1,xm,ym,stepx,stepy;
 K1 = (int)args[2];
 K2 = (int)args[1];			//FIXME points are inverted in my structures :-(
 M1 = (int)args[3];			
 M2 = (int)args[4];
 N1 = (int)K1-M1+1;
 N2 = (int)K2-M2+1;
 A = N2+2*M1;				//FIXME, SO A<->B also, this causes problems only...
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
 /*for (i=0;i<K2;i++)
 for (j=0;j<K1;j++)
 for (x=0;x<K1-(j+1);x++)
   {
    if (igs_array[idx].P[i][j].y > igs_array[idx].P[i][j+1].y)
      {
       U0 = igs_array[idx].P[i][j].x;
       igs_array[idx].P[i][j].x = igs_array[idx].P[i][j+1].x;
       igs_array[idx].P[i][j+1].x = U0;
       U0 = igs_array[idx].P[i][j].y;
       igs_array[idx].P[i][j].y = igs_array[idx].P[i][j+1].y;
       igs_array[idx].P[i][j+1].y = U0;
       U0 = igs_array[idx].P[i][j].z;
       igs_array[idx].P[i][j].z = igs_array[idx].P[i][j+1].z;
       igs_array[idx].P[i][j+1].z = U0;
      }
   }*/
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
 double t1,t2,t3;
 for (xm=U0;xm<U1;xm+=stepx)
 for (ym=V0;ym<V1;ym+=stepy)
    {
     if (xm>U1-stepy-eps || ym>V1-stepy-eps) continue;
     if (xm<U0+eps || ym<V0+eps) continue;
     t1 = bspline_CALKA_X(xm, ym, K1, K2, idx);
     t2 = bspline_CALKA_Y(xm, ym, K1, K2, idx);
     t3 = bspline_CALKA_Z(xm, ym, K1, K2, idx);
     debug("XYZ> %f,%f,%f",t1,t2,t3);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: %f,%f\n", xm,ym);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_CALKA_X(xm+stepx, ym+stepy, K1, K2, idx);
     t2 = bspline_CALKA_Y(xm+stepx, ym+stepy, K1, K2, idx);
     t3 = bspline_CALKA_Z(xm+stepx, ym+stepy, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: +%f,+%f\n", xm+stepx,ym+stepy);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_CALKA_X(xm+stepx, ym, K1, K2, idx);
     t2 = bspline_CALKA_Y(xm+stepx, ym, K1, K2, idx);
     t3 = bspline_CALKA_Z(xm+stepx, ym, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: +%f,%f\n", xm+stepx,ym);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_CALKA_X(xm, ym, K1, K2, idx);
     t2 = bspline_CALKA_Y(xm, ym, K1, K2, idx);
     t3 = bspline_CALKA_Z(xm, ym, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: %f,%f\n", xm,ym);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_CALKA_X(xm, ym+stepy, K1, K2, idx);
     t2 = bspline_CALKA_Y(xm, ym+stepy, K1, K2, idx);
     t3 = bspline_CALKA_Z(xm, ym+stepy, K1, K2, idx);
     if (t1==0.0 && t2==0.0 && t3==0.0) printf("ZERO: %f,+%f\n", xm,ym+stepy);
     igs_array[idx].triangles[x++] = t1;
     igs_array[idx].triangles[x++] = t2;
     igs_array[idx].triangles[x++] = t3;
     t1 = bspline_CALKA_X(xm+stepx, ym+stepy, K1, K2, idx);
     t2 = bspline_CALKA_Y(xm+stepx, ym+stepy, K1, K2, idx);
     t3 = bspline_CALKA_Z(xm+stepx, ym+stepy, K1, K2, idx);
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


void load_126(FILE* f, int n)
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
 debug("load_126");
 lines = (int*)malloc(igs6_num*sizeof(int));
 if (!lines) error("malloc failed");
 for (i=0;i<igs6_num;i++) lines[i] = 0;
 fseek(f, 0, SEEK_SET);
 while (readline(f, line)!=EOF)
   {
    llen = strlen(line);
    if (llen!=81 && llen!=80) error("BAD IGS file, linelen = %d", strlen(line));
    if (line[72] == 'D')
      {
       if (sscanf(line, "%d", &entry)<1) error("scanning entry error");
       if (entry==126)
         {
	  if (ping) igs6_array[np].Did = atoi(line+73);
	  if (!ping) np++;
	  ping = !ping;
         }
      }
    if (line[72] == 'P')
    	{
	 pmid = param_id(line);
	 if (igs6_array[nr].Did==pmid ||
	    (igs6_array[nr].Did!=pmid && nr<igs6_num-1 && igs6_array[nr+1].Did==pmid))
	      {
	       if (igs6_array[nr].Did!=pmid) nr++;
	       lines[nr]++;
	      }
	 else if (nr==igs6_num) { debug("NOT intrested in rest DATA\n"); goto readprm; }
	}
   }
 readprm:
 for (i=0;i<igs6_num;i++) debug("lines[%d] = %d\n", i, lines[i]);
 nlines = (int*)malloc(igs6_num*sizeof(int));
 if (!nlines) error("malloc failed");
 for (i=0;i<igs6_num;i++) nlines[i] = 0;
 param_data = (char***)malloc(igs6_num*sizeof(char**));
 if (!param_data) error("malloc char***");
 for (i=0;i<igs6_num;i++)
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
	 if (igs6_array[nr].Did==pmid ||
	    (igs6_array[nr].Did!=pmid && nr<igs6_num-1 && igs6_array[nr+1].Did==pmid))
	      {
	       if (igs6_array[nr].Did!=pmid) nr++;
	       strcpy(param_data[nr][nlines[nr]], line);
	       nlines[nr]++;
	      }
	 else if (nr==igs6_num) { debug("IGS READ successfully.\n"); goto rsdis; }
	}
   }
 rsdis:
 nargs = (int*)malloc(igs6_num*sizeof(int));
 if (!nargs) error("malloc int*");
 for (i=0;i<igs6_num;i++) nargs[i] = 0;
 for (i=0;i<igs6_num;i++)
    for (j=0;j<lines[i];j++) nargs[i] += get_word_count(param_data[i][j]);
 for (i=0;i<igs6_num;i++) debug("nargs[%d] = %d\n", i, nargs[i]);
 args = (double**)malloc(igs6_num*sizeof(double*));
 if (!args) error("malloc double**");
 for (i=0;i<igs6_num;i++)
   {
    args[i] = (double*)malloc(nargs[i]*sizeof(double));
    if (!args[i]) error("malloc int*");
   }
 for (i=0;i<igs6_num;i++)
   {
    fill_args(args[i], param_data[i], lines[i], nargs[i]);
    generate_structures_126(args[i], nargs[i], i);
   }
 for (i=0;i<igs6_num;i++)
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
 char line[82];
 int n;
 int n126;
 int entry;
 int llen;
 if (interval<=0 || interval>MAX_INTERVAL) error("bad value: interval: %d", interval);
 debug("load_IGSfile: %s", fpath);
 f = fopen(fpath,"r");
 if (!f) error("cannot open file: %s", fpath);
 n=0;
 n126=0;
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
       if (entry==126) n126++;
      }
    if (line[72] == 'P') goto parameter;
   }
 error("No PARAMETER SECTION found.");
 parameter:
 debug("N128=%d\n", n);
 debug("N126=%d\n", n126);
 if (n<=0 && n126<=0) error("no 128 and 126 entires found, try with -E option");
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
 if (n126)
   {
    if (n126%2) error("bad lines count");
    n126 /= 2;
    fseek(f, 0, SEEK_SET);
    igs6_array = (struct IGS6Object*)malloc(n126*sizeof(struct IGS6Object));
    if (!igs6_array) error("malloc igs_array");
    igs6_num = n126;
    debug("SIZEOF IGS6Object is %d\n", sizeof(struct IGS6Object));
    load_126(f,n);
   }
 debug("IGS file read.");
}


void construct_tool(double r, double h)
{
 double kat;
 double dx = 0.0;
 double dy = 0.0;
 double dz = 0.0;
 int i;
 debug("construct_tool");
 debug("Scanned tool size is: R=%f, H=%f\n", r, h);
 mill = (double*)malloc(MILL_POINTS*sizeof(double));
 if (!mill) error("malloc mill structure");
 mill[0] = dx;
 mill[1] = dy;
 mill[2] = dz;
 mill[3] = dx;
 mill[4] = dy+h;
 mill[5] = dz;
 for (i=0;i<36;i++)
   {
    kat = (double)i*10.0;
    mill[6+(3*i)]   = cos(kat*PI/180.0)*r+dx;
    mill[6+(3*i)+1] = dy;
    mill[6+(3*i)+2] = sin(kat*PI/180.0)*r+dz;
   }
 for (i=0;i<36;i++)
   {
    kat = (double)i*10.0;
    mill[114+(3*i)]   = cos(kat*PI/180.0)*r+dx;
    mill[114+(3*i)+1] = dy+h;
    mill[114+(3*i)+2] = sin(kat*PI/180.0)*r+dz;
   }
 mill_buff = (double*)malloc(MILL_POINTS*sizeof(double));
 if (!mill_buff) error("malloc mill_buff");
}


void rotate_tool()
{
 int i;
 double x,y;
 debug("rotate_tool");
 debug("rotations: %d,%d,%d\n", rotate_tool_x, rotate_tool_y, rotate_tool_z);
 for (i=0;i<74;i++)
   {
    y             = mill[3*i+1]*cosines[rotate_tool_x] - mill[3*i+2]*sines[rotate_tool_x];
    mill[3*i+2] = mill[3*i+1]*sines[rotate_tool_x] + mill[3*i+2]*cosines[rotate_tool_x];
    mill[3*i+1] = y;
   }
 for (i=0;i<74;i++)
   {
    x             = mill[3*i]*cosines[rotate_tool_y] - mill[3*i+2]*sines[rotate_tool_y];
    mill[3*i+2] = mill[3*i]*sines[rotate_tool_y]   + mill[3*i+2]*cosines[rotate_tool_y];
    mill[3*i]   = x;
   }
 for (i=0;i<74;i++)
   {
    x             = mill[3*i]*cosines[rotate_tool_z] - mill[3*i+1]*sines[rotate_tool_z];
    mill[3*i+1] = mill[3*i]*sines[rotate_tool_z]   + mill[3*i+1]*cosines[rotate_tool_z];
    mill[3*i] = x;
   }
}


void load_tool(char* fpath)
{
 FILE* f;
 char line[MAX_LINE];
 int found;
 double d1,d2,d3,d4,d5;
 debug("load_tool");
 found = 0;
 f = fopen(fpath,"r");
 if (!f) error("cannot open file: %s", fpath);
 strcpy(line,"");
 while (readline(f, line)!=EOF)
   if ((sscanf(line,"TLDATA/MILL,%lf,%lf,%lf,%lf,%lf",&d1,&d2,&d3,&d4,&d5))==5)
     { found=1; break; }
 if (!found) error("tool MILL not found");
 fclose(f);
 construct_tool((double)d1,(double)d3);
 rotate_tool();
}


void load_camfile(char* fpath)
{
 FILE* f;
 char line[MAX_LINE];
 int n;
 double d1,d2,d3;
 double ix,ax,iy,ay,iz,az;
 debug("load_camfile");
 debug("loading CAM file: %s", fpath);
 load_tool(fpath);
 f = fopen(fpath,"r");
 if (!f) error("cannot open file: %s", fpath);
 n=0;
 strcpy(line,"");
 debug("determining how much lines read");
 while (readline(f, line)!=EOF) if ((sscanf(line,"GOTO/%lf,%lf,%lf",&d1,&d2,&d3))==3) n++;
 if (n<=0) error("no GOTO entires found.");
 npoints = n;
 path = (struct Point*)malloc(n*sizeof(struct Point));
 if (!path) error("malloc path");
 path_buff = (struct Point*)malloc(n*sizeof(struct Point));
 if (!path_buff) error("malloc path_buff");
 fseek(f,0,SEEK_SET);
 strcpy(line,"");
 n=-1;
 ix=iy=iz =  1e10;
 ax=ay=az = -1e10;
 while (readline(f, line)!=EOF) if ((sscanf(line,"GOTO/%lf,%lf,%lf",&d1,&d2,&d3))==3)
   {
    n++;
    path[n].x = d1 + XOFFSET;
    path[n].y = d2 + YOFFSET;
    path[n].z = d3 + ZOFFSET;
    if (path[n].x > ax) ax = path[n].x;
    if (path[n].y > ay) ay = path[n].y;
    if (path[n].z > az) az = path[n].z;
    if (path[n].x < ix) ix = path[n].x;
    if (path[n].y < iy) iy = path[n].y;
    if (path[n].z < iz) iz = path[n].z;
   }
 debug("Read %d GOTO entires.\n",n+1);
 debug("MIN = (%f,%f,%f), MAX = (%f,%f,%f)\n", ix,iy,iz,ax,ay,az);
 debug("THESE ARE MIN/MAXES BY AXIS ONLY, NOT BY 3D POINT.\n");
 fclose(f);
 debug("CAM file read.");
}


void init_mutex()
{
 pthread_mutexattr_t mutex_attr;
 debug("inintalizing MUTEX routine");
 if (pthread_mutexattr_init(&mutex_attr)) error("creating MUTEX attributes.\n");
 if (pthread_mutex_init(&mutex, (const pthread_mutexattr_t*)(&mutex_attr))) error("createing MUTEX");
}


void send_signal()
{
 XEvent an_event;
 int err;
 debug("TIMER_THREAD, sending refresh signal");
 an_event.xclient.type = Expose;
 if (!Xup) { pthread_yield(); return; }
 debug("TIMER_THREAD: waiting on MUTEX");
 WAIT_MUTEX
 debug("TIMER_THREAD: in CRITICAL SECTION, sending ASYNC signal");
 err=XSendEvent(dsp,win,False,0,&an_event);
 XFlush(dsp);
 debug("TIMER_THREAD: freeing MUTEX...");
 FREE_MUTEX
 debug("TIMER_PROC: MUTEX freed.");
}


void* timer_thr_proc(void* dummy)
{
 debug("TIMER_THREAD: Started.");
 while (!done)
   {
    while (stopped) { usleep(MIN_TIMER); pthread_yield(); }
    usleep(TIMERINT);
    time_async+=LINESINC;
    send_signal();
    if (time_async>=npoints) { printf("All lines displayed, timer THREAD finished.\n"); return NULL; }
   }
 return NULL;
}


void create_async_timer()
{
 debug("creating asynchronous timer thread");
 Xup = 0;
 done = 0;
 stopped = 0;
 init_mutex();
 pthread_create(&timer_thr,NULL,timer_thr_proc,NULL);
}


int error_handler(Display*d,  XErrorEvent* xer)
{
 debug("error processing ASYNC signal");
 error("XServer, ASYNC ERROR, PTHREADS ARE GETTING ANGRY.");
 return 0;
}


void set_error_handler()
{
 debug("setting xserver error handler");
 XSetErrorHandler(error_handler);
}


void create_lights()
{
 debug("create_lights");
 lx = LX;
 ly = LY;
 lz = LZ;
}


void set_vars()
{
 debug("set_vars");
 rotate_tool_x = TX_STEP;
 rotate_tool_y = TY_STEP;
 rotate_tool_z = TZ_STEP;
}


void offset_igs()
{
 debug("offset_igs");
 int i,j;
 for (i=0;i<igs_num;i++)
 for (j=0;j<ntr[i];j++)
   {
    if ((j%3)==0) tr[i][j] += XOFFSET;
    if ((j%3)==1) tr[i][j] += YOFFSET;
    if ((j%3)==2) tr[i][j] += ZOFFSET;
   }
}


void offset_igs6()
{
 debug("offset_igs6");
 int i,j;
 for (i=0;i<igs6_num;i++)
 for (j=0;j<npt[i];j++)
   {
    if ((j%3)==0) pt[i][j] += XOFFSET;
    if ((j%3)==1) pt[i][j] += YOFFSET;
    if ((j%3)==2) pt[i][j] += ZOFFSET;
   }
}


void generate_bspline_buffers()
{
 int i,j;
 debug("generate_bspline_buffers");
 if (igs_num>0)
   {
    ntr = (int*)malloc(igs_num*sizeof(int));
    if (!ntr) error("malloc ntr");
    tr = (double**)malloc(igs_num*sizeof(double*));
    if (!tr) error("malloc tr");
    tr_buff = (double**)malloc(igs_num*sizeof(double*));
    if (!tr_buff) error("malloc tr");
    for (i=0;i<igs_num;i++)
      {
       ntr[i] = igs_array[i].ntriangles;
       debug("ntr[%d] = %d\n", i, ntr[i]);
      }
    for (i=0;i<igs_num;i++)
      {
       tr[i] = (double*)malloc(ntr[i]*sizeof(double));
       if (!tr) error("malloc tr_buff");
       for (j=0;j<ntr[i];j++) tr[i][j] = igs_array[i].triangles[j];
       tr_buff[i] = (double*)malloc(ntr[i]*sizeof(double));
       if (!tr_buff) error("malloc tr_buff");
       for (j=0;j<ntr[i];j++) tr_buff[i][j] = tr[i][j];
      }
    offset_igs();
    for (i=0;i<igs_num;i++)
      {
       vdebug("OBJECT #%d\n", i);
       for (j=0;j<ntr[i]/3;j++)
        {
         vdebug("[%f,%f,%f] - ", tr[i][3*j], tr[i][3*j+1], tr[i][3*j+2]);
         vdebug("[%f,%f,%f]\n", tr_buff[i][3*j], tr_buff[i][3*j+1], tr_buff[i][3*j+2]);
        }
      }
   }
 else { ntr = 0; tr = 0; tr_buff = 0; }
 if (igs6_num>0)
   {
    npt = (int*)malloc(igs6_num*sizeof(int));
    if (!npt) error("malloc npt");
    pt = (double**)malloc(igs6_num*sizeof(double*));
    if (!pt) error("malloc pt");
    pt_buff = (double**)malloc(igs6_num*sizeof(double*));
    if (!pt_buff) error("malloc tr");
    for (i=0;i<igs6_num;i++)
      {
       npt[i] = igs6_array[i].npoints;
       debug("npt[%d] = %d\n", i, npt[i]);
      }
    for (i=0;i<igs6_num;i++)
      {
       pt[i] = (double*)malloc(npt[i]*sizeof(double));
       if (!pt) error("malloc tr_buff");
       for (j=0;j<npt[i];j++) pt[i][j] = igs6_array[i].points[j];
       pt_buff[i] = (double*)malloc(npt[i]*sizeof(double));
       if (!pt_buff) error("malloc pt_buff");
       for (j=0;j<npt[i];j++) pt_buff[i][j] = pt[i][j];
      }
    offset_igs6();
   }
 else { npt = 0; pt = 0; pt_buff = 0; }
 debug("structures_IGS: destoyed");
}


void free_igsobj(struct IGSObject* obj)
{
 debug("free_IGS_obj");
 int i;
 if (!obj) return;
 free(obj->triangles);
 free(obj->T);
 free(obj->S);
 for (i=0;i<=obj->K1;i++) { free(obj->P[i]); free(obj->W[i]); }
 free(obj->P);
 free(obj->W);
 obj = NULL;
}


void free_igs6obj(struct IGS6Object* obj)
{
 debug("free_IGS6_obj");
 if (!obj) return;
 free(obj->points);
 free(obj->T);
 free(obj->P);
 free(obj->W);
 obj = NULL;
}


void free_temporary_structures()
{
 debug("free_temporary_structures");
 int i;
 if (igs_num>0)
   {
    if (!igs_array) goto free_next;
    for (i=0;i<igs_num;i++) free_igsobj(&igs_array[i]);
    free(igs_array);
    igs_array = NULL;
   }
 free_next:
 if (igs6_num>0)
   {
    if (!igs6_array) return;
    for (i=0;i<igs6_num;i++) free_igs6obj(&igs6_array[i]);
    free(igs6_array);
    igs6_array = NULL;
   }
}


void X(char* camfile, char* igsfile)
{
 int s_num,i;
 XTextProperty window_name_property;
 XEvent an_event;
 char* window_name;
 debug("initializing XWindow");
 set_vars();
 create_f_table();
 if (want_cam) load_camfile(camfile);
 if (want_igs) load_igsfile(igsfile);
 create_lights();
 if (want_igs) generate_bspline_buffers();
 if (want_igs) free_temporary_structures();
 copy_to_buffers();
 zbuff=NULL;
 if (!want_triangles)
   {
    create_zbuff();
    clear_zbuff();
   }
 if (want_cam) create_async_timer();
 info_out();
 debug("initialize complete, creating window.");
 time_async=0;
 if ((dsp=XOpenDisplay(NULL))==NULL) error("connect to X");
 s_num = DefaultScreen(dsp);
 win = XCreateSimpleWindow(dsp, RootWindow(dsp, s_num),0, 0, WX, WY, 1,WhitePixel(dsp, s_num),BlackPixel(dsp, s_num));
 XMapWindow(dsp, win);
 XFlush(dsp);
 gc = XCreateGC(dsp, win, 0, NULL);
 if ((int)gc<0) error("create GC");
 XSetForeground(dsp, gc, WhitePixel(dsp, s_num));
 XSetBackground(dsp, gc, BlackPixel(dsp, s_num));
 XSelectInput(dsp, win, ExposureMask | KeyPressMask | StructureNotifyMask);
 window_name = (char*)malloc(0x100);
 if (!window_name) error("malloc window_name");
 strcpy(window_name, "CAD/CAM Zaliczenie - Lukasz Gryglicki M1");
 if (XStringListToTextProperty(&window_name,1,&window_name_property))
    XSetWMName(dsp, win, &window_name_property);
 free(window_name);
 debug("window created.");
 set_error_handler();
 refresh_window(0);
 Xup=1;
 while (!done)
    {
     debug("waiting for X event...");
     XNextEvent(dsp, &an_event);
     debug("X event occured: (%d), waiting on MUTEX",an_event.type);
     WAIT_MUTEX
     debug("X_THREAD in CRITICAL SECTION");
     switch (an_event.type)
        {
         case Expose:
             debug("Expose event");
	     if (an_event.xexpose.send_event)
	       {
	        refresh_window(1);
	        refresh_window(0);
	       }
	     else refresh_window(0);
             break;
        case ConfigureNotify:
             debug("Configure event");
/*             XResizeWindow(dsp,win,0x200,0x200);*/
	     if (!want_triangles) delete_zbuff();
             WX = an_event.xconfigure.width;
             WY = an_event.xconfigure.height;
	     if (!want_triangles)
	       {
	        create_zbuff();
	        clear_zbuff();
	       }
	     refresh_window(1);
	     refresh_window(0);
             break;
        case KeyPress:
             debug("KeyPress event");
	     debug("Angles (%d,%d,%d), offsets (%f,%f,%f), light (%f,%f,%f)\n", angle_x,angle_y,angle_z,ox,oy,oz,lx,ly,lz);
	     debug("KeyCode: 0x%x", an_event.xkey.keycode);
	     if  (an_event.xkey.keycode == KEY_Q) done=1;
	     if  (an_event.xkey.keycode == KEY_R)
	       {
		want_triangles = !want_triangles;
		if (!want_triangles)			//INVERTED CONDICTION 
		  {
		   create_zbuff();
		   clear_zbuff();
		  }
		else delete_zbuff();
	        refresh_window(1);
	        refresh_window(0);
	       }
	     if  (an_event.xkey.keycode == KEY_S) stopped = !stopped;
	     if  (an_event.xkey.keycode == KEY_UP)    { refresh_window(1); angle_x+=AX_STEP; refresh_window(0); }
	     if  (an_event.xkey.keycode == KEY_DOWN)  { refresh_window(1); angle_x-=AX_STEP; refresh_window(0); }
	     if  (an_event.xkey.keycode == KEY_LEFT)  { refresh_window(1); angle_y+=AY_STEP; refresh_window(0); }
	     if  (an_event.xkey.keycode == KEY_RIGHT) { refresh_window(1); angle_y-=AY_STEP; refresh_window(0); }
	     if  (an_event.xkey.keycode == KEY_PGUP)  { refresh_window(1); angle_z+=AZ_STEP; refresh_window(0); }
	     if  (an_event.xkey.keycode == KEY_PGDWN) { refresh_window(1); angle_z-=AZ_STEP; refresh_window(0); }
	     if (an_event.xkey.keycode == KEY_O )
                { refresh_window(1);  lz+=OINT; lzb=lz; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_M )
                { refresh_window(1);  lz-=OINT; lzb=lz; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_J )
                { refresh_window(1);  lx-=OINT; lxb=lx; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_L )
                { refresh_window(1);  lx+=OINT; lxb=lx; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_I )
                { refresh_window(1);  ly+=OINT; lyb=ly; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_K )
                { refresh_window(1);  ly-=OINT; lyb=ly; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_E )
                { refresh_window(1);  oz+=OINT; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_Z )
                { refresh_window(1);  oz-=OINT; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_A )
                { refresh_window(1);  ox-=OINT; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_D )
                { refresh_window(1);  ox+=OINT; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_W )
                { refresh_window(1);  oy+=OINT; refresh_window(0);}
	     if (an_event.xkey.keycode == KEY_S )
                { refresh_window(1);  oy-=OINT; refresh_window(0);}
	     if  (an_event.xkey.keycode == KEY_SPC  )
	       {
		refresh_window(1);
		angle_x = angle_y = angle_z = 0;
		ox = oy = oz = 0.0;
		lx = ly = 0.0;
		lz = -3.0;
		refresh_window(0);
	       }
             break;
        default:
	     debug("Unknown or Unhandled event, default action token.");
             break;
        }
     debug("X_THREAD: freeing MUTEX...");
     FREE_MUTEX
     debug("X_THREAD: MUTEX freed.");
   }
 debug("Waiting for TIMER_THREAD to finish...");
 if (!want_triangles)  delete_zbuff();
 if (want_cam)
   {
    debug("destroying CAM structures");
    if (!stopped)
      pthread_join(timer_thr, NULL);
    if (path) free(path);
    path = 0;
    if (path_buff) free(path_buff);
    path_buff = 0;
    if (mill) free(mill);
    mill = 0;
    if (mill_buff) free(mill_buff);
    mill_buff = 0;
   }
 if (want_igs)
   {
    debug("destroying IGS structures");
    if (igs_array) free(igs_array);
    igs_array = NULL;
    if (igs6_array) free(igs6_array);
    igs6_array = NULL;
    if (tr)
      {
       for (i=0;i<igs_num;i++) if (tr[i]) free(tr[i]);
       free(tr);
       for (i=0;i<igs_num;i++) if (tr_buff[i]) free(tr_buff[i]);
       free(tr_buff);
      }
    if (pt)
      {
       for (i=0;i<igs6_num;i++) if (pt[i]) free(pt[i]);
       free(pt);
       for (i=0;i<igs6_num;i++) if (pt_buff[i]) free(pt_buff[i]);
       free(pt_buff);
      }
   }
 debug("X-Main Loop finished, TIMER_THREAD stopped.");
 debug("MEMORY freed");
}


void help()
{
 debug("help");
 printf("xcam is a program displaying CAM menufacturing file\nAvailable options are: (if not written the integer is expected)\n");
 printf("\t-A Win_X_Size (DEFAULT=512)\n");
 printf("\t-S Win_Y-Sizxe (DEFAULT=512)\n");
 printf("\t-X Start_X_Angle (DEFAULT=0)\n");
 printf("\t-Y Start_Y_Angle (DEFAULT=0)\n");
 printf("\t-Z Start_Z_Angle (DEFAULT=0)\n");
 printf("\t-i Step_X_Angle (DEFAULT=15)\n");
 printf("\t-j Step_Y_Angle (DEFAULT=15)\n");
 printf("\t-k Step_Z_Angle (DEFAULT=15)\n");
 printf("\t-I Rotate_X_Tool (DEFAULT=0)\n");
 printf("\t-J Rotate_Y_Tool (DEFAULT=0)\n");
 printf("\t-K Rotate_Z_Tool (DEFAULT=0)\n");
 printf("\t-x X_Offset (DEFAULT=0.0, double value)\n");
 printf("\t-y Y_Offset (DEFAULT=0.0, double value)\n");
 printf("\t-z Z_Offset (DEFAULT=0.0, double value)\n");
 printf("\t-m Max_Lines (DEFAULT=1024)\n");
 printf("\t-t Timer_Interval (DEFAULT=90000, in miliseconds)\n");
 printf("\t-l How_Much_Lines_Add_Per_Timer_Tick (DEFAULT=10)\n");
 printf("\t-a Stretch_X (DEFAULT=1.41, resize x by 1.41)\n");
 printf("\t-b Stretch_Y (DEFAULT=1.41, resize y by 1.41)\n");
 printf("\t-o Distance (DEFAULT=12.0, double value)\n");
 printf("\t-O Object_move_step (DEFAULT=0.2, double value)\n");
 printf("\t-L LightDiff (DEFAULT=10, between dark and light)\n");
 printf("\t-r Triangle_Interval (DEFAULT=10, positive)\n");
 printf("\t-C LineColor (DEFAULT=00:FF:00, format HH:HH:HH, H=hex_digit)\n");
 printf("\t-U MillColor (DEFAULT=96:64:0A, format HH:HH:HH, H=hex_digit)\n");
 printf("\t-V SurfColor (DEFAULT=32:32:DC, format HH:HH:HH, H=hex_digit)\n");
 printf("\t-W CurvColor (DEFAULT=32:FA:78, format HH:HH:HH, H=hex_digit)\n");
 printf("\t-P LightPos (DEFAULT=0.0:0.0:-3.0 , format lf:lf:lf, lf=double)\n");
 printf("\t-f CAM_file(DEFAULT='camfile.dat', string value)\n");
 printf("\t-g IGS_file(DEFAULT='igsfile.igs', string value)\n");
 printf("\t-M Do not use CAM file\n");
 printf("\t-E Do not use IGS file\n");
 printf("\t-d (enables debug)\n");
 printf("\t-h (displays this help message)\n");
 printf("\t-w enables RAW triangles and disables ZBUFF\n");
 printf("Controls:\n");
 printf("\ts -stop/start animation\n");
 printf("\tr -toggle RAW_TRIANGLES\n");
 printf("\tq -quit animation\n");
 printf("\tUP/DOWN ARROW    +/- ANGLE_X\n");
 printf("\tLEFT/RIGHT ARROW +/- ANGLE_Y\n");
 printf("\tPGUP/PGDOWN      +/- ANGLE_Z\n");
 printf("\tLIGHT  MOVING: X,Y,Z: J:L,K:I,M:O\n");
 printf("\tOBJECT MOVING: X,Y,Z: A:D,S:W,Z:E\n");
 printf("Written by MorgothDBMA, morgothdbma@o2.pl, +48693582014\n");
 printf("Lukasz Gryglicki MiNI M1\n\t2004: %s:%d, %s %s\n",__FILE__,__LINE__,__DATE__,__TIME__);
 debug("help done.");
}


int main(int lb, char** par)
{
 char u;
 char camfile[MAXPATH+1];
 char igsfile[MAXPATH+1];
 char linecol[0x100];
 char millcol[0x100];
 char igscol[0x100];
 char igs6col[0x100];
 char lightpos[0x100];
 printf("NIE DO WIARY - PROGRAM DZIALA DOBRZE ALE Z DEBUG JUZ NIE!!!\n");
 printf("BLAD SEGMENTACJI !!!! PAMIECI !!!! CZAI SIE GDZIES WEWNATRZ ----\n");
 printf("THIS PROGRAM IS LUCKY OWNER OF MEMORY SEGMENTATION ERROR SOMEWHERE\n"
		 "THIS ERROR APPEARS ONLY IN SPECIAL DEBUG MODE\n"
		 "BUT CAN CAUSE SOME IRRATIONAL BEHAVIOUR OF PROGRAM\n"
		 "*** YOU WERE WARNED! ***\n");
 strcpy(camfile,"camfile.dat");
 strcpy(igsfile,"igsfile.igs");
 strcpy(linecol, "00:FF:00");
 strcpy(millcol, "96:64:0A");
 strcpy(igscol,  "32:32:DC");
 strcpy(igs6col, "32:FA:78");
 debug("main::getopt");
 while ((u = getopt(lb,par,"A:S:X:i:Y:j:Z:k:x:y:z:m:t:l:a:b:f:I:J:K:o:L:C:g:P:O:r:iU:V:W:EMhwd"))!=-1)
 {
  switch (u)
   {
    case 'A': WX = atoi(optarg); break;
    case 'S': WY = atoi(optarg); break;
    case 'X': START_AX = atoi(optarg); break;
    case 'Y': START_AY = atoi(optarg); break;
    case 'Z': START_AZ = atoi(optarg); break;
    case 'i': AX_STEP  = atoi(optarg); break;
    case 'j': AY_STEP  = atoi(optarg); break;
    case 'k': AZ_STEP  = atoi(optarg); break;
    case 'I': TX_STEP  = atoi(optarg); break;
    case 'J': TY_STEP  = atoi(optarg); break;
    case 'K': TZ_STEP  = atoi(optarg); break;
    case 'x': XOFFSET  = atof(optarg); break;
    case 'y': YOFFSET  = atof(optarg); break;
    case 'z': ZOFFSET  = atof(optarg); break;
    case 'o': zoff     = atof(optarg); break;
    case 'O': OINT     = atof(optarg); break;
    case 'm': MAXLINES = atoi(optarg); break;
    case 't': TIMERINT = atoi(optarg); break;
    case 'l': LINESINC = atoi(optarg); break;
    case 'a': SCALEX   = atoi(optarg); break;
    case 'b': SCALEY   = atoi(optarg); break;
    case 'f': if (strlen(optarg)<MAXPATH) strcpy(camfile, optarg); break;
    case 'g': if (strlen(optarg)<MAXPATH) strcpy(igsfile, optarg); break;
    case 'L': ldiff    = atoi(optarg); break;
    case 'r': interval = atoi(optarg); break;
    case 'C': if (strlen(optarg)<0x40) strcpy(linecol,  optarg); break;
    case 'U': if (strlen(optarg)<0x40) strcpy(millcol,  optarg); break;
    case 'V': if (strlen(optarg)<0x40) strcpy(igscol,   optarg); break;
    case 'W': if (strlen(optarg)<0x40) strcpy(igs6col,  optarg); break;
    case 'P': if (strlen(optarg)<0x80) strcpy(lightpos, optarg); break;
    case 'w': want_triangles=1; break;
    case 'h': help(); return 0;
    case 'd': dbg++; break;
    case 'M': want_cam = 0; break;
    case 'E': want_igs = 0; break;
    default: printf("%s: Unrecognized option\n",par[0]); return 1;
   }
 }
 line_r =    0;
 line_g = 0xFF;
 line_b =    0;
 mill_r =  150;
 mill_g =  100;
 mill_b =   10;
 igs_r  =   50;
 igs_g  =   50;
 igs_b  =  220;
 igs6_r =   50;
 igs6_g =  250;
 igs6_b =  120;
 LX = -0.0;
 LY = -0.0;
 LZ = -3.0;
 sscanf(linecol, "%02x:%02x:%02x", &line_r,&line_g,&line_b);
 sscanf(millcol, "%02x:%02x:%02x", &mill_r,&mill_g,&mill_b);
 sscanf(igscol, "%02x:%02x:%02x",  &igs_r,&igs_g,&igs_b);
 sscanf(igs6col, "%02x:%02x:%02x", &igs6_r,&igs6_g,&igs6_b);
 sscanf(lightpos,"%lf:%lf:%lf", &LX,&LY,&LZ);
 X(camfile, igsfile);
 debug("exiting MAIN.");
 return 0;
}

