#include <stdlib.h>
#include <stdio.h>

void table(char* cubecmd, int sidx, int tid1, int tid2, int tid3, double big)
{
 char syscmd[1024];
 double dx,dy,dz,rx,ry,rz;
 dx = dy = dz = rx = ry = rz = 0.;
 if (!cubecmd) return;
 /* list transforms */
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 1",
	 cubecmd, sidx, tid1, big, big*10., big, rx,ry,rz, -7.5*big+dx, dy, -7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 "
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 1",
	 cubecmd, sidx+12, tid1, big, big*10., big, rx,ry,rz, 7.5*big+dx, dy, -7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 1",
	 cubecmd, sidx+24, tid1, big, big*10., big, rx,ry,rz, -7.5*big+dx, dy, 7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 1",
	 cubecmd, sidx+36, tid1, big, big*10., big, rx,ry,rz, 7.5*big+dx, dy, 7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.4 0.4 0.4 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 100 %f %f %f %f %f %f %f %f %f 0 1 1",
	 cubecmd, sidx+48, tid2, big*10., big, big*10., rx,ry,rz, dx, big*10.+dy, dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.85 0.85 0.75 0.2 0.2 0.5 0.4 0.4 0.5 " 
	 "%d 1.2 1.4 1.6 150 %f %f %f %f %f %f %f %f %f 0 1 1",
	 cubecmd, sidx+60, tid3, big*10.5, big*.5, big*10.5, rx,ry,rz, dx, big*11.+dy, dz);
 system(syscmd);
 /* triangles */
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 0",
	 cubecmd, sidx, tid1, big, big*10., big, rx,ry,rz, -7.5*big+dx, dy, -7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 0",
	 cubecmd, sidx+12, tid1, big, big*10., big, rx,ry,rz, 7.5*big+dx, dy, -7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 0",
	 cubecmd, sidx+24, tid1, big, big*10., big, rx,ry,rz, -7.5*big+dx, dy, 7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.1 0.1 0.1 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 60 %f %f %f %f %f %f %f %f %f 0 1 0",
	 cubecmd, sidx+36, tid1, big, big*10., big, rx,ry,rz, 7.5*big+dx, dy, 7.5*big+dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.0 0.0 0.0 0.4 0.4 0.4 1.0 1.0 1.0 " 
	 "%d 1.2 1.2 1.2 100 %f %f %f %f %f %f %f %f %f 0 1 0",
	 cubecmd, sidx+48, tid2, big*10., big, big*10., rx,ry,rz, dx, big*10.+dy, dz);
 system(syscmd);
 sprintf(syscmd, "%s %d 0.85 0.85 0.75 0.2 0.2 0.5 0.4 0.4 0.5 " 
	 "%d 1.2 1.4 1.6 150 %f %f %f %f %f %f %f %f %f 0 1 0",
	 cubecmd, sidx+60, tid3, big*10.5, big*.5, big*10.5, rx,ry,rz, dx, big*11.+dy, dz);
 system(syscmd);
}

void help()
{
 printf("params: idx, texid(123), size, cubecmd\n");
 printf("suggestion: texid(1,2,3) = 8 4 9\n");
}

int main(int lb, char** par)
{
 int idx, tid1, tid2, tid3;
 double big;
 if (lb != 7) { help(); exit(1); }
 idx = atoi(par[1]);
 tid1 = atoi(par[2]);
 tid2 = atoi(par[3]);
 tid3 = atoi(par[4]);
 big = atof(par[5]);
 table(par[6], idx, tid1, tid2, tid3, big);
 return 0;
}
