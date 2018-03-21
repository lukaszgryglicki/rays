#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void outf(char* par)
{
  char str[1024];
  sprintf(str, "echo \"%s\" >> prog.dat", par);
  system( str );
}

int ran(int num)
{
    return rand() % num;
}

void mkball(int n, int lt)
{
 char str[1024];
 int nt, tid;
 double tr, tg, tb;
 double sr, sg, sb;
 double cr, cg, cb;
 double sfr, sfg, sfb;
 double sx, sy, sz, sf;
 double rx, ry, rz;
 double tx, ty, tz;

 nt = 2000 * n;
 tr = (double)ran(10);
 tg = (double)ran(10);
 tb = (double)ran(10);
 sr = (double)ran(10);
 sg = (double)ran(10);
 sb = (double)ran(10);
 cr = (double)ran(10);
 cg = (double)ran(10);
 cb = (double)ran(10);
 cr = 0.;
 sfr = 1.2 + (double)ran(20) / 100.;
 sfg = 1.2 + (double)ran(20) / 100.;
 sfb = 1.2 + (double)ran(20) / 100.;
 tid = ran(14) + 1;
 sx = 35. + (double)ran(25);
 sy = 35. + (double)ran(25);
 sz = 35. + (double)ran(25);
 rx = (double)ran(360);
 ry = (double)ran(360);
 rz = (double)ran(360);
 tx = (double)ran(250) - 125.;
 ty = (double)ran(250) - 125.;
 tz = (double)ran(250) - 125.;
 sf = 10. + (double)ran(90);

 sprintf(str, "../ball %d %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f %f 0 1 %d 25 20 >> prog.dat",
	 nt, tr, tr, tr, sr, sr, sr, cr, cr, cr, tid, sfr, sfg, sfb, sf, sx, sy, sz, rx, ry, rz, tx, ty, tz, lt);
 system( str );
}

void cptri(int n)
{
  char str[1024];
  sprintf(str, "CopyTriangles: dst=%d,src=0,num=2000", n * 2000);
  outf( str );
}

void go(int n)
{
  char str[1024];  
  time_t tm;
  int r, r2, r3, nt, i;
  time(&tm);
  srand((int)tm);
  system("rm -f prog.dat");
  outf("Screen: (1920,1200)");
  r = ran(14) + 1;
  sprintf(str, "Background: %d", r);
  outf( str );
  outf("MaxRecurse: 100");
  outf("Backup: 500");
  outf("Observer: Vertex: (0,0,-250)");
  r = 25 + ran(20);
  sprintf(str, "LookZ: %d%%", r);
  outf( str );
  r = ran(1000) - 500;
  r2 = ran(1000) - 500;
  r3 = ran(1000) - 500;
  sprintf(str, "Light: Vertex: (%d,%d,%d)", r, r2, r3);	
  outf( str );
  r = ran(0x100);
  r2 = ran(0x100);
  r3 = ran(0x100);
  sprintf(str, "LightColor: %02x%02x%02x", r, r2, r3);
  outf( str );
  outf("TexDirectory: tex");
  outf("NumTextures: 14");	
  nt = n * 2000;
  sprintf(str, "nTriangles: %d", nt);
  outf( str );
  sprintf(str, "ListTransform: [0,%d]", nt-1);
  outf( str );
  outf("{");
  outf("}");
  for (i=0;i<n;i++)
  {
      mkball(i, 1);
  }
  for (i=0;i<n;i++)
  {
      mkball(i, 0);
  }
/*  mkball(0, 0);*/
/*  for (i=1;i<n;i++) cptri(i);*/
}

int main(int lb, char** par)
{
    if (lb < 2)
    {
	printf("%s num\n", par[0]);
	return 1;
    }
    go(atoi(par[1]));
    return 0;
}

