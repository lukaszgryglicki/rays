#include <stdlib.h>
#include <stdio.h>
void cube_gen(int s_idx, 
double t_r, double t_g, double t_b, double s_r, double s_g, double s_b, double c_r, double c_g, double c_b, int t_id, double tr_r, double tr_g, double tr_b, double s_fac,
double x_x, double y_y, double z_z, double r_x, double r_y, double r_z, double d_x, double d_y, double d_z, int inv_n, int faces, int lt)
{
 int i;
 if (lt == 1 || lt == 2)
 {   
 printf("ListTransform: [%d,%d]\n", s_idx, s_idx+11);
 printf("{\nTranslate: (%f,%f,%f)\n", d_x, d_y,d_z);
 printf("RotateX: %f\n", r_x);
 printf("RotateY: %f\n", r_y);
 printf("RotateZ: %f\n", r_z);
 printf("Scale: (%f,%f,%f)\n", x_x, y_y, z_z);
 if (inv_n) printf("NegateN:\n");
 printf("}\n");
 }
 if (lt == 0 || lt == 2)
 {
 for (i=s_idx;i<=s_idx+11;i++)
   {
printf("Triangle: %d\n", i);
printf("{\n");
if (((i-s_idx)%2) == 0)
 {
  printf(" a: Vertex: (-1,-1,1)\n");
  printf(" b: Vertex: (1,1,1)\n");
  printf(" c: Vertex: (1,-1,1)\n");
  printf(" texA: TexCoord: (0,0)\n");
  printf(" texB: TexCoord: (1,1)\n");
  printf(" texC: TexCoord: (1,0)\n");
 }
else
 {
  printf(" a: Vertex: (-1,-1,1)\n");
  printf(" b: Vertex: (-1,1,1)\n");
  printf(" c: Vertex: (1,1,1)\n");
  printf(" texA: TexCoord: (0,0)\n");
  printf(" texB: TexCoord: (0,1)\n");
  printf(" texC: TexCoord: (1,1)\n");
 }
printf(" na: Vector: (0,0,1)\n");
printf(" nb: Vector: (0,0,1)\n");
printf(" nc: Vector: (0,0,1)\n");
printf(" transparency: RGB: (%f,%f,%f)\n", t_r, t_g, t_b);
printf(" specular: RGB: (%f,%f,%f)\n", s_r, s_g, s_b);
printf(" diffuse: RGB: (%f,%f,%f)\n", c_r,c_g,c_b);
printf(" transparencyFactR: (1,%f)\n", tr_r);
printf(" transparencyFactG: (1,%f)\n", tr_g);
printf(" transparencyFactB: (1,%f)\n", tr_b);
printf(" specularFact: %f\n", s_fac);
printf(" faces: %d\n", faces);
printf(" texture: %d\n", t_id);
printf(" Transform:\n");
printf(" {\n");
switch (i-s_idx)
  {
   case 0: 
   case 1: 
       break;
   case 2: 
   case 3: 
  printf(" RotateY: 90\n");
       break;
   case 4: 
   case 5: 
   printf(" RotateY: 180\n");
       break;
   case 6: 
   case 7: 
  printf("  RotateY: 270\n");
       break;
   case 8: 
   case 9: 
  printf("  RotateX: 90\n");
       break;
   case 10: 
   case 11: 
  printf("  RotateX: 270\n");
       break;
   default: printf("Error. idx = %d, t_ct = %d\n", i, i-s_idx);
  }
printf(" }\n");
printf("}\n");
   }
 }
}
void help()
{
 printf("params: t_idx (start idx [t_idx,t_idx+11]\n");
 printf("trans(rgb), spec(rgb), color(rgb), texid, transFact(rgb), specFact\n");
 printf("scale(xyz), rotate(xyz), translate(xyz) invN, faces, lt\n");
 printf("lt = 0, only tr, lt = 1 only listtransf, lt=2 both\n");
}

int main(int lb, char** par)
{
 double t_r,t_g,t_b,s_r,s_g,s_b,c_r,c_g,c_b,tr_r,tr_g,tr_b,s_fac;
 double x_x,y_y,z_z,r_x,r_y,r_z,d_x,d_y,d_z;
 int s_idx, t_id, inv_n, faces, lt;
 if (lb != 28) { printf("lbpar = %d\n", lb); help(); exit(1); }
  s_idx = atoi(par[1]);
  t_r = atof(par[2]);
  t_g = atof(par[3]);
  t_b = atof(par[4]);
  s_r = atof(par[5]);
  s_g = atof(par[6]);
  s_b = atof(par[7]);
  c_r = atof(par[8]);
  c_g = atof(par[9]);
  c_b = atof(par[10]);
  t_id = atoi(par[11]);
  tr_r = atof(par[12]);
  tr_g = atof(par[13]);
  tr_b = atof(par[14]);
  s_fac = atof(par[15]);
  x_x = atof(par[16]);
  y_y = atof(par[17]);
  z_z = atof(par[18]);
  r_x = atof(par[19]);
  r_y= atof(par[20]);
  r_z = atof(par[21]);
  d_x = atof(par[22]);
  d_y = atof(par[23]);
  d_z = atof(par[24]);
  inv_n = atoi(par[25]);
  faces = atoi(par[26]);
  lt = atoi(par[27]);
 cube_gen(s_idx, 
 t_r, t_g, t_b, s_r, s_g, s_b, c_r, c_g, c_b, t_id, tr_r, tr_g, tr_b, s_fac,
 x_x, y_y, z_z, r_x, r_y, r_z, d_x, d_y, d_z, inv_n, faces, lt);
 return 0;
}
