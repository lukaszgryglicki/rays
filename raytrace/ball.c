#include <stdlib.h>
#include <stdio.h>
#include <math.h>
void ball_gen(int s_idx, 
double t_r, double t_g, double t_b, double s_r, double s_g, double s_b, double c_r, double c_g, double c_b, int t_id, double tr_r, double tr_g, double tr_b, double s_fac,
double x_x, double y_y, double z_z, double r_x, double r_y, double r_z, double d_x, double d_y, double d_z, int inv_n, int faces, int lt, int pa, int pb)
{
 int i,j;
 
 float alfa1, alfa2;
 float beta1, beta2;
 
 float fx11,fy11,fz11;
 float fx12,fy12,fz12;
 float fx21,fy21,fz21;
 float fx22,fy22,fz22;
 
 float tx1,ty1;
 float tx2,ty2;
 
 if (lt == 1 || lt == 2)
 {   
 printf("ListTransform: [%d,%d]\n", s_idx, s_idx+pa*pb*2-1);
 printf("{\nTranslate: (%f,%f,%f)\n", d_x, d_y,d_z);
 printf("RotateX: %f\n", r_x);
 printf("RotateY: %f\n", r_y);
 printf("RotateZ: %f\n", r_z);
 printf("Scale: (%f,%f,%f)\n", x_x, y_y, z_z);
 if (inv_n) printf("NegateN:\n");
 printf("}\n");
 printf("ListTransform: [%d,%d]\n", pa*pb*2+s_idx, pa*pb*2+s_idx+pa*pb*2-1);
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
 for (i=0;i<pa;i++)
 for (j=0;j<pb;j++)
 {
  alfa1 =  ((float)i / (float)pa) * 2. * M_PI;
  beta1 =  ((float)j / (float)pb) * .5 * M_PI;
  alfa2 = ((float)(i+1) / (float)pa) * 2. * M_PI;
  beta2 = ((float)(j+1) / (float)pb) * .5 * M_PI;

  fx11 = cos(alfa1) * cos(beta1);
  fy11 = sin(alfa1) * cos(beta1);
  fz11 = sin(beta1);

  fx12 = cos(alfa1) * cos(beta2);
  fy12 = sin(alfa1) * cos(beta2);
  fz12 = sin(beta2);

  fx21 = cos(alfa2) * cos(beta1);
  fy21 = sin(alfa2) * cos(beta1);
  fz21 = sin(beta1);

  fx22 = cos(alfa2) * cos(beta2);
  fy22 = sin(alfa2) * cos(beta2);
  fz22 = sin(beta2);

  tx1 = (float)i / (float)pa;
  ty1 = (float)j / (float)pb;
  tx2 = (float)(i+1) / (float)pa;
  ty2 = (float)(j+1) / (float)pb;
  
printf("Triangle: %d\n", s_idx+i*pa*2+2*j);
printf("{\n");
  printf(" a: Vertex: (%f,%f,%f)\n", fx11, fy11, fz11);
  printf(" b: Vertex: (%f,%f,%f)\n", fx21, fy21, fz21);
  printf(" c: Vertex: (%f,%f,%f)\n", fx22, fy22, fz22);
  printf(" texA: TexCoord: (%f,%f)\n",tx1,ty1);
  printf(" texB: TexCoord: (%f,%f)\n",tx2,ty1);
  printf(" texC: TexCoord: (%f,%f)\n",tx2,ty2);
printf(" na: Vector: (%f,%f,%f)\n", fx11, fy11, fz11);
printf(" nb: Vector: (%f,%f,%f)\n", fx21, fy21, fz21);
printf(" nc: Vector: (%f,%f,%f)\n", fx22, fy22, fz22);
printf(" transparency: RGB: (%f,%f,%f)\n", t_r, t_g, t_b);
printf(" specular: RGB: (%f,%f,%f)\n", s_r, s_g, s_b);
printf(" diffuse: RGB: (%f,%f,%f)\n", c_r,c_g,c_b);
printf(" transparencyFactR: (1,%f)\n", tr_r);
printf(" transparencyFactG: (1,%f)\n", tr_g);
printf(" transparencyFactB: (1,%f)\n", tr_b);
printf(" specularFact: %f\n", s_fac);
printf(" faces: %d\n", faces);
printf(" texture: %d\n", t_id);
printf("}\n");
printf("Triangle: %d\n", s_idx+i*pa*2+2*j+1);
printf("{\n");
  printf(" a: Vertex: (%f,%f,%f)\n", fx11, fy11, fz11);
  printf(" b: Vertex: (%f,%f,%f)\n", fx22, fy22, fz22);
  printf(" c: Vertex: (%f,%f,%f)\n", fx12, fy12, fz12);
  printf(" texA: TexCoord: (%f,%f)\n",tx1,ty1);
  printf(" texB: TexCoord: (%f,%f)\n",tx2,ty2);
  printf(" texC: TexCoord: (%f,%f)\n",tx1,ty2);
printf(" na: Vector: (%f,%f,%f)\n", fx11, fy11, fz11);
printf(" nb: Vector: (%f,%f,%f)\n", fx22, fy22, fz22);
printf(" nc: Vector: (%f,%f,%f)\n", fx12, fy12, fz12);
printf(" transparency: RGB: (%f,%f,%f)\n", t_r, t_g, t_b);
printf(" specular: RGB: (%f,%f,%f)\n", s_r, s_g, s_b);
printf(" diffuse: RGB: (%f,%f,%f)\n", c_r,c_g,c_b);
printf(" transparencyFactR: (1,%f)\n", tr_r);
printf(" transparencyFactG: (1,%f)\n", tr_g);
printf(" transparencyFactB: (1,%f)\n", tr_b);
printf(" specularFact: %f\n", s_fac);
printf(" faces: %d\n", faces);
printf(" texture: %d\n", t_id);
printf("}\n");
   }
 s_idx += pa*pb*2;
 for (i=0;i<pa;i++)
 for (j=0;j<pb;j++)
 {
  alfa1 =  ((float)i / (float)pa) * 2. * M_PI;
  beta1 =  ((float)j / (float)pb) * .5 * M_PI;
  alfa2 = ((float)(i+1) / (float)pa) * 2. * M_PI;
  beta2 = ((float)(j+1) / (float)pb) * .5 * M_PI;

  fx11 = cos(alfa1) * cos(beta1);
  fy11 = sin(alfa1) * cos(beta1);
  fz11 = -sin(beta1);

  fx12 = cos(alfa1) * cos(beta2);
  fy12 = sin(alfa1) * cos(beta2);
  fz12 = -sin(beta2);

  fx21 = cos(alfa2) * cos(beta1);
  fy21 = sin(alfa2) * cos(beta1);
  fz21 = -sin(beta1);

  fx22 = cos(alfa2) * cos(beta2);
  fy22 = sin(alfa2) * cos(beta2);
  fz22 = -sin(beta2);

  tx1 = (float)i / (float)pa;
  ty1 = (float)j / (float)pb;
  tx2 = (float)(i+1) / (float)pa;
  ty2 = (float)(j+1) / (float)pb;
  
printf("Triangle: %d\n", s_idx+i*pa*2+2*j);
printf("{\n");
  printf(" a: Vertex: (%f,%f,%f)\n", fx11, fy11, fz11);
  printf(" b: Vertex: (%f,%f,%f)\n", fx21, fy21, fz21);
  printf(" c: Vertex: (%f,%f,%f)\n", fx22, fy22, fz22);
  printf(" texA: TexCoord: (%f,%f)\n",tx1,ty1);
  printf(" texB: TexCoord: (%f,%f)\n",tx2,ty1);
  printf(" texC: TexCoord: (%f,%f)\n",tx2,ty2);
printf(" na: Vector: (%f,%f,%f)\n", fx11, fy11, fz11);
printf(" nb: Vector: (%f,%f,%f)\n", fx21, fy21, fz21);
printf(" nc: Vector: (%f,%f,%f)\n", fx22, fy22, fz22);
printf(" transparency: RGB: (%f,%f,%f)\n", t_r, t_g, t_b);
printf(" specular: RGB: (%f,%f,%f)\n", s_r, s_g, s_b);
printf(" diffuse: RGB: (%f,%f,%f)\n", c_r,c_g,c_b);
printf(" transparencyFactR: (1,%f)\n", tr_r);
printf(" transparencyFactG: (1,%f)\n", tr_g);
printf(" transparencyFactB: (1,%f)\n", tr_b);
printf(" specularFact: %f\n", s_fac);
printf(" faces: %d\n", faces);
printf(" texture: %d\n", t_id);
printf("}\n");
printf("Triangle: %d\n", s_idx+i*pa*2+2*j+1);
printf("{\n");
  printf(" a: Vertex: (%f,%f,%f)\n", fx11, fy11, fz11);
  printf(" b: Vertex: (%f,%f,%f)\n", fx22, fy22, fz22);
  printf(" c: Vertex: (%f,%f,%f)\n", fx12, fy12, fz12);
  printf(" texA: TexCoord: (%f,%f)\n",tx1,ty1);
  printf(" texB: TexCoord: (%f,%f)\n",tx2,ty2);
  printf(" texC: TexCoord: (%f,%f)\n",tx1,ty2);
printf(" na: Vector: (%f,%f,%f)\n", fx11, fy11, fz11);
printf(" nb: Vector: (%f,%f,%f)\n", fx22, fy22, fz22);
printf(" nc: Vector: (%f,%f,%f)\n", fx12, fy12, fz12);
printf(" transparency: RGB: (%f,%f,%f)\n", t_r, t_g, t_b);
printf(" specular: RGB: (%f,%f,%f)\n", s_r, s_g, s_b);
printf(" diffuse: RGB: (%f,%f,%f)\n", c_r,c_g,c_b);
printf(" transparencyFactR: (1,%f)\n", tr_r);
printf(" transparencyFactG: (1,%f)\n", tr_g);
printf(" transparencyFactB: (1,%f)\n", tr_b);
printf(" specularFact: %f\n", s_fac);
printf(" faces: %d\n", faces);
printf(" texture: %d\n", t_id);
printf("}\n");
   }
 }
}
void help()
{
 printf("params: t_idx (start idx)\n");
 printf("trans(rgb), spec(rgb), color(rgb), texid, transFact(rgb), specFact\n");
 printf("scale(xyz), rotate(xyz), translate(xyz) invN, faces, lt, pa, pb\n");
 printf("lt = 0, only tr, lt = 1 only listtransf, lt=2 both\n");
}

int main(int lb, char** par)
{
 double t_r,t_g,t_b,s_r,s_g,s_b,c_r,c_g,c_b,tr_r,tr_g,tr_b,s_fac;
 double x_x,y_y,z_z,r_x,r_y,r_z,d_x,d_y,d_z;
 int s_idx, t_id, inv_n, faces, lt, pa, pb;
 if (lb != 30) { printf("lbpar = %d\n", lb); help(); exit(1); }
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
  pa = atoi(par[28]);
  pb = atoi(par[29]);
 ball_gen(s_idx, 
 t_r, t_g, t_b, s_r, s_g, s_b, c_r, c_g, c_b, t_id, tr_r, tr_g, tr_b, s_fac,
 x_x, y_y, z_z, r_x, r_y, r_z, d_x, d_y, d_z, inv_n, faces, lt, pa, pb);
 return 0;
}
