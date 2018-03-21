#include "3dsLoader.h"
#include "ImageIO.h"
#ifdef UNIX
#include <GL/glut.h>
#else
#include "GL/glut.h"
#endif

extern Loader3D loader;
double angX, angY, angZ;
double tX, tY, tZ;
double water_level;
double scale;
bool want_water;
bool want_normals;

bool Loader3D::Load3DS(char* filename, int lb, char** par)
{
 char syscmd[511];
 int tidd;
 //printf("Load3DS\n");
 tidd = 0;
 dargc = lb;
 dargv = par;
 ParseDatOptions();
 strncpy(inputFileName, filename, 254);
 if (!mLoad3ds.Import3DS(&m3DModel, filename))
   {
    printf("Import failed.\n");
	return false;
   }
 //printf("Import Succeded.\n");
 printf("Number of materials: %d\n", m3DModel.numOfMaterials);
 dinfo.datTexturesNum = m3DModel.numOfMaterials;
 printf("Number of objects: %d\n", m3DModel.numOfObjects);
 GLUT_PreInit();
 Texture_3ds(&def_tex, "default.jpeg", 0);
 for(int i = 0; i < m3DModel.numOfMaterials; i++)				
	{
	 //printf("Processing material: %d\n", i);
	 //printf("Color ref is: %x %x %x\n", m3DModel.pMaterials[i].color[0],m3DModel.pMaterials[i].color[1], m3DModel.pMaterials[i].color[0]);
	 //printf("About getting texture: %s(%d)\n", m3DModel.pMaterials[i].strName, m3DModel.pMaterials[i].texureId);
	 TextureDatCopied[i] = 0;
	 if (strlen(m3DModel.pMaterials[i].strFile) > 0)				
		{
		 printf("Getting texture from file: %s\n", m3DModel.pMaterials[i].strFile);
		 char downcas[1024];
		 strcpy(downcas, m3DModel.pMaterials[i].strFile);
		 int len = strlen(downcas);
		 for (int ii=0;ii<len;ii++) if (downcas[ii] >= 'A' && downcas[ii] <= 'Z') downcas[ii] += 0x20;
		 sprintf(syscmd, "cp %s/%s %s/%d.jpeg", dinfo.datTextureDir, downcas, dinfo.datTextureDir, i+1);
		 //printf("SCMD: %s\n", syscmd);
		 system(syscmd);
		 tidd = i;
		 if (!Texture_3ds(TextureArray3ds, downcas, i))
		    {
				printf("Cannot load texture: %s\n", downcas);
				tidd = 0;
				printf("\n\n\t\t--- NO (%d) TEXTURE DEFINITION ---\n\n",i);
		//		return false;
		    }
		 else tidd = i;
		 TextureDatCopied[i] = 1;
		}
	 //else printf("SCMD: no texture for (%d)\n", i);
	 m3DModel.pMaterials[i].texureId = tidd;						
	}
 return true;
}

void GLUT_DrawWater(double x, double y, double water_h)
{
 glDisable(GL_TEXTURE_2D);
 glBegin(GL_TRIANGLES);
 glColor4f(0.0f, 0.f, 1.0f, 0.33f);
 glNormal3d(0., -1., 0.);
 glVertex3d(-x, water_h, -y);
 glVertex3d(-x, water_h, y);
 glVertex3d(x, water_h, -y);

 glNormal3d(0., -1., 0.);
 glVertex3d(x, water_h, -y);
 glVertex3d(x, water_h, y);
 glVertex3d(-x, water_h, y);
 glEnd();
 glEnable(GL_TEXTURE_2D);
}

void GLUT_display()
{
 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 glLoadIdentity();
 gluLookAt(0.0, 0.0, -20.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
 glTranslated(tX, tY, tZ);
 glRotated(angX, 1., 0., 0.);
 glRotated(angY, 0., 1., 0.);
 glRotated(angZ, 0., 0., 1.);
 glScaled(scale, scale, scale);
 loader.SetDatTransformation();
 loader.CallbackRender3DS();
 if (want_water) GLUT_DrawWater(1000., 1000., water_level);
 glFlush();
 glutSwapBuffers();
 glutPostRedisplay();
}

void GLUT_reshape(int w, int h)
{
 glViewport(0, 0, (GLsizei) w, (GLsizei) h);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 10000.0);
 glMatrixMode(GL_MODELVIEW);
}

void GLUT_ResetPosition()
{
 angX = angY = angZ = 0.;
 tX = tY = tZ = 0.;
 water_level = -0.4;
 scale = 1.;
}

void GLUT_keyboard(unsigned char key, int x, int y)
{
 switch (key)
   {
    case 'q': exit(1); break;
	case 'n': want_normals = !want_normals; break;
	case 'w': angX += 10.; break;
	case 's': angX -= 10.; break;
	case 'd': angY += 10.; break;
	case 'a': angY -= 10.; break;
	case 'e': angZ += 10.; break;
	case 'x': angZ -= 10.; break;
	case 'i': tY += 1.5; break;
	case 'k': tY -= 1.5; break;
	case 'j': tX += 1.5; break;
	case 'l': tX -= 1.5; break;
	case 'o': tZ += 1.5; break;
	case 'm': tZ -= 1.5; break;
	case 'W': angX += 45.; break;
	case 'S': angX -= 45.; break;
	case 'D': angY += 45.; break;
	case 'A': angY -= 45.; break;
	case 'E': angZ += 45.; break;
	case 'X': angZ -= 45.; break;
	case 'I': tY += 100; break;
	case 'K': tY -= 100; break;
	case 'J': tX += 100; break;
	case 'L': tX -= 100; break;
	case 'O': tZ += 100; break;
	case 'M': tZ -= 100; break;
	case '=': water_level += .3; break;
	case '-': water_level -= .3; break;
	case 'y': scale /= 2.5; break;
	case 'u': scale *= 2.5; break;
	case 'Y': scale /= 2.5; break;
	case 'U': scale *= 2.5; break;
	case '\\': 
		want_water = ! want_water; 
		if (!want_water) glDisable(GL_BLEND);
		else 
		  {
		   glEnable(GL_BLEND);
		   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		  }
		break;
	case ' ': GLUT_ResetPosition(); break;
	case '1': loader.WriteTris(false); break;
	case '2': loader.WriteTris(true); break;
	case '3': loader.WriteDat(); break;
   }
}

void Loader3D::SetDatTransformation()
{
 dinfo.datTranslationX = tX;
 dinfo.datTranslationY = tY;
 dinfo.datTranslationZ = tZ;
 dinfo.datRotationX = angX;
 dinfo.datRotationY = angY;
 dinfo.datRotationZ = angZ;
 dinfo.datScale = scale;

}

void Loader3D::SetDatDefaults()
{
 dinfo.datScreenX = 1024;
 dinfo.datScreenY = 768;
 dinfo.datBackup = 1024;
 dinfo.datObserverX = 0.;
 dinfo.datObserverY = 0.;
 dinfo.datObserverZ = -400.;
 dinfo.datLightX = 50.;
 dinfo.datLightY = 50.;
 dinfo.datLightZ = -200.;
 dinfo.datTranslationX = 0.;
 dinfo.datTranslationY = 0.;
 dinfo.datTranslationZ = 0.;
 dinfo.datRotationX = 0.;
 dinfo.datRotationY = 0.;
 dinfo.datRotationZ = 0.;
 dinfo.datScale= 1.;
 dinfo.datNoTex = 0;
 strcpy(dinfo.datTextureDir, "texture");
 dinfo.datColorRed = 3.;
 dinfo.datColorGreen = 3.;
 dinfo.datColorBlue = 3.;
 dinfo.datSpecularRed = 3;
 dinfo.datSpecularGreen = 2;
 dinfo.datSpecularBlue = 1;
 dinfo.datTransRed = 1;
 dinfo.datTransGreen = 2;
 dinfo.datTransBlue = 3;
 dinfo.datTransFactorRed = 2.4;
 dinfo.datTransFactorGreen = 2.6;
 dinfo.datTransFactorBlue = 2.9;
 dinfo.datFaces = 2;
 dinfo.datSpecularFactor = 32.;
}

void Loader3D::ParseDatOptions()
{
 char opt[512], value[512];
 SetDatDefaults();
 for (int i=1;i<dargc;i++)
   {
    if (!dargv[i]) continue;
    if (sscanf(dargv[i], "%s = %s", opt, value) != 2) continue;
    if (!strcmp(opt, "cr")) sscanf(value, "%lf",&dinfo.datColorRed);
    if (!strcmp(opt, "cg")) sscanf(value, "%lf",&dinfo.datColorGreen);
    if (!strcmp(opt, "cb")) sscanf(value, "%lf",&dinfo.datColorBlue);
    if (!strcmp(opt, "tfr")) sscanf(value, "%lf",&dinfo.datTransFactorRed);
    if (!strcmp(opt, "tfg")) sscanf(value, "%lf",&dinfo.datTransFactorGreen);
    if (!strcmp(opt, "tfb")) sscanf(value, "%lf",&dinfo.datTransFactorBlue);
    if (!strcmp(opt, "tr")) sscanf(value, "%lf",&dinfo.datTransRed);
    if (!strcmp(opt, "tg")) sscanf(value, "%lf",&dinfo.datTransGreen);
    if (!strcmp(opt, "tb")) sscanf(value, "%lf",&dinfo.datTransBlue);
    if (!strcmp(opt, "sr")) sscanf(value, "%lf",&dinfo.datSpecularRed);
    if (!strcmp(opt, "sg")) sscanf(value, "%lf",&dinfo.datSpecularGreen);
    if (!strcmp(opt, "sb")) sscanf(value, "%lf",&dinfo.datSpecularBlue);
    if (!strcmp(opt, "sf")) sscanf(value, "%lf",&dinfo.datSpecularFactor);
    if (!strcmp(opt, "faces")) sscanf(value, "%d",&dinfo.datFaces);
    if (!strcmp(opt, "notex")) sscanf(value, "%d",&dinfo.datNoTex);
    if (!strcmp(opt, "tx")) sscanf(value, "%lf",&dinfo.datTranslationX);
    if (!strcmp(opt, "ty")) sscanf(value, "%lf",&dinfo.datTranslationY);
    if (!strcmp(opt, "tz")) sscanf(value, "%lf",&dinfo.datTranslationZ);
    if (!strcmp(opt, "rx")) sscanf(value, "%lf",&dinfo.datRotationX);
    if (!strcmp(opt, "ry")) sscanf(value, "%lf",&dinfo.datRotationY);
    if (!strcmp(opt, "rz")) sscanf(value, "%lf",&dinfo.datRotationZ);
    if (!strcmp(opt, "scale")) sscanf(value, "%lf",&dinfo.datScale);
    printf("setting '%s' to '%s'\n", opt, value);
   }
 angX = dinfo.datRotationX;
 angY = dinfo.datRotationY;
 angZ = dinfo.datRotationZ;
 tX = dinfo.datTranslationX;
 tY = dinfo.datTranslationY;
 tZ = dinfo.datTranslationZ;
 scale = dinfo.datScale;
}

void Loader3D::WriteDatHeader(FILE* f)
{
 fprintf(f, "Screen: (%d,%d)\n", dinfo.datScreenX, dinfo.datScreenY);	
 fprintf(f, "Backup: %d\n", dinfo.datBackup);
 fprintf(f, "Observer: Vertex: (%f,%f,%f)\n", dinfo.datObserverX, dinfo.datObserverY, dinfo.datObserverZ);
 fprintf(f, "Light: Vertex: (%f,%f,%f)\n", dinfo.datLightX, dinfo.datLightY, dinfo.datLightZ);
 fprintf(f, "TexDirectory: %s\n", dinfo.datTextureDir);
 fprintf(f, "NumTextures: %d\n", dinfo.datTexturesNum);
 fprintf(f, "nTriangles: %d\n", dinfo.datTrisNum);
/* fprintf(f, "nNURBS: 0\n", 0);*/
 fprintf(f, "ListTransform: [0,%d]\n", dinfo.datTrisNum-1);
 fprintf(f, "{\n");
 fprintf(f, " Translate: (%lf,%lf,%lf)\n", dinfo.datTranslationX, dinfo.datTranslationY, dinfo.datTranslationZ);
 fprintf(f, " Rotate: (%lf,%lf,%lf)\n", dinfo.datRotationX, dinfo.datRotationY, -dinfo.datRotationZ);
 fprintf(f, " Scale: (%lf,%lf,%lf)\n", dinfo.datScale*3.6, dinfo.datScale*3.6, dinfo.datScale*3.6);
 fprintf(f, "}\n");
 
}

void Loader3D::TextureMapping(double& tc)
{
 if (tc < 1e-4) tc = 0.;
 if (tc > 1.-1e-4) tc = 1.;
}

void Loader3D::WriteDat()
{
 //int texid;
 FILE* plik;
 char name[512];
 strcpy(name, inputFileName);
 strcat(name, ".dat");
 plik = fopen(name, "w");
 if (!plik) return;
 int* idxt = new int[m3DModel.numOfObjects];
 int k;
 k = 0;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
	{
	 idxt[i] = k;
	 k += m3DModel.pObject[i].numOfFaces;
	 //printf("i = %d, idxt[i] = %d, k = %d, nV = %d\n", i, idxt[i], k, m3DModel.pObject[i].numOfVerts);
	}
 dinfo.datTrisNum = k;
 WriteDatHeader(plik);
 
 int i1,i2,i3;
 int tid;
 double x1,x2,x3;
 double nx1,nx2,nx3;
 double y1,y2,y3;
 double ny1,ny2,ny3;
 double z1,z2,z3;
 double nz1,nz2,nz3;
 double u1,u2,u3;
 double v1,v2,v3;

 k = 0;

 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {
    Object3D *pObject = &m3DModel.pObject[i];
    for(int j = 0; j < pObject->numOfFaces; j++)
	   {		
	    if (pObject->bHasTexture)
	      {
	       tid = pObject->materialID+1;
	       if (!TextureDatCopied[pObject->materialID]) tid = 0;
	//	  printf("matID: %d, TID=%d\n", pObject->materialID, tid);
	      }
	   else tid = 0;
	   
	   i1 = pObject->pFaces[j].vertIndex[0];
	   i2 = pObject->pFaces[j].vertIndex[1];
	   i3 = pObject->pFaces[j].vertIndex[2];
	   
	   x1 = pObject->pVerts[i1].x;
	   x2 = pObject->pVerts[i2].x;
	   x3 = pObject->pVerts[i3].x;
	   
	   y1 = pObject->pVerts[i1].y;
	   y2 = pObject->pVerts[i2].y;
	   y3 = pObject->pVerts[i3].y;
	   
	   z1 = pObject->pVerts[i1].z;
	   z2 = pObject->pVerts[i2].z;
	   z3 = pObject->pVerts[i3].z;

	   if (want_normals)
	     {
	      nx1 = pObject->pNormals[i1].x;
	      nx2 = pObject->pNormals[i2].x;
	      nx3 = pObject->pNormals[i3].x;
	   
	      ny1 = pObject->pNormals[i1].y;
	      ny2 = pObject->pNormals[i2].y;
	      ny3 = pObject->pNormals[i3].y;
	   
	      nz1 = pObject->pNormals[i1].z;
	      nz2 = pObject->pNormals[i2].z;
	      nz3 = pObject->pNormals[i3].z;
	     }
	   else
	     {
	      nx1 = nx2 = nx3 = 0.;
	      ny1 = ny2 = ny3 = 0.;
	      nz1 = nz2 = nz3 = 0.;
	     }

	   if (pObject->pTexVerts)
	     {
	      /*u1 = pObject->pTexVerts[i1].x - floor(pObject->pTexVerts[i1].x);
	      u2 = pObject->pTexVerts[i2].x - floor(pObject->pTexVerts[i2].x);
	      u3 = pObject->pTexVerts[i3].x - floor(pObject->pTexVerts[i3].x);
	   
	      v1 = pObject->pTexVerts[i1].y - floor(pObject->pTexVerts[i1].y);
	      v2 = pObject->pTexVerts[i2].y - floor(pObject->pTexVerts[i2].y);
	      v3 = pObject->pTexVerts[i3].y - floor(pObject->pTexVerts[i3].y);*/
	      u1 = pObject->pTexVerts[i1].x;
	      u2 = pObject->pTexVerts[i2].x;
	      u3 = pObject->pTexVerts[i3].x;
	   
	      v1 = pObject->pTexVerts[i1].y;
	      v2 = pObject->pTexVerts[i2].y;
	      v3 = pObject->pTexVerts[i3].y;

	     }
	   else
	     { 
              tid = 0;
	      u1 = u2 = u3 = v1 = v2 = v3 = 0.;
	     }
	   /*if (tid == 7 && j < 2)
	     {
              printf("%d: %f %f, %f %f, %f %f\n", j, u1, v1, u2, v2, u3, v3);
	     }*/
	   TextureMapping(u1);
	   TextureMapping(u2);
	   TextureMapping(u3);
	   TextureMapping(v1);
	   TextureMapping(v2);
	   TextureMapping(v3);

	   u1 = 1. - u1;
	   u2 = 1. - u2;
	   u3 = 1. - u3;
	   //v1 = 1. - v1;
	   //v2 = 1. - v2;
	   //v3 = 1. - v3;
           fprintf(plik, "Triangle: %d\n", k);
           fprintf(plik, "{\n");
           fprintf(plik, "a: Vertex: (%lf,%lf,%lf)\n",x1,y1,z1);
           fprintf(plik, "b: Vertex: (%lf,%lf,%lf)\n",x2,y2,z2);
           fprintf(plik, "c: Vertex: (%lf,%lf,%lf)\n",x3,y3,z3);
           fprintf(plik, "texA: TexCoord: (%lf,%lf)\n", u1, v1);
           fprintf(plik, "texB: TexCoord: (%lf,%lf)\n", u2, v2);
           fprintf(plik, "texC: TexCoord: (%lf,%lf)\n", u3, v3);
           fprintf(plik, "na: Vector: (%lf,%lf,%lf)\n", nx1, ny1, nz1);
           fprintf(plik, "nb: Vector: (%lf,%lf,%lf)\n", nx2, ny2, nz2);
           fprintf(plik, "nc: Vector: (%lf,%lf,%lf)\n", nx3, ny3, nz3);
           /*fprintf(plik, "na: Vector: (%lf,%lf,%lf)\n", .0, .0, .0);
           fprintf(plik, "nb: Vector: (%lf,%lf,%lf)\n", .0, .0, .0);
           fprintf(plik, "nc: Vector: (%lf,%lf,%lf)\n", .0, .0, .0);*/
           fprintf(plik, "transparency: RGB: (%lf,%lf,%lf)\n", dinfo.datTransRed, dinfo.datTransGreen, dinfo.datTransBlue);
           fprintf(plik, "specular: RGB: (%lf,%lf,%lf)\n", dinfo.datSpecularRed, dinfo.datSpecularGreen, dinfo.datSpecularBlue);
           if (pObject->bHasTexture) 
             fprintf(plik, "diffuse: RGB: (%lf,%lf,%lf)\n", dinfo.datColorRed, dinfo.datColorGreen, dinfo.datColorBlue);
	    else 
              {			
	       if (m3DModel.pMaterials.size() < (unsigned int)pObject->materialID) 
	          {
		   unsigned char *pColor = m3DModel.pMaterials[pObject->materialID].color;
                   fprintf(plik, "diffuse: RGB: (%lf,%lf,%lf)\n", (double)pColor[0]/255., (double)pColor[1]/255., (double)pColor[2]/255.);
		  }
	       else 
                  fprintf(plik, "diffuse: RGB: (%lf,%lf,%lf)\n", dinfo.datColorRed, dinfo.datColorGreen, dinfo.datColorBlue);
	      }		
           fprintf(plik, "surface: ABCD: (0,0,0,0)\n");
           fprintf(plik, "transparencyFactR: (1,%lf)\n", dinfo.datTransFactorRed);
           fprintf(plik, "transparencyFactG: (1,%lf)\n", dinfo.datTransFactorGreen);
           fprintf(plik, "transparencyFactB: (1,%lf)\n", dinfo.datTransFactorBlue);
           fprintf(plik, "specularFact: %lf\n", dinfo.datSpecularFactor);
           fprintf(plik, "faces: %d\n",dinfo.datFaces);
	   if (dinfo.datNoTex) fprintf(plik, "texture: %d\n", 0);
	   else fprintf(plik, "texture: %d\n", tid);
           fprintf(plik, "}\n");
	   k ++;
	  }		
   }
/* return;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {		
	 if(m3DModel.pObject.size() <= 0) break;							
	 Object3D *pObject = &m3DModel.pObject[i];					
	 for (int j=0;j<pObject->numOfVerts;j++)
		{
		 if (pObject->pTexVerts)
			{
                         if (!binary)				
			        fprintf(plik, "%1.12f %1.12f %1.12f %1.12f %1.12f %1.12f %1.12f %1.12f\n",
				pObject->pVerts[j].x, pObject->pVerts[j].y, pObject->pVerts[j].z,
				pObject->pNormals[j].x, pObject->pNormals[j].y, pObject->pNormals[j].z, 
				pObject->pTexVerts[j].x - floor(pObject->pTexVerts[j].x),
				pObject->pTexVerts[j].y - floor(pObject->pTexVerts[j].y));
			 else
			   {
                         tmp = pObject->pVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pTexVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pTexVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			   }
			}
		 else
			{
			 if (!binary)
			        fprintf(plik, "%1.12f %1.12f %1.12f %1.12f %1.12f %1.12f 0.00001 0.00001\n",
				pObject->pVerts[j].x, pObject->pVerts[j].y, pObject->pVerts[j].z,
				pObject->pNormals[j].x, pObject->pNormals[j].y, pObject->pNormals[j].z);
                         else
			   {
                         tmp = pObject->pVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = 0.;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = 0.;
			 fwrite(&tmp, sizeof(float), 1, plik);
			   }
			}
		}
   }
 k = 0;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {
	   k += m3DModel.pObject[i].numOfFaces;
   }
 if (!binary) fprintf(plik, "%d\n", k);
 else         fwrite(&k, sizeof(int), 1, plik);
 int itmp;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {
    k = idxt[i];
	Object3D *pObject = &m3DModel.pObject[i];
    for(int j = 0; j < pObject->numOfFaces; j++)
	   {		
	    if (!binary) fprintf(plik, "%d %d %d\n", 
			pObject->pFaces[j].vertIndex[0]+k+1,
			pObject->pFaces[j].vertIndex[1]+k+1,
			pObject->pFaces[j].vertIndex[2]+k+1);
	    else
	      {
	   itmp = pObject->pFaces[j].vertIndex[0]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	   itmp = pObject->pFaces[j].vertIndex[1]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	   itmp = pObject->pFaces[j].vertIndex[2]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	      }
	  }		
   }*/
 fclose(plik);
 printf("File: %s written.\n", name);
}

void Loader3D::WriteTris(bool binary)
{
 //int texid;
 FILE* plik;
 char name[512];
 strcpy(name, inputFileName);
 strcat(name, ".tri");
 if (binary)  plik = fopen(name, "wb");
 else         plik = fopen(name, "w");
 if (!plik) return;
 int* idxt = new int[m3DModel.numOfObjects];
 int k;
 k = 0;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
	{
	 idxt[i] = k;
	 k += m3DModel.pObject[i].numOfVerts;
	 //printf("i = %d, idxt[i] = %d, k = %d, nV = %d\n", i, idxt[i], k, m3DModel.pObject[i].numOfVerts);
	}
float tmp;
if (!binary) fprintf(plik, "%d\n", k);
else         fwrite(&k, sizeof(int), 1, plik);
for(int i = 0; i < m3DModel.numOfObjects; i++)
   {		
	 if(m3DModel.pObject.size() <= 0) break;							
	 Object3D *pObject = &m3DModel.pObject[i];					
	 for (int j=0;j<pObject->numOfVerts;j++)
		{
		 if (pObject->pTexVerts)
			{
                         if (!binary)				
			        fprintf(plik, "%1.12f %1.12f %1.12f %1.12f %1.12f %1.12f %1.12f %1.12f\n",
				pObject->pVerts[j].x, pObject->pVerts[j].y, pObject->pVerts[j].z,
				pObject->pNormals[j].x, pObject->pNormals[j].y, pObject->pNormals[j].z, 
				pObject->pTexVerts[j].x - floor(pObject->pTexVerts[j].x),
				pObject->pTexVerts[j].y - floor(pObject->pTexVerts[j].y));
			 else
			   {
                         tmp = pObject->pVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pTexVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pTexVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			   }
			}
		 else
			{
			 if (!binary)
			        fprintf(plik, "%1.12f %1.12f %1.12f %1.12f %1.12f %1.12f 0.00001 0.00001\n",
				pObject->pVerts[j].x, pObject->pVerts[j].y, pObject->pVerts[j].z,
				pObject->pNormals[j].x, pObject->pNormals[j].y, pObject->pNormals[j].z);
                         else
			   {
                         tmp = pObject->pVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = 0.;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = 0.;
			 fwrite(&tmp, sizeof(float), 1, plik);
			   }
			}
		}
   }
 k = 0;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {
	   k += m3DModel.pObject[i].numOfFaces;
   }
 if (!binary) fprintf(plik, "%d\n", k);
 else         fwrite(&k, sizeof(int), 1, plik);
 int itmp;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {
    k = idxt[i];
	Object3D *pObject = &m3DModel.pObject[i];
    for(int j = 0; j < pObject->numOfFaces; j++)
	   {		
	    if (!binary) fprintf(plik, "%d %d %d\n", 
			pObject->pFaces[j].vertIndex[0]+k+1,
			pObject->pFaces[j].vertIndex[1]+k+1,
			pObject->pFaces[j].vertIndex[2]+k+1);
	    else
	      {
	   itmp = pObject->pFaces[j].vertIndex[0]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	   itmp = pObject->pFaces[j].vertIndex[1]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	   itmp = pObject->pFaces[j].vertIndex[2]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	      }
	  }		
   }
 fclose(plik);
}

/*void Loader3D::WriteTrisBinary()
{
 //int texid;
 FILE* plik;
 char name[512];
 strcpy(name, inputFileName);
 strcat(name, ".tri");
 plik = fopen(name, "wb");
 if (!plik) return;
 int* idxt = new int[m3DModel.numOfObjects];
 int k;
 k = 0;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
	{
	 idxt[i] = k;
	 k += m3DModel.pObject[i].numOfVerts;
	 printf("i = %d, idxt[i] = %d, k = %d, nV = %d\n", i, idxt[i], k, m3DModel.pObject[i].numOfVerts);
	}
float tmp;
// fprintf(plik, "%d\n", k);
fwrite(&k, sizeof(int), 1, plik);
for(int i = 0; i < m3DModel.numOfObjects; i++)
   {		
	 if(m3DModel.pObject.size() <= 0) break;							
	 Object3D *pObject = &m3DModel.pObject[i];					
	 for (int j=0;j<pObject->numOfVerts;j++)
		{
		 if (pObject->pTexVerts)
			{
             tmp = pObject->pVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pTexVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pTexVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			}
		 else
			{
             tmp = pObject->pVerts[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pVerts[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].x;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].y;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = pObject->pNormals[j].z;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = 0.;
			 fwrite(&tmp, sizeof(float), 1, plik);
			 tmp = 0.;
			 fwrite(&tmp, sizeof(float), 1, plik);
			}
		}
   }
 k = 0;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {
	   k += m3DModel.pObject[i].numOfFaces;
   }
 //fprintf(plik, "%d\n", k);
 fwrite(&k, sizeof(int), 1, plik);
 int itmp;
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {
    k = idxt[i];
	Object3D *pObject = &m3DModel.pObject[i];
    for(int j = 0; j < pObject->numOfFaces; j++)
	   {		
	   itmp = pObject->pFaces[j].vertIndex[0]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	   itmp = pObject->pFaces[j].vertIndex[1]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	   itmp = pObject->pFaces[j].vertIndex[2]+k+1;
	   fwrite(&itmp, sizeof(int), 1, plik);
	  }		
   }
 fclose(plik);
}*/

void Loader3D::CallbackRender3DS()
{
/* int tid;*/
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {		
	 if(m3DModel.pObject.size() <= 0) break;							
	 Object3D *pObject = &m3DModel.pObject[i];					
/*	 tid = -1;*/
	 if(pObject->bHasTexture)									
		{	
		 glEnable(GL_TEXTURE_2D);									
		 glColor3ub(255, 255, 255);									
		 glBindTexture(GL_TEXTURE_2D, TextureArray3ds[pObject->materialID]); 
/*		 tid = TextureArray3ds[pObject->materialID]; */
		} 
	 else if (pObject->bHasTexture && pObject->pTexVerts)
		{
		 glEnable(GL_TEXTURE_2D);									
		 glColor3ub(255, 255, 255);								
		 glBindTexture(GL_TEXTURE_2D, def_tex); 
		}
	 else
	    {
		 glDisable(GL_TEXTURE_2D);									
		 glColor3ub(235, 235, 235);								
	    }
    glBegin(GL_TRIANGLES);												
    for(int j = 0; j < pObject->numOfFaces; j++)
	   {		
	    for (int whichVertex = 0; whichVertex < 3; whichVertex++)
		   {
		    int index = pObject->pFaces[j].vertIndex[whichVertex];		
		    glNormal3f(pObject->pNormals[ index ].x, pObject->pNormals[ index ].y, pObject->pNormals[ index ].z);				
		    if (pObject->bHasTexture) 
		      {			
			   if(pObject->pTexVerts) 
			     {
			      glTexCoord2f(pObject->pTexVerts[ index ].x, pObject->pTexVerts[ index ].y);
/*			      if (tid == 6 && j < 2) printf("%d: %f,%f\n", j, pObject->pTexVerts[ index ].x, pObject->pTexVerts[ index ].y);*/
			     }
		      } 
		   else 
		     {			
			  if (m3DModel.pMaterials.size() < (unsigned int)pObject->materialID) 
			    {
			     unsigned char *pColor = m3DModel.pMaterials[pObject->materialID].color;
			     glColor3ub(pColor[0], pColor[1], pColor[2]);
			    }
		     }		
		   glVertex3f(pObject->pVerts[ index ].x, pObject->pVerts[ index ].y, pObject->pVerts[ index ].z);
		  }
	  }		
    glEnd();
   }
}

void GLUT_init()
{
 GLfloat ambientLight[] = {0.05f, 0.1f, 0.1f, 1.0f};
 GLfloat diffuseLight[] = {0.9f, 0.3f, 0.1f, 1.0f};
 GLfloat lightPos[] = {0.0f, 0.0f, 1.0f, 0.0f};

 GLfloat ambientLight2[] = {0.1f, 0.1f, 0.05f, 1.0f};
 GLfloat diffuseLight2[] = {0.1f, 0.3f, 0.9f, 1.0f};
 GLfloat lightPos2[] = {0.0f, 0.0f, -1.0f, 0.0f};

 GLfloat ambientLight3[] = {0.1f, 0.05f, 0.1f, 1.0f};
 GLfloat diffuseLight3[] = {0.1f, 0.9f, 0.3f, 1.0f};
 GLfloat lightPos3[] = {1.0f, 0.0f, 0.0f, 0.0f};

 GLfloat ambientLight4[] = {0.05f, 0.1f, 0.05f, 1.0f};
 GLfloat diffuseLight4[] = {0.9f, 0.1f, 0.9f, 1.0f};
 GLfloat lightPos4[] = {0.0f, -1.0f, 0.0f, 0.0f};

 /*angX = 30.;
 angY = 0.;
 angZ = 0;
 scale = 1.;
 tX = tY = tZ = 0.;*/
 water_level = -0.4;
 want_water = false;
 want_normals = true;
 glShadeModel(GL_SMOOTH);
 //glEnable(GL_BLEND);
 glEnable(GL_DEPTH_TEST);
 glDepthFunc(GL_LEQUAL);
 glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
 glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight2);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPos2);
	glEnable(GL_LIGHT1);
	glLightfv(GL_LIGHT2, GL_AMBIENT, ambientLight3);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseLight3);
	glLightfv(GL_LIGHT2, GL_POSITION, lightPos3);
	glEnable(GL_LIGHT2);
	glLightfv(GL_LIGHT3, GL_AMBIENT, ambientLight4);
	glLightfv(GL_LIGHT3, GL_DIFFUSE, diffuseLight4);
	glLightfv(GL_LIGHT3, GL_POSITION, lightPos4);
	glEnable(GL_LIGHT3);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
 glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
}

void Loader3D::GLUT_PreInit()
{
 int argc;
 char** argv;
 argc = 1;
 argv = new char*[1];
 argv[0] = new char[strlen(inputFileName)+1];
 strcpy(argv[0], inputFileName);
 glutInit(&argc, argv);
 glutInitDisplayMode(GLUT_DOUBLE);
 glutInitWindowSize(600, 600);
 glutInitWindowPosition(10, 10);
 glutCreateWindow(argv[0]);
 glutDisplayFunc(GLUT_display);
 glutReshapeFunc(GLUT_reshape);
 glutKeyboardFunc(GLUT_keyboard);
 GLUT_init();
}

void Loader3D::GLUT_Render3DS()
{
 glutMainLoop();
}

bool Loader3D::Texture_3ds(unsigned int textureArray[], char* strFileName, int ID)
{
 //FILE *pFile = NULL;									
 unsigned char *pJpeg = NULL;
 unsigned long width, height; 
 int type;
 JpegLoader jpeg;
 char tempstring[5] = {0};	
 strncpy(tempstring, strFileName + strlen(strFileName)-4, 4);	
 char FilePath[255];
 sprintf(FilePath, "%s/%s", dinfo.datTextureDir, strFileName);
 if(!strFileName) return false;	
 if (!strcmp(tempstring, ".JPG") || !strcmp(tempstring, "JPEG") ||
		!strcmp(tempstring, ".jpg") || !strcmp(tempstring, "jpeg")) 
	     {
		  printf("Loading JPEG texture: %s\n", FilePath);
		  jpeg.GetSize(FilePath, width, height);
		  printf("Size is: (%dx%d)\n", (int)width, (int)height);
		  pJpeg = new unsigned char[width*height*3];
		  if (!jpeg.Load(FilePath, 24, pJpeg)) 
		    {
				printf("Cannot load JPEG texture: %s\n", FilePath);
				return false;
		    }
		  printf("JPEG texture succesfully loaded.\n");
		  type = 2;
	     }
 else if (!strcmp(tempstring, ".BMP") || !strcmp(tempstring, ".bmp"))
    {
	 printf("Loading BMP texture: %s\n", FilePath);
	 if (!LoadBMP(FilePath, pJpeg, width, height)) 
		    {
				printf("Cannot load BMP texture: %s\n", FilePath);
				return false;
		    }
	 printf("BMP texture succesfully loaded.\n");
	 type = 1;
   }
 else if (!strcmp(tempstring, ".TGA") || !strcmp(tempstring, ".tga"))
    {
	 printf("Loading TGA texture: %s\n", FilePath);
	 if (!LoadTGA(FilePath, pJpeg, width, height)) 
		    {
				printf("Cannot load TGA texture: %s\n", FilePath);
				return false;
		    }
	 printf("TGA texture succesfully loaded.\n");
	 type = 3;
   }
 else return false;
 glGenTextures(1, &textureArray[ID]);
 printf("type = %d\n", type);
 glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
 glBindTexture(GL_TEXTURE_2D, textureArray[ID]);
 gluBuild2DMipmaps(GL_TEXTURE_2D, 3, width, height, GL_RGB, GL_UNSIGNED_BYTE, pJpeg);
 //gluBuild2DMipmaps(GL_TEXTURE_2D, 3, width, height, GL_RGB, GL_UNSIGNED_BYTE, pJpeg);
 glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
 glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	
 if (pJpeg) delete [] pJpeg;
 return true;
}

void Loader3D::Unload3DS()
{
 //printf("Freeying memory.\n");
 //printf("Number of objects: %d\n", m3DModel.numOfObjects);
 for(int i = 0; i < m3DModel.numOfObjects; i++)	
	{
	 //printf("Freeying object #%d\n", i);
	 delete [] m3DModel.pObject[i].pFaces;
	 delete [] m3DModel.pObject[i].pNormals;
	 delete [] m3DModel.pObject[i].pVerts;
	 delete [] m3DModel.pObject[i].pTexVerts;
	}
}
