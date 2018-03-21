#include "3dsLoader.h"
//#include "ImageIO.h"
//#ifdef UNIX
//#include <GL/glut.h>
//#else
//#include "GL/glut.h"
//#endif

//extern Loader3D loader;
//double angX, angY, angZ;
//double tX, tY, tZ;
//double water_level;
//double scale;
//bool want_water;

bool Loader3D::Load3DS(char* filename)
{
 //printf("Load3DS\n");
 strncpy(inputFileName, filename, 254);
 if (!mLoad3ds.Import3DS(&m3DModel, filename))
   {
    printf("Import failed.\n");
	return false;
   }
 //printf("Import Succeded.\n");
 printf("Number of materials: %d\n", m3DModel.numOfMaterials);
 printf("Number of objects: %d\n", m3DModel.numOfObjects);
 //GLUT_PreInit();
 /*Texture_3ds(&def_tex, "default.jpeg", 0);
 for(int i = 0; i < m3DModel.numOfMaterials; i++)				
	{
	 //printf("Processing material: %d\n", i);
	 //printf("Color ref is: %x %x %x\n", m3DModel.pMaterials[i].color[0],m3DModel.pMaterials[i].color[1], m3DModel.pMaterials[i].color[0]);
	 //printf("About getting texture: %s(%d)\n", m3DModel.pMaterials[i].strName, m3DModel.pMaterials[i].texureId);
	 if(strlen(m3DModel.pMaterials[i].strFile) > 0)				
		{
		 printf("Getting texture from file: %s\n", m3DModel.pMaterials[i].strFile);
		 char downcas[1024];
		 strcpy(downcas, m3DModel.pMaterials[i].strFile);
		 int len = strlen(downcas);
		 for (int ii=0;ii<len;ii++) if (downcas[ii] >= 'A' && downcas[ii] <= 'Z') downcas[ii] += 0x20;
		 if (!Texture_3ds(TextureArray3ds, downcas, i))
		    {
				printf("Cannot load texture: %s\n", downcas);
				return false;
		    }
		}
		m3DModel.pMaterials[i].texureId = i;						
	}*/
 return true;
}

/*void GLUT_DrawWater(double x, double y, double water_h)
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
}*/

/*void GLUT_display()
{
 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 glLoadIdentity();
 gluLookAt(0.0, 0.0, -20.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
 glTranslated(tX, tY, tZ);
 glRotated(angX, 1., 0., 0.);
 glRotated(angY, 0., 1., 0.);
 glRotated(angZ, 0., 0., 1.);
 glScaled(scale, scale, scale);
 loader.CallbackRender3DS();
 if (want_water) GLUT_DrawWater(1000., 1000., water_level);
 glFlush();
 glutSwapBuffers();
 glutPostRedisplay();
}*/

/*void GLUT_reshape(int w, int h)
{
 glViewport(0, 0, (GLsizei) w, (GLsizei) h);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 glFrustum(-1.0, 1.0, -1.0, 1.0, 1.5, 10000.0);
 glMatrixMode(GL_MODELVIEW);
}*/

/*void GLUT_ResetPosition()
{
 angX = angY = angZ = 0.;
 tX = tY = tZ = 0.;
 water_level = -0.4;
 scale = 1.;
}*/

/*void GLUT_keyboard(unsigned char key, int x, int y)
{
 switch (key)
   {
    case 'q': exit(1); break;
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
	case ' ': GLUT_ResetPosition();
	case '1': loader.WriteTris();
   }
}*/

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

/*void Loader3D::CallbackRender3DS()
{
 for(int i = 0; i < m3DModel.numOfObjects; i++)
   {		
	 if(m3DModel.pObject.size() <= 0) break;							
	 Object3D *pObject = &m3DModel.pObject[i];					
	 if(pObject->bHasTexture)									
		{	
		 glEnable(GL_TEXTURE_2D);									
		 glColor3ub(255, 255, 255);									
		 glBindTexture(GL_TEXTURE_2D, TextureArray3ds[pObject->materialID]); 
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
}*/

/*void GLUT_init()
{
 G1.12float ambientLight[] = {0.05f, 0.1f, 0.1f, 1.0f};
 G1.12float diffuseLight[] = {0.9f, 0.3f, 0.1f, 1.0f};
 G1.12float lightPos[] = {0.0f, 0.0f, 1.0f, 0.0f};

 G1.12float ambientLight2[] = {0.1f, 0.1f, 0.05f, 1.0f};
 G1.12float diffuseLight2[] = {0.1f, 0.3f, 0.9f, 1.0f};
 G1.12float lightPos2[] = {0.0f, 0.0f, -1.0f, 0.0f};

 G1.12float ambientLight3[] = {0.1f, 0.05f, 0.1f, 1.0f};
 G1.12float diffuseLight3[] = {0.1f, 0.9f, 0.3f, 1.0f};
 G1.12float lightPos3[] = {1.0f, 0.0f, 0.0f, 0.0f};

 G1.12float ambientLight4[] = {0.05f, 0.1f, 0.05f, 1.0f};
 G1.12float diffuseLight4[] = {0.9f, 0.1f, 0.9f, 1.0f};
 G1.12float lightPos4[] = {0.0f, -1.0f, 0.0f, 0.0f};

 angX = 30.;
 angY = 0.;
 angZ = 0;
 tX = tY = tZ = 0.;
 water_level = -0.4;
 scale = 1.;
 want_water = false;
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
}*/

/*void Loader3D::GLUT_PreInit()
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
}*/

/*void Loader3D::GLUT_Render3DS()
{
 glutMainLoop();
}*/

/*bool Loader3D::Texture_3ds(unsigned int textureArray[], char* strFileName, int ID)
{
 //FILE *pFile = NULL;									
 unsigned char *pJpeg = NULL;
 unsigned long width, height; 
 int type;
 JpegLoader jpeg;
 char tempstring[5] = {0};	
 strncpy(tempstring, strFileName + strlen(strFileName)-4, 4);	
 char FilePath[255];
 sprintf(FilePath, "texture/%s", strFileName);
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
}*/

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

