#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <GL/glut.h>
#define PI 3.1415926

enum enumTexture { PCX, BMP};

time_t t1,t2;
int fps;
int hdr;
char strModel[256];
char strTexture[256];
float angle = 0.0f;				
float radians = 0.0f;			
float percent;
float angx, angy, angz, scale;
int idx;
int fac, tid;
float sf, nd;
float fr,fg,fb;
float tr,tg,tb;
float sr,sg,sb;
float cr,cg,cb;

typedef struct _Texture
{
	int textureType;
	int width;					
	int height;					
	long int scaledWidth;
	long int scaledHeight;

	unsigned int texID;			
	unsigned char *data;		
	unsigned char *palette;
} Texture;


typedef struct _PcxHeader
{
	unsigned char manufacturer;
	unsigned char version;
	unsigned char encoding;
	unsigned char bits;
	unsigned char xMin;
	unsigned char yMin;
	unsigned char xMax;
	unsigned char yMax;
	unsigned char *palette;
} PcxHeader;

typedef struct _BMPTag
{
 char ident[2];
 int fsize;
 int dummy;
 int offset;
 int dummy2;
 int bm_x;
 int bm_y;
 short planes;
 short bpp;
 int compress;
 int nbytes;
 int no_matter[4];
} BMPTag;

typedef struct _Vector
{
   float point[3];
} Vector;


typedef struct _TexCoord
{
   float s;
   float t;
} TexCoord;


typedef struct _StIndex
{
   short s;
   short t;
} StIndex;


typedef struct _CustomVertex
{
   unsigned char v[3];
   unsigned char normalIndex;	
} CustomVertex;


typedef struct _Frame
{
   float scale[3];
   float translate[3];
   char name[16];
   CustomVertex fp[1];
} Frame;


typedef struct _Mesh
{
   unsigned short meshIndex[3];		
   unsigned short stIndex[3];		
} Mesh;


typedef struct _ModelData
{
	int numFrames;				
	int numPoints;				
	int numTriangles;			
    	int numST;					
	int frameSize;				
	int texWidth, texHeight;	
	int currentFrame;			
	int nextFrame;				
	float interpol;				
	Mesh *triIndex;			
	TexCoord *st;				
	Vector *pointList;		
	Texture *modelTex;		
} ModelData;

typedef struct _ModelHeader
{
   int ident;		
   int version;		
   int skinwidth;    
   int skinheight;   
   int framesize;    
   int numSkins;     
   int numXYZ;       
   int numST;        
   int numTris;      
   int numGLcmds;
   int numFrames;    
   int offsetSkins;  
   int offsetST;     
   int offsetTris;   
   int offsetFrames; 
   int offsetGLcmds; 
   int offsetEnd;    
} ModelHeader;

ModelData *myModel;
float interValue = 0.0f;			

void free_model(ModelData *model)
{
 if (model->triIndex != NULL)  free(model->triIndex);
 if (model->pointList != NULL) free(model->pointList);
 if (model->st != NULL)        free(model->st);
 if (model != NULL)            free(model);
}

void init_bmp(BMPTag* b)
{
 int i;
 b->ident[0]='B';
 b->ident[1]='M';
 b->fsize=0;
 b->dummy=0;
 b->offset=sizeof(BMPTag);
 b->bm_x=b->bm_y=0x20;
 b->dummy2=40;
 b->bpp=0x18;
 b->planes=1;
 b->compress=0;
 b->nbytes=3*32*32;
 for (i=0;i<4;i++) b->no_matter[i]=0;
}




void write_as_bmp(char* fn, unsigned char* pix, int x, int y)
{
 char temp[512];
 int i,j;
 FILE* plik;
 BMPTag bm_handle;
 strcpy(temp, fn);
 if (strstr(temp, ".")) 
   {
    i = 0;
    while (temp[i] != '.') i++;
    temp[i] = 0;
    strcat(temp, ".bmp");
   }
 else strcat(temp, ".bmp");
 plik = fopen(temp, "w");
 if (!plik) return;
 init_bmp(&bm_handle);
 fprintf(plik,"%c%c",'B', 'M');
 bm_handle.bm_y = y;
 bm_handle.bm_x = x;
 bm_handle.fsize = sizeof(BMPTag)+(bm_handle.bm_y*bm_handle.bm_x*3);
 fwrite(&bm_handle.fsize,4,1,plik);
 fwrite(&bm_handle.dummy,4,1,plik);
 bm_handle.offset=sizeof(BMPTag);
 bm_handle.planes=1;
 bm_handle.bpp=24;
 fwrite(&bm_handle.offset,4,1,plik);
 fwrite(&bm_handle.dummy2,4,1,plik);
 fwrite(&bm_handle.bm_x,4,1,plik);
 fwrite(&bm_handle.bm_y,4,1,plik);
 fwrite(&bm_handle.planes,2,1,plik);
 fwrite(&bm_handle.bpp,2,1,plik);
 fwrite(&bm_handle.compress,4,1,plik);
 fwrite(&bm_handle.nbytes,4,1,plik);
 for (i=0;i<4;i++)  fwrite(&bm_handle.no_matter[i],4,1,plik);
 fseek(plik,bm_handle.offset,SEEK_SET);
 for (i=0;i<bm_handle.bm_y;i++)  for (j=0;j<bm_handle.bm_x;j++)
    fprintf(plik,"%c%c%c", 
		    pix[(y*j+i)+2], pix[(y*j+i)+1], pix[(y*j+i)]);
 fclose(plik);
 printf("Bitmap: %s written. (currently broken conversion)\n", temp);
}

unsigned char *load_pcx_file(char *filename, PcxHeader *pcxHeader)
{
 int idx;
 int c;                             
 int i;                             
 int numRepeat;      
 FILE *filePtr;                
 int width;                         
 int height;                        
 unsigned char *pixelData;     
 unsigned char *paletteData;   
 idx = 0;
 filePtr = fopen(filename, "rb"); 
 if (!filePtr) return NULL;	
 c = getc(filePtr);
 if (c != 10)
   {
    fclose(filePtr);
    return NULL;
   } 
 c = getc(filePtr);
 if (c != 5)
   {
    fclose(filePtr);
    return NULL;
   } 
 rewind(filePtr);	
 fgetc(filePtr);
 fgetc(filePtr);
 fgetc(filePtr);
 fgetc(filePtr);     
 pcxHeader->xMin = fgetc(filePtr);       
 pcxHeader->xMin |= fgetc(filePtr) << 8; 
 pcxHeader->yMin = fgetc(filePtr);       
 pcxHeader->yMin |= fgetc(filePtr) << 8; 
 pcxHeader->xMax = fgetc(filePtr);       
 pcxHeader->xMax |= fgetc(filePtr) << 8; 
 pcxHeader->yMax = fgetc(filePtr);       
 pcxHeader->yMax |= fgetc(filePtr) << 8; 
 width = pcxHeader->xMax - pcxHeader->xMin + 1;
 height = pcxHeader->yMax - pcxHeader->yMin + 1;
 pixelData = (unsigned char*)malloc(width*height);
 fseek(filePtr, 128, SEEK_SET);
 while (idx < (width*height))
   {
    c = getc(filePtr);
    if (c > 0xbf)
      {
	numRepeat = 0x3f & c;
	c = getc(filePtr);
	for (i = 0; i < numRepeat; i++) pixelData[idx++] = c;
      }
    else pixelData[idx++] = c;
    fflush(stdout);
   }
 paletteData = (unsigned char*)malloc(768);
 fseek(filePtr, -769, SEEK_END);     
 c = getc(filePtr);
 if (c != 12)
   {
    fclose(filePtr);
    return NULL;
   } 
 for (i = 0; i < 768; i++)
   {
    c = getc(filePtr);
    paletteData[i] = c;
    }	
 fclose(filePtr);
 pcxHeader->palette = paletteData;
/* write_as_bmp(filename, pixelData, width, height);*/
 return pixelData;
}



Texture *load_pcx_texture(char *filename)
{
     PcxHeader texInfo;            
     Texture *thisTexture;      
     unsigned char *unscaledData;
     int i;                             
     int j;                             
     int width;                         
     int height;                        
     thisTexture = (Texture*)malloc(sizeof(Texture));
     if (thisTexture == NULL)
          return NULL;
     thisTexture->data = load_pcx_file(filename, &texInfo);
     if (thisTexture->data == NULL)
     {
          free(thisTexture->data);
          return NULL;
     }
     thisTexture->palette = texInfo.palette;
     thisTexture->width = texInfo.xMax - texInfo.xMin + 1;
     thisTexture->height = texInfo.yMax - texInfo.yMin + 1;
     thisTexture->textureType = PCX;

     
     unscaledData = (unsigned char*)malloc(thisTexture->width*thisTexture->height*4);

     
     for (j = 0; j < thisTexture->height; j++) 
     {
          for (i = 0; i < thisTexture->width; i++) 
          {
               unscaledData[4*(j*thisTexture->width+i)+0] = (unsigned char)thisTexture->palette[3*thisTexture->data[j*thisTexture->width+i]+0];
               unscaledData[4*(j*thisTexture->width+i)+1] = (unsigned char)thisTexture->palette[3*thisTexture->data[j*thisTexture->width+i]+1];
               unscaledData[4*(j*thisTexture->width+i)+2] = (unsigned char)thisTexture->palette[3*thisTexture->data[j*thisTexture->width+i]+2];
               unscaledData[4*(j*thisTexture->width+i)+3] = (unsigned char)255;
          }
     }

     
     width = thisTexture->width;
     height = thisTexture->height;

     
     i = 0;
     while (width)
     {
          width /= 2;
          i++;
     }
     thisTexture->scaledHeight = (long)pow(2, i-1);

     
     i = 0;
     while (height)
     {
          height /= 2;
          i++;
     }
     thisTexture->scaledWidth = (long)pow(2, i-1);

     
     if (thisTexture->data != NULL)
     {
          free(thisTexture->data);
          thisTexture->data = NULL;
     }

     
     thisTexture->data = (unsigned char*)malloc(thisTexture->scaledWidth*thisTexture->scaledHeight*4);
     
     
     gluScaleImage (GL_RGBA, thisTexture->width, thisTexture->height, GL_UNSIGNED_BYTE, unscaledData, thisTexture->scaledWidth, thisTexture->scaledHeight, GL_UNSIGNED_BYTE, thisTexture->data);

     return thisTexture;
}


Texture *load_bmp_texture(char *filename)
{
 FILE* plik;
 int i,j;
 char b,g,r;
 char B,M;
 BMPTag bm_handle;
 Texture* t;
 plik = fopen(filename, "rb");
 if (!plik) return NULL;
 init_bmp(&bm_handle);
 i = fscanf(plik,"%c%c",&B,&M);
 if (i != 2) return NULL;
 if (B != 'B' || M != 'M') return NULL;
/* printf("reading bitmap...\n");*/
 fread(&bm_handle.fsize,4,1,plik);
 fread(&bm_handle.dummy,4,1,plik);
 fread(&bm_handle.offset,4,1,plik);
 fread(&bm_handle.dummy2,4,1,plik);
 fread(&bm_handle.bm_x,4,1,plik);
 fread(&bm_handle.bm_y,4,1,plik);
 fread(&bm_handle.planes,2,1,plik);
 fread(&bm_handle.bpp,2,1,plik);
 if (bm_handle.bpp != 24) return NULL;
 fseek(plik,bm_handle.offset,SEEK_SET);
 t = (Texture*)malloc(sizeof(Texture));
 t->data   = (unsigned char*)malloc(3*bm_handle.bm_y*bm_handle.bm_x*sizeof(unsigned char));
 t->width  = bm_handle.bm_x;
 t->height = bm_handle.bm_y;
 for (i=0;i<bm_handle.bm_y;i++)  for (j=0;j<bm_handle.bm_x;j++)
    {
     fscanf(plik,"%c%c%c", &b,&g,&r);
     t->data[3*(t->height * i + j)   ] = r;
     t->data[3*(t->height * i + j) +1] = g;
     t->data[3*(t->height * i + j) +2] = b;
    }
 fclose(plik);
 t->palette = NULL;
 t->scaledHeight = 0;
 t->scaledWidth = 0;
 t->textureType = BMP;
 return t;
}



Texture *load_texture(char *filename)
{
	Texture *thisTexture;
	char *extStr;
        thisTexture = NULL;
	
	extStr = strchr(filename, '.');
	extStr++;

	
	if ((strcmp(extStr, "BMP") == 0) || (strcmp(extStr, "bmp") == 0))
		thisTexture = load_bmp_texture(filename);
	else if ((strcmp(extStr, "PCX") == 0) || (strcmp(extStr, "pcx") == 0) )
		thisTexture = load_pcx_texture(filename);
	return thisTexture;
}


void setup_md2_texture(Texture *thisTexture)
{
	glGenTextures(1, &thisTexture->texID);
	glBindTexture(GL_TEXTURE_2D, thisTexture->texID);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);

	switch (thisTexture->textureType)
	{
		case BMP:
			gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, thisTexture->width, thisTexture->height, 
							  GL_RGB, GL_UNSIGNED_BYTE, thisTexture->data);
			break;
		case PCX:
			gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, thisTexture->width, thisTexture->height,
							  GL_RGBA, GL_UNSIGNED_BYTE, thisTexture->data);
		default:
			break;
	}
}




ModelData *load_md2(char *filename, char *textureName)
{
	FILE *filePtr;						
	int fileLen;						
    char *buffer;						
		
	ModelData *model;					
	ModelHeader *modelHeader;			
	Texture *md2Texture;				

	StIndex *stPtr;					
    Frame *frame;						
	Vector *pointListPtr;				
    Mesh *triIndex, *bufIndexPtr;		
    int i, j;							

	
	filePtr = fopen(filename, "rb");
	if (filePtr == NULL)
		return NULL;

	
    fseek(filePtr, 0, SEEK_END);
    fileLen = ftell(filePtr);
    fseek(filePtr, 0, SEEK_SET);
	
    
    buffer = (char*)malloc(fileLen + 2);
    fread(buffer, sizeof(char), fileLen, filePtr);

	
    modelHeader = (ModelHeader*)buffer;

	
   	model = (ModelData*)malloc(sizeof(ModelData));
	if (model == NULL)
		return NULL;

	
    model->pointList = (Vector*)malloc(sizeof(Vector)*modelHeader->numXYZ * modelHeader->numFrames+1);

	
    model->numPoints = modelHeader->numXYZ;
    model->numFrames = modelHeader->numFrames;
    model->frameSize = modelHeader->framesize;

/*    printf("model initialized.\n");    */
    for(j = 0; j < modelHeader->numFrames; j++)
    {
       
       frame = (Frame*)(&buffer[modelHeader->offsetFrames + modelHeader->framesize * j]);
       pointListPtr = (Vector*)(&model->pointList[modelHeader->numXYZ * j]);
       for(i = 0; i < modelHeader->numXYZ; i++)
       {
          pointListPtr[i].point[0] = frame->scale[0] * frame->fp[i].v[0] + frame->translate[0];
          pointListPtr[i].point[1] = frame->scale[1] * frame->fp[i].v[1] + frame->translate[1];
          pointListPtr[i].point[2] = frame->scale[2] * frame->fp[i].v[2] + frame->translate[2];
       }
    }
			 
/*   printf("frames read.\n");*/
	md2Texture = load_texture(textureName);
	if (md2Texture != NULL)
	{
		
		setup_md2_texture(md2Texture);
		model->modelTex = md2Texture;
	}
	else
		return NULL;
/*	printf("texture loaded.\n");*/

    
    model->st = (TexCoord*)malloc(sizeof(TexCoord)*modelHeader->numST);

	
    model->numST = modelHeader->numST;

	
    stPtr = (StIndex*)&buffer[modelHeader->offsetST];

	
    for (i = 0; i < modelHeader->numST; i++)
    {
		model->st[i].s = (float)stPtr[i].s / (float)md2Texture->width;
        model->st[i].t = (float)stPtr[i].t / (float)md2Texture->height;
    }

	
	triIndex = (Mesh*)malloc(sizeof(Mesh) * modelHeader->numTris);

	
	model->numTriangles = modelHeader->numTris;
	model->triIndex = triIndex;
	
	
	bufIndexPtr = (Mesh*)&buffer[modelHeader->offsetTris];

	
	for (j = 0; j < model->numFrames; j++)		
	{
		
		for(i = 0; i < modelHeader->numTris; i++)
		{
		   triIndex[i].meshIndex[0] = bufIndexPtr[i].meshIndex[0];
		   triIndex[i].meshIndex[1] = bufIndexPtr[i].meshIndex[1];
		   triIndex[i].meshIndex[2] = bufIndexPtr[i].meshIndex[2];
		   triIndex[i].stIndex[0] = bufIndexPtr[i].stIndex[0];
		   triIndex[i].stIndex[1] = bufIndexPtr[i].stIndex[1];
		   triIndex[i].stIndex[2] = bufIndexPtr[i].stIndex[2];
		}
	}

	
	fclose(filePtr);
    free(buffer);

	model->currentFrame = 0;
	model->nextFrame = 1;
	model->interpol = 0.0;

    return model;
}



void calc_normal( float *p1, float *p2, float *p3 )
{
   float a[3], b[3], result[3];
   float length;

   a[0] = p1[0] - p2[0];
   a[1] = p1[1] - p2[1];
   a[2] = p1[2] - p2[2];

   b[0] = p1[0] - p3[0];
   b[1] = p1[1] - p3[1];
   b[2] = p1[2] - p3[2];

   result[0] = a[1] * b[2] - b[1] * a[2];
   result[1] = b[0] * a[2] - a[0] * b[2];
   result[2] = a[0] * b[1] - b[0] * a[1];

   
   length = (float)sqrt(result[0]*result[0] + result[1]*result[1] + result[2]*result[2]);

   
   glNormal3f(result[0]/length, result[1]/length, result[2]/length);
}

FILE* get_next_dat_file()
{
 char temp[512];
 char outf[512];
 int i,done;
 FILE* out;
 strcpy(temp, strModel);
 if (strstr(temp, "."))
   {
    i = 0;
    while (temp[i] != '.') i++;
    temp[i] = 0;
   }
 done = 0;
 i = 0;
 do
  {
   i ++;
   sprintf(outf, "%s_%d.dat", temp, i);
   out = fopen(outf, "r");
   if (out) { fclose(out);  continue; }
   out = fopen(outf, "w");
   if (!out) return NULL;
   done = 1;
  }
 while (!done);
 printf("Writing to: %s\n", outf);
 return out;
}

void write_dat_header(FILE* p, int nt)
{
 int ntex, r, b;
 float mins,maxs,amb;
 float lookz, obsz;
 float lx,ly,lz;
 ntex = 36;
 mins = .1;
 maxs = .6;
 amb  = .3;
 lookz = 35.;
 obsz = -300.;
 r = 6;
 b = 128;
 lx = 50.;
 ly = 100.;
 lz = -250.;
 if (hdr)
  {
   fprintf(p,"Screen: (800,600)\n");
/*printf(" NormalDistorber: %f%%\n",d_rand());*/ 
   fprintf(p,"MinShadow: %f\n", mins);
   fprintf(p,"MaxShadow: %f\n", maxs);
   fprintf(p,"Ambient: %f\n", amb);
   fprintf(p,"MaxRecurse: %d\n", r);
   fprintf(p,"Backup: %d\n", b);
   fprintf(p,"Observer: Vertex: (%f,%f,%f)\n", 0., 0., obsz);
   fprintf(p,"LookZ: %f%%\n", lookz);	
   fprintf(p,"Light: Vertex: (%f,%f,%f)\n", lx, ly, lz);	
   fprintf(p,"TexDirectory: textures\n");
   fprintf(p,"NumTextures: %d\n",ntex);
   fprintf(p,"nTriangles: %d\n",nt);
   fprintf(p,"ListTransform: [%d,%d]\n",0, nt-1);
   fprintf(p,"{\n");
   fprintf(p," RotateX: -90\n");
   fprintf(p," RotateZ: 90\n");
   fprintf(p," Scale: (2,2,2)\n");
   fprintf(p,"\n");
   fprintf(p,"} \n");
  }
}

void write_triangle(FILE* p, int idx, float s1, float t1, float s2, float t2, float s3, float t3, 
		float x1, float y1, float z1, 
		float x2, float y2, float z2, 
		float x3, float y3, float z3)
{
    fprintf(p, "Triangle: %d\n",idx);
    fprintf(p, "{\n");
    fprintf(p, " a: Vertex: (%f,%f,%f)\n", x1, y1, z1);
    fprintf(p, " b: Vertex: (%f,%f,%f)\n", x2, y2, z2);
    fprintf(p, " c: Vertex: (%f,%f,%f)\n", x3, y3, z3);
    fprintf(p, " texA: TexCoord: (%f,%f)\n", s1, t1);
    fprintf(p, " texB: TexCoord: (%f,%f)\n", s2, t2);
    fprintf(p, " texC: TexCoord: (%f,%f)\n", s3, t3);
    fprintf(p, " na: Vector: (%f,%f,%f)\n", 0., 0., 0.);
    fprintf(p, " nb: Vector: (%f,%f,%f)\n", 0., 0., 0.);
    fprintf(p, " nc: Vector: (%f,%f,%f)\n", 0., 0., 0.);
    fprintf(p, " transparency: RGB: (%f,%f,%f)\n",tr,tg,tb);
    fprintf(p, " specular: RGB: (%f,%f,%f)\n",sr,sg,sb);
    fprintf(p, " diffuse: RGB: (%f,%f,%f)\n",cr,cg,cb);
    fprintf(p, " transparencyFactR: (1,%f)\n", fr);
    fprintf(p, " transparencyFactG: (1,%f)\n", fg);
    fprintf(p, " transparencyFactB: (1,%f)\n", fb);
    fprintf(p, " normalDist: %f%%\n",nd);
    fprintf(p, " specularFact: %f\n", sf);
    fprintf(p, " faces: %d\n",fac);		
    fprintf(p, " texture: %d\n", tid);
}

void convert_md2(ModelData *model, int startFrame, int endFrame, float percent)
{
	Vector *pointList;			
	Vector *nextPointList;		
	int i;							
	float x1, y1, z1;				
	float x2, y2, z2;				
	FILE* p;
	Vector vertex[3];	
	if (model == NULL) return;
	if ((startFrame < 0) || (endFrame < 0)) return;
	if ((startFrame >= model->numFrames) || (endFrame >= model->numFrames)) return;
	if (model->interpol >= 1.0)
	{
	 model->interpol = 0.0f;
	 model->currentFrame++;
	 if (model->currentFrame >= endFrame) model->currentFrame = startFrame; 
	 model->nextFrame = model->currentFrame + 1;
	 if (model->nextFrame >= endFrame) model->nextFrame = startFrame;
	}
	pointList = &model->pointList[model->numPoints*model->currentFrame];
	nextPointList = &model->pointList[model->numPoints*model->nextFrame];
	p = get_next_dat_file();
	if (!p) return;
/*	glBindTexture(GL_TEXTURE_2D, model->modelTex->texID);*/
/*	glBegin(GL_TRIANGLES);*/
	write_dat_header(p, model->numTriangles);
	for (i = 0; i < model->numTriangles; i++)
	  {
/*           printf("idx = %d\n", model->triIndex[i].meshIndex[0]);*/
	  }
	
	for (i = 0; i < model->numTriangles; i++)
		{
		 x1 = pointList[model->triIndex[i].meshIndex[0]].point[0];
		 y1 = pointList[model->triIndex[i].meshIndex[0]].point[1];
		 z1 = pointList[model->triIndex[i].meshIndex[0]].point[2];
		 x2 = nextPointList[model->triIndex[i].meshIndex[0]].point[0];
		 y2 = nextPointList[model->triIndex[i].meshIndex[0]].point[1];
		 z2 = nextPointList[model->triIndex[i].meshIndex[0]].point[2];			
		 vertex[0].point[0] = x1 + model->interpol * (x2 - x1);
		 vertex[0].point[1] = y1 + model->interpol * (y2 - y1);
		 vertex[0].point[2] = z1 + model->interpol * (z2 - z1);
		
			
		 x1 = pointList[model->triIndex[i].meshIndex[2]].point[0];
		 y1 = pointList[model->triIndex[i].meshIndex[2]].point[1];
		 z1 = pointList[model->triIndex[i].meshIndex[2]].point[2];
		 x2 = nextPointList[model->triIndex[i].meshIndex[2]].point[0];
		 y2 = nextPointList[model->triIndex[i].meshIndex[2]].point[1];
		 z2 = nextPointList[model->triIndex[i].meshIndex[2]].point[2];

			
		 vertex[2].point[0] = x1 + model->interpol * (x2 - x1);
		 vertex[2].point[1] = y1 + model->interpol * (y2 - y1);
		 vertex[2].point[2] = z1 + model->interpol * (z2 - z1);	
	
			
		 x1 = pointList[model->triIndex[i].meshIndex[1]].point[0];
		 y1 = pointList[model->triIndex[i].meshIndex[1]].point[1];
		 z1 = pointList[model->triIndex[i].meshIndex[1]].point[2];
		 x2 = nextPointList[model->triIndex[i].meshIndex[1]].point[0];
		 y2 = nextPointList[model->triIndex[i].meshIndex[1]].point[1];
		 z2 = nextPointList[model->triIndex[i].meshIndex[1]].point[2];

			
		 vertex[1].point[0] = x1 + model->interpol * (x2 - x1);
		 vertex[1].point[1] = y1 + model->interpol * (y2 - y1);
		 vertex[1].point[2] = z1 + model->interpol * (z2 - z1);

			
		 calc_normal(vertex[0].point, vertex[2].point, vertex[1].point);
/*		 Here normal is computed, can write 0 to .dat file*/

			
/*		glTexCoord2f(model->st[model->triIndex[i].stIndex[0]].s, model->st[model->triIndex[i].stIndex[0]].t);*/
/*		glVertex3fv(vertex[0].point);*/
		write_triangle(p, idx, 
				model->st[model->triIndex[i].stIndex[0]].s, 
				model->st[model->triIndex[i].stIndex[0]].t,
				model->st[model->triIndex[i].stIndex[2]].s, 
				model->st[model->triIndex[i].stIndex[2]].t,
				model->st[model->triIndex[i].stIndex[1]].s, 
				model->st[model->triIndex[i].stIndex[1]].t,
				vertex[0].point[0], vertex[0].point[1], vertex[0].point[2],
				vertex[2].point[0], vertex[2].point[1], vertex[2].point[2],
				vertex[1].point[0], vertex[1].point[1], vertex[1].point[2]
				);
		idx++;

/*		glTexCoord2f(model->st[model->triIndex[i].stIndex[2]].s ,model->st[model->triIndex[i].stIndex[2]].t);*/
/*		glVertex3fv(vertex[2].point);*/

/*		glTexCoord2f(model->st[model->triIndex[i].stIndex[1]].s,model->st[model->triIndex[i].stIndex[1]].t);*/
/*		glVertex3fv(vertex[1].point);*/
		}
/*	glEnd();*/
	model->interpol += percent;	
	fclose(p);
}


void display_md2(ModelData *model, int startFrame, int endFrame, float percent)
{
	Vector *pointList;			
	Vector *nextPointList;		
	int i;							
	float x1, y1, z1;				
	float x2, y2, z2;				

	Vector vertex[3];	

	if (model == NULL)
		return;
	
	if ( (startFrame < 0) || (endFrame < 0) )
		return;

	if ( (startFrame >= model->numFrames) || (endFrame >= model->numFrames) )
		return;

	if (model->interpol >= 1.0)
	{
		model->interpol = 0.0f;
		model->currentFrame++;
		if (model->currentFrame >= endFrame)
			model->currentFrame = startFrame; 

		model->nextFrame = model->currentFrame + 1;

		if (model->nextFrame >= endFrame)
			model->nextFrame = startFrame;
	}

	pointList = &model->pointList[model->numPoints*model->currentFrame];
	nextPointList = &model->pointList[model->numPoints*model->nextFrame];
	
	glBindTexture(GL_TEXTURE_2D, model->modelTex->texID);
	glBegin(GL_TRIANGLES);
		for (i = 0; i < model->numTriangles; i++)
		{
			
			x1 = pointList[model->triIndex[i].meshIndex[0]].point[0];
			y1 = pointList[model->triIndex[i].meshIndex[0]].point[1];
			z1 = pointList[model->triIndex[i].meshIndex[0]].point[2];
			x2 = nextPointList[model->triIndex[i].meshIndex[0]].point[0];
			y2 = nextPointList[model->triIndex[i].meshIndex[0]].point[1];
			z2 = nextPointList[model->triIndex[i].meshIndex[0]].point[2];

			
			vertex[0].point[0] = x1 + model->interpol * (x2 - x1);
			vertex[0].point[1] = y1 + model->interpol * (y2 - y1);
			vertex[0].point[2] = z1 + model->interpol * (z2 - z1);
		
			
			x1 = pointList[model->triIndex[i].meshIndex[2]].point[0];
			y1 = pointList[model->triIndex[i].meshIndex[2]].point[1];
			z1 = pointList[model->triIndex[i].meshIndex[2]].point[2];
			x2 = nextPointList[model->triIndex[i].meshIndex[2]].point[0];
			y2 = nextPointList[model->triIndex[i].meshIndex[2]].point[1];
			z2 = nextPointList[model->triIndex[i].meshIndex[2]].point[2];

			
			vertex[2].point[0] = x1 + model->interpol * (x2 - x1);
			vertex[2].point[1] = y1 + model->interpol * (y2 - y1);
			vertex[2].point[2] = z1 + model->interpol * (z2 - z1);	
	
			
			x1 = pointList[model->triIndex[i].meshIndex[1]].point[0];
			y1 = pointList[model->triIndex[i].meshIndex[1]].point[1];
			z1 = pointList[model->triIndex[i].meshIndex[1]].point[2];
			x2 = nextPointList[model->triIndex[i].meshIndex[1]].point[0];
			y2 = nextPointList[model->triIndex[i].meshIndex[1]].point[1];
			z2 = nextPointList[model->triIndex[i].meshIndex[1]].point[2];

			
			vertex[1].point[0] = x1 + model->interpol * (x2 - x1);
			vertex[1].point[1] = y1 + model->interpol * (y2 - y1);
			vertex[1].point[2] = z1 + model->interpol * (z2 - z1);

			
			calc_normal(vertex[0].point, vertex[2].point, vertex[1].point);

			
			glTexCoord2f(model->st[model->triIndex[i].stIndex[0]].s,
						 model->st[model->triIndex[i].stIndex[0]].t);
			glVertex3fv(vertex[0].point);

			glTexCoord2f(model->st[model->triIndex[i].stIndex[2]].s ,
					     model->st[model->triIndex[i].stIndex[2]].t);
			glVertex3fv(vertex[2].point);

			glTexCoord2f(model->st[model->triIndex[i].stIndex[1]].s,
					     model->st[model->triIndex[i].stIndex[1]].t);
			glVertex3fv(vertex[1].point);
		}
	glEnd();

	model->interpol += percent;	

}

void clean_up()
{
 free_model(myModel);
}



void init()
{
 glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
 glShadeModel(GL_SMOOTH);					
 glEnable(GL_DEPTH_TEST);					
 glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	
 glEnable(GL_LIGHTING);					
 glEnable(GL_LIGHT0);
 glEnable(GL_COLOR_MATERIAL);
 glEnable(GL_TEXTURE_2D);	
 printf("loading model..\n");
 myModel = load_md2(strModel, strTexture);
 if (!myModel) { printf("Cant load model.\n"); exit(1); }
 printf("model loaded.\n");
 angx = 270.;
 angy = angz = 0.;
 scale = 1.;
 percent = 0.03;
}

void time_counter()
{
 char tstr[64];
 if (t1 == (time_t)0)
   {
    time(&t1);
    time(&t2);
    return;
   }
 fps++;
 time(&t2);
 if (t2 > t1)
   {
    sprintf(tstr, "Display MD2 %d FPS", fps);
    t1 = t2;
    glutSetWindowTitle(tstr);
    fps = 0;
   }
}


void render_scene()
{
 radians =  (float)(PI*(angle-90.0f)/180.0f);
/* percent = 0.005;	*/
/* percent = 0.03;	*/
 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		
 glLoadIdentity();
 gluLookAt(0., 0., -100., 0., 0., 0., 0.0, 1.0, 0.0);
 glPushMatrix();
 glRotatef(angx, 1.0f, 0.0f, 0.0f);
 glRotatef(angy, 0.0f, 1.0f, 0.0f);
 glRotatef(angz, 0.0f, 0.0f, 1.0f);
 glScalef(scale, scale, scale);
 display_md2(myModel, 0, 40, percent);
 glColor3f(1.0, 1.0, 1.0);
 glPopMatrix();
 time_counter();
 glFlush();
 glutSwapBuffers();
}

void resize_scene(int w, int h)
{
 glViewport(0, 0, w, h);		
 glMatrixMode(GL_PROJECTION);			
 glLoadIdentity();					
 gluPerspective(54.0f,(GLfloat)w/(GLfloat)h,1.0f,1000.0f);
 glMatrixMode(GL_MODELVIEW);				
 glLoadIdentity();						
}

void anim(void)
{
    glutPostRedisplay();
}


void visible(int vis)
{
    if (vis == GLUT_VISIBLE) glutIdleFunc(anim);
    else                     glutIdleFunc(NULL);
}

void help()
{
 printf("usage: md22dat file.md2 file.{bmp|pcx}\n");
 printf("optional: f(rgb), t(rgb), s(rgb), c(rgb), sf\n");
 printf("nd[in %%], fac, tid, idx, hdr\n");
 printf("keys:\n");
 printf("SPACE \t write current frame to .DAT file\n");
 printf("Q/ESC \t quit\n");
 printf("-=    \t manipulate speed\n");
 printf("12    \t manipulate scale\n");
 printf("xXyYzZ\t manipulate angles\n");
 printf("d\t write header toggle\n");
}

void keyboard(unsigned char key, int x, int y)
{
/* percent = 0.005;*/
    printf("percent = %f, angle = (%f,%f,%f), scale = %f\n", percent, angx, angy, angz, scale);
 switch (key)
   {
    case 27:  case 'q': clean_up(); exit(0); break;
    case ' ': convert_md2(myModel, 0, 40, percent); break;
    case 'h': help(); break;
    case 'd': hdr = ! hdr; break;
    case '=': percent *= 1.05; break;
    case '-': percent /= 1.05; break;
    case 'x': angx -= 3.; break;
    case 'y': angy -= 3.; break;
    case 'z': angz -= 3.; break;
    case 'X': angx += 3.; break;
    case 'Y': angy += 3.; break;
    case 'Z': angz += 3.; break;
    case '1': scale *= 1.05; break;
    case '2': scale /= 1.05; break;
   }
}

void set_defaults()
{
 fr = 1.05;
 fg = 1.08;
 fb = 1.10;
 tr = 1.;
 tg = 1.;
 tb = 1.;
 sr = 1.;
 sg = 1.;
 sb = 1.;
 cr = 1.;
 cg = 1.;
 cb = 1.;
 sf = 150.;
 nd = 0.;
 fac = 1;
 tid = 2;
 idx = 0;
 hdr = 1;
}

void load_from_cmdline(int lb, char** par)
{
 fr = atof(par[3]);
 fg = atof(par[4]);
 fb = atof(par[5]);
 tr = atof(par[6]);
 tg = atof(par[7]);
 tb = atof(par[8]);
 sr = atof(par[9]);
 sg = atof(par[10]);
 sb = atof(par[11]);
 cr = atof(par[12]);
 cg = atof(par[13]);
 cb = atof(par[14]);
 sf = atof(par[15]);
 nd = atof(par[16]);
 fac = atoi(par[17]);
 tid = atoi(par[18]);
 idx = atoi(par[19]);
 hdr = atoi(par[20]);
}


int main(int lb, char** par)
{
 if (lb < 2) { help(); exit(1); }
 else if (lb == 2) strcpy(strTexture, "default.pcx");
 else if (lb == 3) set_defaults();
 else if (lb == 21) load_from_cmdline(lb, par);
 else { help(); exit(1); }
 strcpy(strModel, par[1]);
 if (lb >= 3) strcpy(strTexture, par[2]);
 glutInit(&lb, par);
 glutInitDisplayMode(GLUT_DOUBLE);
 glutInitWindowSize(800, 600);
 glutInitWindowPosition(100, 100);
 glutCreateWindow(par[0]);
 init();
 glutDisplayFunc(render_scene);
 glutReshapeFunc(resize_scene);
 glutKeyboardFunc(keyboard);
 glutVisibilityFunc(visible);
 glutMainLoop();
 return 0;
}

