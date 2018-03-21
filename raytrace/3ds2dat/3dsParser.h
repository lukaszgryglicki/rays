#ifndef _3DS_PARSER__MDBMA_
#define _3DS_PARSER__MDBMA_
#include "Main.h"
#include "Vector.h"

#define PRIMARY       0x4D4D
#define OBJECTINFO    0x3D3D				 
#define VERSION       0x0002				
#define EDITKEYFRAME  0xB000				
#define MATERIAL	  0xAFFF				
#define OBJECT		  0x4000				
#define MATNAME       0xA000				
#define MATDIFFUSE    0xA020				
#define MATMAP        0xA200				
#define MATMAPFILE    0xA300				
#define OBJECT_MESH   0x4100				
#define OBJECT_VERTICES     0x4110			
#define OBJECT_FACES		0x4120			
#define OBJECT_MATERIAL		0x4130			
#define OBJECT_UV			0x4140			

typedef struct Face
{
 int vertIndex[3];							
 int coordIndex[3];							
} Face;

typedef struct MaterialInfo
{
 char  strName[255];							
 char  strFile[255];							
 unsigned char  color[3];	/* BYTE */
 int   texureId;								
 float uTile;								
 float vTile;								
 float uOffset;								
 float vOffset;									
} MaterialInfo;

typedef struct Object3D
{
	int  numOfVerts;			
	int  numOfFaces;			
	int  numTexVertex;			
	int  materialID;			
	bool bHasTexture;			
	char strName[255];			
	unsigned int *pIndices;		
	Vector3  *pVerts;			
	Vector3  *pNormals;		
	Vector2  *pTexVerts;		
	Face *pFaces;				
} Object3D;

typedef struct Model3D
{
	int numOfObjects;							
	int numOfMaterials;							
	vector<MaterialInfo> pMaterials;			
	vector<Object3D> pObject;					
} Model3D;

typedef struct Indices 
{							
 unsigned short a, b, c, bVisible;		
} Indices;

typedef struct Chunk
{
 unsigned short int ID;					
 unsigned int length;					
 unsigned int bytesRead;					
} Chunk;


class Parser3DS
{
 public:
	Parser3DS();										
	bool Import3DS(Model3D *pModel, char *strFileName);

private:
	int GetString(char *);
	void ReadChunk(Chunk *);
	void ProcessNextChunk(Model3D *pModel, Chunk *);
	void ProcessNextObjectChunk(Model3D *pModel, Object3D *pObject, Chunk *);
	void ProcessNextMaterialChunk(Model3D *pModel, Chunk *);
	void ReadColorChunk(MaterialInfo *pMaterial, Chunk *pChunk);
	void ReadVertices(Object3D *pObject, Chunk *);
	void ReadVertexIndices(Object3D *pObject, Chunk *);
	void ReadUVCoordinates(Object3D *pObject, Chunk *);
	void ReadObjectMaterial(Model3D *pModel, Object3D *pObject, Chunk *pPreviousChunk);
	void ComputeNormals(Model3D *pModel);
	void CleanUp();
	
	FILE *m_FilePointer;
	Chunk *m_CurrentChunk;
	Chunk *m_TempChunk;
	int* ibuffer;
};


#endif

