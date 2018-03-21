#include "Vector.h"
#include "3dsParser.h"

Parser3DS::Parser3DS()
{
 //printf("Constructor::\n");
 m_CurrentChunk = new Chunk;				 
 m_TempChunk    = new Chunk;	
 ibuffer = new int[UNKNOWN_CHUNK_BUFSIZE];
}


bool Parser3DS::Import3DS(Model3D *pModel, char *strFileName)
{
 printf("Opening modelfile %s\n", strFileName);
 //char strMessage[255] = {0};
 m_FilePointer = fopen(strFileName, "rb");
 if(!m_FilePointer) 
	{
		printf("Unable to find the file: %s\n", strFileName);
		return false;
	}
 //printf("File opened.\n");
 //printf("reading first chunk.\n");
 ReadChunk(m_CurrentChunk);	
 //printf("First chunk read.\n");
 if (m_CurrentChunk->ID != PRIMARY)
	{
		printf("Unable to load PRIMARY chuck from file: %s\n", strFileName);
		return false;
	}	
 //printf("First chunk ok.\nProcessing next chunk(s)\n");
 ProcessNextChunk(pModel, m_CurrentChunk);	
 printf("All chunks processed.\nComputing normals\n");
 ComputeNormals(pModel);	
 printf("Normals computed.\nCleaning up.\n");
 CleanUp();	
 printf("Model imported.\n");
 return true;
}

void Parser3DS::CleanUp()
{
	fclose(m_FilePointer);						 
	delete m_CurrentChunk;						
	delete m_TempChunk;			
	delete[] ibuffer;
}


void Parser3DS::ProcessNextChunk(Model3D *pModel, Chunk *pPreviousChunk)
{
	//printf("Processing next chunk: model: %p, chunk: %p\n", pModel, pPreviousChunk);
	Object3D newObject = {0};					 /* FIXME */
	MaterialInfo newTexture = {0};				
	unsigned int version = 0;
	//printf("Previous chunk ID: %d\n", pPreviousChunk->ID);
	m_CurrentChunk = new Chunk;		
	//printf("New chunk pointer: %p\n", m_CurrentChunk);
	//printf("While not read all bytes.\n");
	while (pPreviousChunk->bytesRead < pPreviousChunk->length)
	{
		//printf("Read %d of %d bytes\n", pPreviousChunk->bytesRead, pPreviousChunk->length);
		ReadChunk(m_CurrentChunk);
		//printf("Case ID of %x\n", m_CurrentChunk->ID);
		switch (m_CurrentChunk->ID)
		{
		case VERSION:							
			//printf("Got Version Chunk.\n");
			//printf("will read %d bytes\n", m_CurrentChunk->length - m_CurrentChunk->bytesRead);
			m_CurrentChunk->bytesRead += (unsigned int)fread(&version, 1, m_CurrentChunk->length - m_CurrentChunk->bytesRead, m_FilePointer);			
			if (version > 0x03)
				printf("This 3DS file is over version 3 so it may load incorrectly\n");
			break;
			
		case OBJECTINFO:						
			ReadChunk(m_TempChunk);			
			m_TempChunk->bytesRead += (unsigned int)fread(&version, 1, m_TempChunk->length - m_TempChunk->bytesRead, m_FilePointer);			
			m_CurrentChunk->bytesRead += m_TempChunk->bytesRead;
			ProcessNextChunk(pModel, m_CurrentChunk);
			break;

		case MATERIAL:						
			pModel->numOfMaterials++;
			pModel->pMaterials.push_back(newTexture);
			ProcessNextMaterialChunk(pModel, m_CurrentChunk);
			break;			

		case OBJECT:							
			pModel->numOfObjects++;
			pModel->pObject.push_back(newObject);
			memset(&(pModel->pObject[pModel->numOfObjects - 1]), 0, sizeof(Object3D));
			m_CurrentChunk->bytesRead += GetString(pModel->pObject[pModel->numOfObjects - 1].strName);
			ProcessNextObjectChunk(pModel, &(pModel->pObject[pModel->numOfObjects - 1]), m_CurrentChunk);
			break;
			
		case EDITKEYFRAME:			
			m_CurrentChunk->bytesRead += (unsigned int)fread(ibuffer, 1, m_CurrentChunk->length - m_CurrentChunk->bytesRead, m_FilePointer);
			break;
		default: 
			
			m_CurrentChunk->bytesRead += (unsigned int)fread(ibuffer, 1, m_CurrentChunk->length - m_CurrentChunk->bytesRead, m_FilePointer);
			break;
		}	
		pPreviousChunk->bytesRead += m_CurrentChunk->bytesRead;
	}
	//printf("Finally Read %d of %d bytes\n", pPreviousChunk->bytesRead, pPreviousChunk->length);
	
	delete m_CurrentChunk;
	m_CurrentChunk = pPreviousChunk;
}


void Parser3DS::ProcessNextObjectChunk(Model3D *pModel, Object3D *pObject, Chunk *pPreviousChunk)
{
 m_CurrentChunk = new Chunk;	
 while (pPreviousChunk->bytesRead < pPreviousChunk->length)
	{
	 ReadChunk(m_CurrentChunk);	
	 switch (m_CurrentChunk->ID)
		{
		case OBJECT_MESH:		
			ProcessNextObjectChunk(pModel, pObject, m_CurrentChunk);
			break;

		case OBJECT_VERTICES:			 
			ReadVertices(pObject, m_CurrentChunk);
			break;

		case OBJECT_FACES:				 
			ReadVertexIndices(pObject, m_CurrentChunk);
			break;

		case OBJECT_MATERIAL:			 
			ReadObjectMaterial(pModel, pObject, m_CurrentChunk);			
			break;
			
		case OBJECT_UV:					
			ReadUVCoordinates(pObject, m_CurrentChunk);
			break;
			
		default: 
			m_CurrentChunk->bytesRead += (unsigned int)fread(ibuffer, 1, m_CurrentChunk->length - m_CurrentChunk->bytesRead, m_FilePointer);
			break;
		}
		
		pPreviousChunk->bytesRead += m_CurrentChunk->bytesRead;
	}
	
	delete m_CurrentChunk;
	m_CurrentChunk = pPreviousChunk;
}


void Parser3DS::ProcessNextMaterialChunk(Model3D *pModel, Chunk *pPreviousChunk)
{
	m_CurrentChunk = new Chunk;
	
	while (pPreviousChunk->bytesRead < pPreviousChunk->length)
	{
		ReadChunk(m_CurrentChunk);
		
		//switch (pModel, m_CurrentChunk->ID)
		switch (m_CurrentChunk->ID)
		{
		case MATNAME:						 
			m_CurrentChunk->bytesRead += (unsigned int)fread(pModel->pMaterials[pModel->numOfMaterials - 1].strName, 1, m_CurrentChunk->length - m_CurrentChunk->bytesRead, m_FilePointer);
			break;
			
		case MATDIFFUSE:					 
			ReadColorChunk(&(pModel->pMaterials[pModel->numOfMaterials - 1]), m_CurrentChunk);
			break;
			
		case MATMAP:						 
			ProcessNextMaterialChunk(pModel, m_CurrentChunk);
			break;
			
		case MATMAPFILE:						 
			m_CurrentChunk->bytesRead += (unsigned int)fread(pModel->pMaterials[pModel->numOfMaterials - 1].strFile, 1, m_CurrentChunk->length - m_CurrentChunk->bytesRead, m_FilePointer);
			break;
			
		default:  
			m_CurrentChunk->bytesRead += (unsigned int)fread(ibuffer, 1, m_CurrentChunk->length - m_CurrentChunk->bytesRead, m_FilePointer);
			break;
		}
		pPreviousChunk->bytesRead += m_CurrentChunk->bytesRead;
	}
	
	delete m_CurrentChunk;
	m_CurrentChunk = pPreviousChunk;
}


void Parser3DS::ReadChunk(Chunk *pChunk)
{
	pChunk->bytesRead = (unsigned int)fread(&pChunk->ID, 1, 2, m_FilePointer);
	pChunk->bytesRead += (unsigned int)fread(&pChunk->length, 1, 4, m_FilePointer);
}


int Parser3DS::GetString(char *pBuffer)
{
	int index = 0;
	fread(pBuffer, 1, 1, m_FilePointer);
	while (*(pBuffer + index++) != 0) 
	{
	 fread(pBuffer + index, 1, 1, m_FilePointer);
	}
	return (int)strlen(pBuffer)+1;
}


void Parser3DS::ReadColorChunk(MaterialInfo *pMaterial, Chunk *pChunk)
{
	ReadChunk(m_TempChunk);
	m_TempChunk->bytesRead += (unsigned int)fread(pMaterial->color, 1, m_TempChunk->length - m_TempChunk->bytesRead, m_FilePointer);
	pChunk->bytesRead += m_TempChunk->bytesRead;
}


void Parser3DS::ReadVertexIndices(Object3D *pObject, Chunk *pPreviousChunk)
{
	unsigned short index = 0;					 
	pPreviousChunk->bytesRead += (unsigned int)fread(&pObject->numOfFaces, 1, 2, m_FilePointer);
	pObject->pFaces = new Face [pObject->numOfFaces];
	memset(pObject->pFaces, 0, sizeof(Face) * pObject->numOfFaces);
	for(int i = 0; i < pObject->numOfFaces; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			pPreviousChunk->bytesRead += (unsigned int)fread(&index, 1, sizeof(index), m_FilePointer);
			if(j < 3)
			{
			 pObject->pFaces[i].vertIndex[j] = index;
			}
		}
	}
}



void Parser3DS::ReadUVCoordinates(Object3D *pObject, Chunk *pPreviousChunk)
{
	pPreviousChunk->bytesRead += (unsigned int)fread(&pObject->numTexVertex, 1, 2, m_FilePointer);
	pObject->pTexVerts = new Vector2[pObject->numTexVertex];
	pPreviousChunk->bytesRead += (unsigned int)fread(pObject->pTexVerts, 1, pPreviousChunk->length - pPreviousChunk->bytesRead, m_FilePointer);
}


void Parser3DS::ReadVertices(Object3D *pObject, Chunk *pPreviousChunk)
{ 
	pPreviousChunk->bytesRead += (unsigned int)fread(&(pObject->numOfVerts), 1, 2, m_FilePointer);
	pObject->pVerts = new Vector3[pObject->numOfVerts];
	memset(pObject->pVerts, 0, sizeof(Vector3) * pObject->numOfVerts);
	pPreviousChunk->bytesRead += (unsigned int)fread(pObject->pVerts, 1, pPreviousChunk->length - pPreviousChunk->bytesRead, m_FilePointer);
}


void Parser3DS::ReadObjectMaterial(Model3D *pModel, Object3D *pObject, Chunk *pPreviousChunk)
{
	char strMaterial[255] = {0};			
	pPreviousChunk->bytesRead += GetString(strMaterial);
	for(int i = 0; i < pModel->numOfMaterials; i++)
	{
		if(strcmp(strMaterial, pModel->pMaterials[i].strName) == 0)
		{
			pObject->materialID = i;
			if(strlen(pModel->pMaterials[i].strFile) > 0) {
				pObject->bHasTexture = true;
			}	
			break;
		}
	}
	pPreviousChunk->bytesRead += (unsigned int)fread(ibuffer, 1, pPreviousChunk->length - pPreviousChunk->bytesRead, m_FilePointer);
}			

void Parser3DS::ComputeNormals(Model3D *pModel)
{
	Vector3 vVector1, vVector2, vNormal, vPoly[3];
	if(pModel->numOfObjects <= 0)
		return;
	printf("Number of objects: %d\n", pModel->numOfObjects);
	for(int index = 0; index < pModel->numOfObjects; index++)
	{
		printf("Object no. #%d\n", index);
		Object3D *pObject = &(pModel->pObject[index]);
		
		Vector3 *pNormals		= new Vector3 [pObject->numOfFaces];
		Vector3 *pTempNormals	= new Vector3 [pObject->numOfFaces];
		pObject->pNormals		= new Vector3 [pObject->numOfVerts];
		printf("Number of vertices: %d\n", pObject->numOfVerts);
		printf("Number of faces: %d\n", pObject->numOfFaces);
		for(int i=0; i < pObject->numOfFaces; i++)
		{												
			//printf("Face no. #%d/%d\n", i, pObject->numOfFaces);
			vPoly[0] = pObject->pVerts[pObject->pFaces[i].vertIndex[0]];
			vPoly[1] = pObject->pVerts[pObject->pFaces[i].vertIndex[1]];
			vPoly[2] = pObject->pVerts[pObject->pFaces[i].vertIndex[2]];
			vVector1 = vPoly[0] - vPoly[2];				
			vVector2 = vPoly[2] - vPoly[1];				
			vNormal  = Cross(vVector1, vVector2);		
			pTempNormals[i] = vNormal;					
			vNormal  = Normalize(vNormal);				
			pNormals[i] = vNormal;						
		}
		printf("Faces done.\nInterpolating.\n");
		Vector3 vSum(0.0, 0.0, 0.0);
		Vector3 vZero = vSum;
		int shared=0;
		//printf("Number of vertices: %d\n", pObject->numOfVerts);
		for (int i = 0; i < pObject->numOfVerts; i++)			
		{
			//printf("Vertex no. #%d/%d\n", i, pObject->numOfVerts);
			for (int j = 0; j < pObject->numOfFaces; j++)	
			{
				//printf("Face no. #%d\n", j);
				if (pObject->pFaces[j].vertIndex[0] == i || 
					pObject->pFaces[j].vertIndex[1] == i || 
					pObject->pFaces[j].vertIndex[2] == i)
				{
					vSum = vSum + pTempNormals[j];			
					shared++;								
				}
			}      
			
			pObject->pNormals[i] = DivideVectorByScaler(vSum, float(-shared));
			pObject->pNormals[i] = Normalize(pObject->pNormals[i]);	
			vSum = vZero;									 
			shared = 0;										
		}
		delete [] pTempNormals;
		delete [] pNormals;
	}
 printf("Interpolated.\n");
}


