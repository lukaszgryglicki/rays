#ifndef _3DS_LOADER_H__MDBMA_
#define _3DS_LOADER_H__MDBMA_
#include "Main.h"
#include "3dsParser.h"
//#include "JpegLoader.h"

class Loader3D
{
	private:
		Parser3DS mLoad3ds;						
		unsigned int TextureArray3ds[MAX_TEXTURE];
		unsigned int def_tex;
		Model3D m3DModel;	
		char inputFileName[255];
		//void GLUT_PreInit();
		bool Texture_3ds(unsigned int textureArray[], char* strFileName, int ID);

	public:
		bool Load3DS(char* filename);
		//void GLUT_Render3DS();
		//void CallbackRender3DS();
		void Unload3DS();
		void WriteTris(bool);
};


#endif
