#ifndef _3DS_LOADER_H__MDBMA_
#define _3DS_LOADER_H__MDBMA_
#include "Main.h"
#include "3dsParser.h"
#include "JpegLoader.h"

struct DatInfo
{
 double datObserverX, datLightX;
 double datObserverY, datLightY;
 double datObserverZ, datLightZ;
 double datRotationX;
 double datRotationY;
 double datRotationZ;
 double datTranslationX;
 double datTranslationY;
 double datTranslationZ;
 double datScale;
 double datColorRed, datColorGreen, datColorBlue;
 double datSpecularRed, datSpecularGreen, datSpecularBlue;
 double datTransRed, datTransGreen, datTransBlue;
 double datTransFactorRed, datTransFactorGreen, datTransFactorBlue;
 double datSpecularFactor;
 int datScreenX, datScreenY, datScreenZ;
 int datBackup, datTexturesNum, datTrisNum;
 int datFaces;
 int datNoTex;
 char datTextureDir[255];
};

class Loader3D
{
	private:
		Parser3DS mLoad3ds;						
		unsigned int TextureArray3ds[MAX_TEXTURE];
		unsigned int TextureDatCopied[MAX_TEXTURE];
		unsigned int def_tex;
		Model3D m3DModel;	
		DatInfo dinfo;
		char inputFileName[255];
		int dargc;
		char **dargv;
		
		void GLUT_PreInit();
		bool Texture_3ds(unsigned int textureArray[], char* strFileName, int ID);
		void SetDatDefaults();
		void ParseDatOptions();
		void WriteDatHeader(FILE*);
		void TextureMapping(double&);

	public:
		bool Load3DS(char* filename, int lb, char** par);
		void GLUT_Render3DS();
		void CallbackRender3DS();
		void Unload3DS();
		void WriteTris(bool = false);
		void WriteDat();
		void SetDatTransformation();
};


#endif
