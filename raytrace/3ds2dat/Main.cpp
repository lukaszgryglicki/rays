#include "Main.h"
#include "Vector.h"
#include "3dsParser.h"
#include "3dsLoader.h"

Loader3D loader;

int main(int argc, char** argv)
{
 bool loaded;
 //printf("Create loader.\n");
 if (argc <= 1) loaded = loader.Load3DS("model.3ds", argc, argv);
 else loaded = loader.Load3DS(argv[1], argc, argv);
 if (!loaded) 
   { 
	printf("Cannot load model.\n"); 
	return 1; 
   }
 loader.GLUT_Render3DS();
 loader.Unload3DS();
 //printf("Loader destroyed.\n");
 return 0;
}

