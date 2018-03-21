#include "Main.h"
#include "Vector.h"
#include "3dsParser.h"
#include "3dsLoader.h"


int main(int argc, char** argv)
{
 bool loaded;
 Loader3D loader;
 printf("Create loader.\n");
 if (argc <= 1) loaded = loader.Load3DS("model.3ds");
 else loaded = loader.Load3DS(argv[1]);
 if (!loaded) 
   { 
	printf("Cannot load model.\n"); 
	return 1; 
   }
 loader.WriteTris(false);
 loader.Unload3DS();
 printf("Loader destroyed.\n");
 return 0;
}

