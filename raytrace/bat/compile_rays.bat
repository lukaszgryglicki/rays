gcc -DREAL=__REAL_32bit -D_REENTRANT -ffast-math -march=pentiumpro -O3 -Wall -ansi -pedantic -O3 -o cyg_rays32.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
gcc -DREAL=__REAL_64bit -D_REENTRANT -ffast-math -march=pentiumpro -O3 -Wall -ansi -pedantic -O3 -o cyg_rays64.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
gcc -DREAL=__REAL_80bit -D_REENTRANT -ffast-math -march=pentiumpro -O3 -Wall -ansi -pedantic -O3 -o cyg_rays96.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
strip -s cyg_rays32.exe
strip -s cyg_rays64.exe
strip -s cyg_rays96.exe

