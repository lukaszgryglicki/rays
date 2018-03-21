gcc -DREAL=__REAL_32bit -DDEBUG -D_REENTRANT -march=pentiumpro -g3 -Wall -ansi -pedantic -o cyg_raysd32.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
gcc -DREAL=__REAL_64bit -DDEBUG -D_REENTRANT -march=pentiumpro -g3 -Wall -ansi -pedantic -o cyg_raysd64.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
gcc -DREAL=__REAL_80bit -DDEBUG -D_REENTRANT -march=pentiumpro -g3 -Wall -ansi -pedantic -o cyg_raysd96.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm

