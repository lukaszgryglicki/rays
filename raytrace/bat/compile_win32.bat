gcc -O3 -o cyg_md22dat.exe md22dat.c -lopengl32 -lglut32 -lglu32
gcc -O3 -o cyg_nurbs.exe nurbs.c libfparser.c -lopengl32 -lglut32 -lglu32
gcc -O2 -o cyg_torrusgen.exe torrusgen.c -lc -lm
gcc -O2 -o cyg_table.exe table.c -lc -lm
gcc -O2 -o cyg_cube.exe cube.c -lc -lm
gcc -O2 -o cyg_convert.exe convert.c -lc -lm
gcc -DREAL=__REAL_32bit -D_REENTRANT -ffast-math -march=pentiumpro -O3 -Wall -ansi -pedantic -O3 -o cyg_rays32.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
gcc -DREAL=__REAL_64bit -D_REENTRANT -ffast-math -march=pentiumpro -O3 -Wall -ansi -pedantic -O3 -o cyg_rays64.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
gcc -DREAL=__REAL_80bit -D_REENTRANT -ffast-math -march=pentiumpro -O3 -Wall -ansi -pedantic -O3 -o cyg_rays96.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc -lm
gcc -O2 -o cyg_wrapper.exe wrapper.c -lc
gcc -O2 -o cyg_60faces.exe 60faces.c -lc
gcc -O2 -o cyg_rtriangle.exe rtriangle.c -lc
gcc -O2 -o cyg_uli2dat.exe uli2dat.c -lc
gcc -O2 -o cyg_sizes.exe sizes.c -lc
gcc -O2 -o cyg_cone.exe cone.c -lc
gcc -O2 -o cyg_cylinder.exe cylinder.c -lc
gcc -O2 -o cyg_nurbs2dat.exe nurbs2dat.c -lc
gcc -O2 -o cyg_randnurb.exe randnurb.c -lc
gcc -O2 -o cyg_randnurbfull.exe randnurbfull.c -lc
gcc -O2 -o cyg_getbmp.exe getbmp.c -lc
gcc -O2 -o cyg_terminal.exe terminal.c -lc
gcc -O2 -o cyg_snum.exe snum.c -lc
strip -s cyg_*.exe

