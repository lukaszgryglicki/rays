CC = gcc
all: unix # windows
unix: rays.debug rays96.debug rays32.debug rays32.fast rays64.fast rays96.fast rays.small torrusgen nurbs ball cube table convert wrapper 60faces rtriangle uli2dat cone sizes cylinder md22dat nurbs2dat randnurb randnurbfull snum getbmp tex terminal  iges2dat btreeconv
windows: cyg_rays.exe cyg_torrusgen.exe cyg_ball.exe cyg_cube.exe cyg_table.exe cyg_convert.exe cyg_60faces.exe cyg_wrapper.exe cyg_rtriangle.exe cyg_uli2dat.exe cyg_cone.exe cyg_cylinder.exe cyg_tex.exe cyg_nurbs2dat.exe cyg_randnurb.exe cyg_randnurbfull.exe cyg_snum.exe cyg_getbmp.exe cyg_terminal.exe cyg_sizes.exe cyg_md22dat.exe cyg_nurbs.exe cyg_iges2dat.exe cyg_btreeconv.exe
rays.efence: rays.c rayslib.c rayslib.h
	$(CC) -ffast-math  -Werror -D_REENTRANT -I/usr/local/include -L/usr/local/lib -DDEBUG -DNOGL -DNOINET -Wall -ansi -pedantic -g3 -o rays.efence rays.c rayslib.c /usr/lib/libefence.a  -ljpeg
rays32.debug: rays.c rayslib.c rayslib.h
	$(CC) -DREAL=__REAL_32bit -ffast-math  -Werror -D_REENTRANT -pthread -I/usr/X11R6/include -L/usr/X11R6/lib -I/usr/local/include -L/usr/local/lib -DDEBUG -Wall -ansi -pedantic -g3 -o rays32.debug rays.c rayslib.c -lglut -lm -ljpeg
rays.debug: rays.c rayslib.c rayslib.h
	$(CC) -DREAL=__REAL_64bit -ffast-math  -Werror -D_REENTRANT -pthread -I/usr/X11R6/include -L/usr/X11R6/lib -I/usr/local/include -L/usr/local/lib -DDEBUG -Wall -ansi -pedantic -g3 -o rays64.debug rays.c rayslib.c -lglut -lm -ljpeg
	ln -s rays64.debug rays.debug
rays96.debug: rays.c rayslib.c rayslib.h
	$(CC) -DREAL=__REAL_80bit -ffast-math  -Werror -D_REENTRANT -pthread -I/usr/X11R6/include -L/usr/X11R6/lib -I/usr/local/include -L/usr/local/lib -DDEBUG -Wall -ansi -pedantic -g3 -o rays96.debug rays.c rayslib.c -lglut -lm -ljpeg
rays.small: rays.c rayslib.c rayslib.h
	$(CC) -O3 -ffast-math -Werror -DNOTHREADS -DNOINET -DNOJPEG -DNOGL -DNOSIGNALS -Wall -ansi -pedantic -g -o rays.small rays.c rayslib.c -lc -lm
	strip -s rays.small
rays.static: rays.c rayslib.c rayslib.h
	$(CC) -funroll-all-loops -fomit-frame-pointer -fexpensive-optimizations -ffast-math -I/usr/local/include -L/usr/local/lib -pthread -O4 -ffast-math -Werror -static -DNOGL -Wall -ansi -pedantic -o rays.static rays.c rayslib.c  -ljpeg
	strip -s rays.static
rays32.fast: rays.c rayslib.c rayslib.h
	$(CC) -DREAL=__REAL_32bit -funroll-all-loops -fomit-frame-pointer -fexpensive-optimizations -ffast-math -Werror -D_REENTRANT -pthread -I/usr/X11R6/include -L/usr/X11R6/lib -I/usr/local/include -L/usr/local/lib -Wall -ansi -pedantic -O4 -o rays32.fast rays.c rayslib.c -lglut -lm -ljpeg
	strip -s rays32.fast
rays64.fast: rays.c rayslib.c rayslib.h
	$(CC) -DREAL=__REAL_64bit -funroll-all-loops -fomit-frame-pointer -fexpensive-optimizations -ffast-math -Werror -D_REENTRANT -pthread -I/usr/X11R6/include -L/usr/X11R6/lib -I/usr/local/include -L/usr/local/lib -Wall -ansi -pedantic -O4 -o rays64.fast rays.c rayslib.c -lglut -lm -ljpeg
	strip -s rays64.fast
	ln -s rays64.fast rays.fast
rays96.fast: rays.c rayslib.c rayslib.h
	$(CC) -DREAL=__REAL_80bit -funroll-all-loops -fomit-frame-pointer -fexpensive-optimizations -ffast-math -Werror -D_REENTRANT -pthread -I/usr/X11R6/include -L/usr/X11R6/lib -I/usr/local/include -L/usr/local/lib -Wall -ansi -pedantic -O4 -o rays96.fast rays.c rayslib.c -lglut -lm -ljpeg
	strip -s rays96.fast
cyg_rays.exe: rays.c rayslib.c rayslib.h
	gcc -D_REENTRANT -ffast-math -march=pentiumpro -O3 -Wall -ansi -pedantic -O3 -o cyg_rays.exe rays.c rayslib.c -lopengl32 -lglut32 -lglu32 -ljpeg -lc 
	strip -s cyg_rays.exe
torrusgen: torrusgen.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o torrusgen torrusgen.c -lc -lm
	strip -s torrusgen
iges2dat: iges2dat.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o iges2dat iges2dat.c -lc 
	strip -s iges2dat
btreeconv: btreeconv.c
	$(CC) -O3 -Wall -ansi -pedantic -g -o btreeconv btreeconv.c -lc 
	strip -s btreeconv
cyg_torrusgen.exe: torrusgen.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_torrusgen.exe torrusgen.c 
	strip -s cyg_torrusgen.exe
cyg_iges2dat.exe: iges2dat.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_iges2dat.exe iges2dat.c 
	strip -s cyg_iges2dat.exe
cyg_btreeconv.exe: btreeconv.c
	gcc -O3 -Wall -ansi -pedantic -g -o cyg_btreeconv.exe btreeconv.c 
	strip -s cyg_btreeconv.exe
convert: convert.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o convert convert.c -lc
	strip -s convert
cyg_convert.exe: convert.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_convert.exe convert.c
	strip -s cyg_convert.exe
rtriangle: rtriangle.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o rtriangle rtriangle.c -lc
	strip -s rtriangle
cyg_rtriangle.exe: rtriangle.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_rtriangle.exe rtriangle.c
	strip -s cyg_rtriangle.exe
cube: cube.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o cube cube.c -lc 
	strip -s cube
cyg_cube.exe: cube.c
	gcc -Wall -ansi -pedantic -g -o cyg_cube.exe cube.c 
	strip -s cyg_cube.exe
ball: ball.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o ball ball.c -lc 
	strip -s ball
cyg_ball.exe: ball.c
	gcc -Wall -ansi -pedantic -g -o cyg_ball.exe ball.c 
	strip -s cyg_ball.exe
60faces: 60faces.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o 60faces 60faces.c -lc 
	strip -s 60faces
cyg_60faces.exe: 60faces.c
	gcc -Wall -ansi -pedantic -g -o cyg_60faces.exe 60faces.c 
	strip -s cyg_60faces.exe
table: table.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o table table.c -lc 
	strip -s table
cyg_table.exe: table.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_table.exe table.c 
	strip -s cyg_table.exe
uli2dat: uli2dat.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o uli2dat uli2dat.c -lc 
	strip -s uli2dat
cyg_uli2dat.exe: uli2dat.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_uli2dat.exe uli2dat.c 
	strip -s cyg_uli2dat.exe
nurbs: nurbs.c libfparser.c libfparser.h
	$(CC) -O2 -pthread -I/usr/X11R6/include -L/usr/X11R6/lib -Wall -ansi -pedantic -g -o nurbs nurbs.c libfparser.c -lc  -lglut
	strip -s nurbs
md22dat: md22dat.c
	$(CC) -pthread -D_REENTRANT -I/usr/X11R6/include -L/usr/X11R6/lib -Wall -ansi -pedantic -O3 -o md22dat md22dat.c -lc  -lglut
	strip -s md22dat
cone: cone.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o cone cone.c -lc -lm
	strip -s cone
cyg_cone.exe: cone.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_cone.exe cone.c 
	strip -s cyg_cone.exe
cylinder: cylinder.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o cylinder cylinder.c -lc -lm
	strip -s cylinder
cyg_cylinder.exe: cylinder.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_cylinder.exe cylinder.c 
	strip -s cyg_cylinder.exe
wrapper: wrapper.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o wrapper wrapper.c -lc 
	strip -s wrapper
cyg_wrapper.exe: wrapper.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_wrapper.exe wrapper.c 
	strip -s cyg_wrapper.exe
sizes: sizes.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o sizes sizes.c -lc 
	strip -s sizes
getbmp: getbmp.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o getbmp getbmp.c -lc 
	strip -s getbmp
cyg_getbmp.exe: getbmp.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_getbmp.exe getbmp.c 
	strip -s cyg_getbmp.exe
tex: tex.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o tex tex.c -lc 
	strip -s tex
tex_tex.exe: tex.c
	gcc -O2 -Wall -ansi -pedantic -g -o tex_tex.exe tex.c 
	strip -s tex_tex.exe
cyg_md22dat.exe: md22dat.c
	gcc -O3 -o cyg_md22dat.exe md22dat.c -lopengl32 -lglut32 -lglu32
	strip -s cyg_md22dat.exe
cyg_sizes.exe: sizes.c
	gcc -O3 -o cyg_sizes.exe sizes.c
	strip -s cyg_sizes.exe
cyg_nurbs.exe: nurbs.c libfparser.c libfparser.h
	gcc -O3 -o cyg_nurbs.exe nurbs.c libfparser.c -lopengl32 -lglut32 -lglu32
	strip -s cyg_nurbs.exe
terminal: terminal.c
	$(CC) -O2 -Wall -ansi -pedantic -g -o terminal terminal.c -lc 
	strip -s terminal
cyg_terminal.exe: terminal.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_terminal.exe terminal.c 
	strip -s cyg_terminal.exe
randnurb: randnurb.c
	$(CC) -Wall -ansi -pedantic -O3 -o randnurb randnurb.c -lc 
	strip -s randnurb
cyg_randnurb.exe: randnurb.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_randnurb.exe randnurb.c 
	strip -s cyg_randnurb.exe
randnurbfull: randnurbfull.c
	$(CC) -Wall -ansi -pedantic -O3 -o randnurbfull randnurbfull.c -lc 
	strip -s randnurbfull
cyg_randnurbfull.exe: randnurbfull.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_randnurbfull.exe randnurbfull.c 
	strip -s cyg_randnurbfull.exe
snum: snum.c
	$(CC) -Wall -ansi -pedantic -O3 -o snum snum.c -lc 
	strip -s snum
cyg_snum.exe: snum.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_snum.exe snum.c 
	strip -s cyg_snum.exe
nurbs2dat: nurbs2dat.c
	$(CC) -Wall -ansi -pedantic -O3 -o nurbs2dat nurbs2dat.c -lc -lm
	strip -s nurbs2dat
cyg_nurbs2dat.exe: nurbs2dat.c
	gcc -O2 -Wall -ansi -pedantic -g -o cyg_nurbs2dat.exe nurbs2dat.c 
	strip -s cyg_nurbs2dat.exe
clean:
	-rm -f terminal getbmp snum randnurbfull randnurb nurbs2dat md22dat cylinder cone uli2dat rtriangle wrapper nurbs rays.debug rays.fast rays32.fast rays64.fast rays96.fast rays.small rays.gc* 60faces convert cube table torrusgen *.bin *.bmp cyg_*.exe *.jpeg debug.out sizes rays.static
