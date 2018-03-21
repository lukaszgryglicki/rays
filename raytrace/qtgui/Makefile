all: unix # windows
windows: cyg_qtmain.exe
unix: qtmain
cyg_qtmain.exe: qtmain.o qtgui.o ../rayslib.c
	g++ -o cyg_qtmain.exe qtmain.o qtgui.o ../rayslib.c -L/usr/lib/qt3/lib -lqt-mt -lopengl32 -lglu32 -lglut32 -ljpeg -lc -lm
	#g++ -o cyg_qtmain.exe rayslib.o qtmain.o qtgui.o -L/usr/lib/qt3/lib -lqt-mt -lopengl32 -lglu32 -lglut32 -ljpeg -lc -lm
qtmain: qtmain.o qtgui.o rayslib.o
	#g++ -pthread -o qtmain qtmain.o qtgui.o rayslib.o ../libjpeg-freebsd-6.0-amd64.a -L/usr/local/lib -L/usr/X11R6/lib -lqt-mt -lglut
	g++ -pthread -o qtmain qtmain.o qtgui.o rayslib.o ../libjpeg-freebsd-5.3-i386.a -L/usr/local/lib -L/usr/X11R6/lib -lqt-mt -lglut
qtmain.o: qtmain.cxx
	g++ -c qtmain.cxx -I/usr/X11R6/include -I/usr/include/qt3
qtgui.o: moc_qtgui.cxx
	g++ -c qtgui.cxx -I/usr/X11R6/include -I/usr/include/qt3
rayslib.o: ../rayslib.c
	g++ -I/usr/X11R6/include -I/usr/local/include -D_REENTRANT -Wall -ansi -pedantic -O3 -c ../rayslib.c -o rayslib.o
moc_qtgui.cxx: qtgui.hxx qtgui.cxx
	moc -o moc_qtgui.cxx qtgui.hxx	
clean:
	-rm -f moc_* *~ *.o cyg_*.exe qtmain
