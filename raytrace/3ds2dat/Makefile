CFLAGS= -Wall -ansi -pedantic -pthread -O3 -DUNIX
LDFLAGS= -L/usr/X11R6/lib -L/usr/local/lib -lglut -ljpeg
all: 3dsLoader.bsd
3dsLoader.bsd: 3dsLoader.o 3dsParser.o JpegLoader.o ImageIO.o Main.o Vector.o libjpeg_freebsd-5.3-i386.a
	#g++ $(CFLAGS) -o 3dsLoader.bsd 3dsLoader.o 3dsParser.o JpegLoader.o Main.o Vector.o ImageIO.o /usr/local/lib/libefence.a $(LDFLAGS)
	g++ $(CFLAGS) -o 3dsLoader.bsd 3dsLoader.o 3dsParser.o JpegLoader.o Main.o Vector.o ImageIO.o $(LDFLAGS)
JpegLoader.o: JpegLoader.cpp JpegLoader.h
	g++ $(CFLAGS) -c JpegLoader.cpp -I/usr/local/include
Main.o: Main.cpp Main.h
	g++ $(CFLAGS) -c Main.cpp
Vector.o: Vector.cpp Vector.h
	g++ $(CFLAGS) -c Vector.cpp
ImageIO.o: ImageIO.cpp ImageIO.h
	g++ $(CFLAGS) -c ImageIO.cpp -I/usr/local/include
3dsParser.o: 3dsParser.cpp 3dsParser.h
	g++ $(CFLAGS) -c 3dsParser.cpp
3dsLoader.o: 3dsLoader.cpp 3dsLoader.h
	g++ $(CFLAGS) -c 3dsLoader.cpp -I/usr/X11R6/include
clean:
	-rm -f *.o 3dsLoader.bsd

