CFLAGS= -Wall -ansi -pedantic -DUNIX
all: 3ds2tri
3ds2tri: 3dsLoader.o 3dsParser.o Main.o Vector.o
	g++ $(CFLAGS) -o 3ds2tri 3dsLoader.o 3dsParser.o Main.o Vector.o
	strip -s 3ds2tri
Main.o: Main.cpp Main.h
	g++ $(CFLAGS) -c Main.cpp
Vector.o: Vector.cpp Vector.h
	g++ $(CFLAGS) -c Vector.cpp
3dsParser.o: 3dsParser.cpp 3dsParser.h
	g++ $(CFLAGS) -c 3dsParser.cpp
3dsLoader.o: 3dsLoader.cpp 3dsLoader.h
	g++ $(CFLAGS) -c 3dsLoader.cpp
clean:
	-rm -f *.o 3ds2tri

