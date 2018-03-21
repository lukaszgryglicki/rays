all: tester
tester: tester.o i386_asm_rays.o Makefile
	gcc -g -o tester tester.o i386_asm_rays.o -lc -lm
tester.o: tester.c Makefile
	gcc -O3 --save-temps -ffast-math -g -ffast-math -D_REENTRANT -pthread -Wall  -c tester.c
i386_asm_rays.o: i386_asm_rays.s Makefile
	as -o i386_asm_rays.o i386_asm_rays.s
clean:
	-rm -f tester *.o *.i
