gcc -O3 -o cyg_md22dat.exe md22dat.c -lopengl32 -lglut32 -lglu32
strip -s cyg_*.exe
rm ogromod_*.dat

