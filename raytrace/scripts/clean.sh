find . -name "*~" -exec rm -f "{}" \;
find . -name "*.[cChH]" -exec chmod -x "{}" \;
find . -name "*.[bB][mM][pP]" -exec chmod -x "{}" \;
find . -name "*.[jJ][pP]*[gG]" -exec chmod -x "{}" \;
find . -name "*.[dD][aA][tT]" -exec chmod -x "{}" \;
find . -name "*.[pP][cC][xX]" -exec chmod -x "{}" \;
find . -name "*.[mM][dD]2" -exec chmod -x "{}" \;
find . -name "Makefile" -exec chmod -x "{}" \;
find . -name "*.btree" -exec chmod -x "{}" \;
find . -name "*.[bB][iI][nN]" -exec chmod -x "{}" \;
find . -name "*.sizes" -exec chmod -x "{}" \;
find . -name "*.info" -exec chmod -x "{}" \;
find . -name "*.uli" -exec chmod -x "{}" \;
find . -name "*.[tT][xX][tT]" -exec chmod -x "{}" \;
find . -name "*.stackdump" -exec rm -f "{}" \;
find . -name "*.core" -exec rm -f "{}" \;
make clean
chmod -x TODO DONE README.TXT tags
rm -f debug.out *.gcov *.efence out *.i *.o

