git rm -r --cached blib/
git rm --cached *.o *.dll LikeR.c Makefile MYMETA.*
dzil clean && dzil build && dzil release
