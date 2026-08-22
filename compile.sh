#git rm -r --cached blib/
#git rm --cached *.o *.dll LikeR.c Makefile MYMETA.*
# OPTIMIZE *replaces* perl's own $Config{optimize} (-O2 here) rather than adding
# to it, so passing bare -Wall built the whole module at -O0.  That cost 1.3x on
# sum/mean/min/max and 2x on cor: 1e6 NVs summed in 4.10ms at -O0 against 2.63ms
# at -O2, measured on perl-5.44.0.  Keep -O2 ahead of the warning flags.
make clean && perl Makefile.PL OPTIMIZE='-O2 -Wall' && make && make test && make install
