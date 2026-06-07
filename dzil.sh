#!/bin/sh
set -e
perl md2pod.pl
git commit -am "Update generated docs" || true   # -a = tracked files only; NOT -A
dzil clean
dzil build
echo "==== tarball contents (verify: no .c/.o/.dll/.bs/.gcda/blib, no Stats-LikeR/) ===="
tar tzf Stats-LikeR-*.tar.gz
echo "If that looks clean, run: dzil release"
