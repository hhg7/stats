#!/bin/sh
set -e
perl md2pod.pl                          # regenerate README/POD
git add -A && git commit -m "Update generated docs" || true
dzil clean
dzil build
echo "==== tarball contents (verify: no .c/.o/.dll/.bs/.gcda/blib) ===="
tar tzf Stats-LikeR-*.tar.gz
echo "If that looks clean, run: dzil release"
