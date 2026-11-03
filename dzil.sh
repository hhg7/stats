#!/bin/sh
set -e
perl md2pod.pl

# The version in progress, used as the commit message, as the history already does.
VERSION=$(perl -ne 'if (/\$VERSION\s*=\s*[\x27"]([^\x27"]+)/) { print $1; exit }' lib/Stats/LikeR.pm)
if [ -z "$VERSION" ]; then
	echo "dzil.sh: could not read \$VERSION from lib/Stats/LikeR.pm" >&2
	exit 1
fi

# [Git::GatherDir] builds the dist only from files git tracks, so a new test left
# untracked is silently missing from the tarball. Add the shipped paths by name
# rather than `git add *`, and exclude *.orig: PruneFiles keeps anything under
# lib/ and t/, so a committed backup there would ship.
git add -- lib t LikeR.xs ppport.h Changes README.md CONTRIBUTING.md SECURITY.md ':!*.orig'
# -a also picks up modified and deleted tracked files outside those paths.
git commit -am "$VERSION" || true	# nothing to commit is not an error

dzil clean
dzil build
# Pushed only after a successful build: set -e stops the script before this on failure.
git push
echo "==== tarball contents (verify: no .c/.o/.dll/.bs/.gcda/blib, no Stats-LikeR/) ===="
tar tzf Stats-LikeR-*.tar.gz
echo "If that looks clean, run: dzil release"
