#!/bin/sh
# Commit, build and push the version in progress, then stop short of uploading:
# `dzil release` is left to be run by hand once the tarball has been checked.
set -e
# Every path below is relative to the repo root, so run from there whatever the
# caller's directory.
cd "$(dirname "$0")"

perl md2pod.pl

# The version in progress, used as the commit message, as the history already does.
VERSION=$(perl -ne 'if (/\$VERSION\s*=\s*[\x27"]([^\x27"]+)/) { print $1; exit }' lib/Stats/LikeR.pm)
if [ -z "$VERSION" ]; then
	echo "dzil.sh: could not read \$VERSION from lib/Stats/LikeR.pm" >&2
	exit 1
fi

# md2pod.pl only checks that Changes parses. A newest section that names some
# other version still parses, and would ship as this release's notes.
CHANGES_VERSION=$(perl -ne 'if (/^(\d\S*)\s/) { print $1; exit }' Changes)
if [ "$CHANGES_VERSION" != "$VERSION" ]; then
	echo "dzil.sh: \$VERSION is $VERSION but the newest section of Changes is '$CHANGES_VERSION'" >&2
	exit 1
fi

# [Git::GatherDir] builds the dist only from files git tracks, so a new test left
# untracked is silently missing from the tarball. Add the shipped paths by name
# rather than `git add *`, and exclude *.orig: PruneFiles keeps anything under
# lib/ and t/, so a committed backup there would ship.
git add -- lib t LikeR.xs ppport.h Changes README.md CONTRIBUTING.md SECURITY.md ':!*.orig'
# -a also picks up modified and deleted tracked files outside those paths.
# Commit only when something changed, rather than `git commit || true`, which
# also swallowed a commit that failed for a real reason.
if ! git diff-index --quiet HEAD --; then
	git commit -am "$VERSION"
fi

# `dzil clean` removes all of .build/, including the logs and work trees that
# test.all.perls.pl keeps in .build/multiperl. That is intentional: the
# maintainer wants a release to start with .build/ empty.
dzil clean
dzil build

TARBALL="Stats-LikeR-$VERSION.tar.gz"
if [ ! -f "$TARBALL" ]; then
	echo "dzil.sh: dzil build did not write $TARBALL" >&2
	exit 1
fi
# Build products, coverage output and editor/patch leftovers must never ship.
# The top-level directory (listed with no trailing slash) is dropped and its
# prefix stripped first, so a Stats-LikeR-* left is a nested one: a dzil build
# that gathered an earlier build into this one.
BAD=$(tar tzf "$TARBALL" | sed -e "/^Stats-LikeR-$VERSION\$/d" -e "s|^Stats-LikeR-$VERSION/||" | grep -E \
	'\.(c|o|obj|so|dll|dylib|bs|def|exp|gcda|gcno|gcov|orig|swp)$|(^|/)(blib|cover_db|\.build|Stats-LikeR-[^/]*)(/|$)|^(Makefile|Makefile\.old|MYMETA\.(json|yml)|pm_to_blib)$' \
	|| true)	# grep exits 1 when nothing matches, which is the good case
if [ -n "$BAD" ]; then
	echo "dzil.sh: $TARBALL contains files that must not ship:" >&2
	echo "$BAD" >&2
	exit 1
fi

# Pushed only once the build and the tarball check have both passed.
git push
echo "==== $TARBALL contents ===="
tar tzf "$TARBALL"
echo "If that looks right, run: dzil release"
