#!/usr/bin/env bash
# Print the PEP 440 version of the checkout, the way setuptools_scm would
# (version_scheme "no-guess-dev", see pyproject.toml), from `git describe`.
# The Docker build has no .git (see .dockerignore), so build-docker-image.sh
# passes this as --build-arg CAL_DISP_VERSION and the image reports the real
# version in product metadata instead of a hard-coded placeholder.
#
# usage: pep440-version.sh              # from `git describe --tags --always --dirty`
#        pep440-version.sh DESCRIBE     # convert a given describe string
#
#   v0.2                   -> 0.2
#   v0.2-27-g88c0c46       -> 0.2.post1.dev27+g88c0c46
#   v0.2-27-g88c0c46-dirty -> 0.2.post1.dev27+g88c0c46.d<today, UTC>
#   v0.2-dirty             -> 0.2+d<today, UTC>
#   88c0c46                -> 0.0.0.post1.dev0+g88c0c46   (no tag reachable)
#
# Anything else (e.g. a tag that is not a version) is an error; pass the
# version explicitly then: build-docker-image.sh --version X.Y.Z
set -o errexit
set -o nounset
set -o pipefail

pep440_from_describe() {
    local describe="$1"
    local dirty=""
    if [[ "$describe" == *-dirty ]]; then
        dirty=".d$(date -u +%Y%m%d)"
        describe="${describe%-dirty}"
    fi
    local tag distance hash
    if [[ "$describe" =~ ^[0-9a-f]{7,40}$ ]]; then
        # No tag reachable: `git describe --always` printed the bare hash
        echo "0.0.0.post1.dev0+g${describe}${dirty}"
    elif [[ "$describe" =~ ^v?([0-9][0-9A-Za-z.]*)-([0-9]+)-g([0-9a-f]+)$ ]]; then
        tag="${BASH_REMATCH[1]}"
        distance="${BASH_REMATCH[2]}"
        hash="${BASH_REMATCH[3]}"
        echo "${tag}.post1.dev${distance}+g${hash}${dirty}"
    elif [[ "$describe" =~ ^v?([0-9][0-9A-Za-z.]*)$ ]]; then
        tag="${BASH_REMATCH[1]}"
        # A clean tagged commit has no local part; a dirty one keeps the date
        echo "${tag}${dirty:++${dirty#.}}"
    else
        echo "cannot convert git describe output '$1' to a PEP 440 version" >&2
        return 1
    fi
}

if [[ "$#" -gt 0 ]]; then
    pep440_from_describe "$1"
else
    # --abbrev=9: the hash length setuptools_scm prints, so the image and a
    # `pip install` of the same clean checkout report the same version
    describe="$(git describe --tags --always --dirty --abbrev=9)"
    pep440_from_describe "$describe"
fi
