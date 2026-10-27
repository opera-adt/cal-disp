#!/usr/bin/env bash

# Enable common error handling options.
set -o errexit
set -o nounset
set -o pipefail

# shellcheck disable=SC2016  # the backticks in the help text are literal
readonly HELP='usage: ./create-lockfile.sh --file ENVFILE [--pkgs PACKAGE ...] [--no-docker] > specfile.txt

Create a conda lockfile (an @EXPLICIT spec with md5 checksums) from an
environment YAML file, for reproducible environments (docker/Dockerfile
installs docker/conda-lock.txt with `micromamba create -f`).

The lock is for linux-64. By default the solve runs in a clean micromamba
Docker container; with --no-docker it runs with the local mamba/micromamba/
conda instead (which must then be on a linux-64 machine).

options:
--file ENVFILE    Specify a YAML file containing package specifications.
--pkgs PACKAGE    Specify additional packages separated by spaces. Example: --pkgs numpy scipy
--no-docker       Solve locally with mamba (or micromamba/conda) instead of Docker.
-h, --help        Show this help message and exit

example (from the repository root):
  docker/create-lockfile.sh --file environment.yml --no-docker > docker/conda-lock.txt
'

# The lockfile captures only the conda layer. Strip any pip editable
# install of the local project (e.g. a `- pip:` block with `- -e .`)
# before locking: the solve runs in a container / scratch prefix that has
# no project source, so pip would fail with "does not appear to be a Python
# project". The cal-disp package itself is installed separately (see
# docker/Dockerfile), so environment.yml can keep `-e .` for local dev
# without breaking locking.
sanitize_envfile() {
    local ENVFILE="$1"
    # Global (not local) so the EXIT trap can still see it after this function
    # returns. A RETURN trap would re-fire on main's return where it is unset.
    SANITIZED=$(mktemp --suffix=.yml)
    trap 'rm -f "${SANITIZED:-}"; rm -rf "${LOCAL_PREFIX:-}"' EXIT
    grep -vE '^[[:space:]]*(- pip:[[:space:]]*|- -e[[:space:]].*)$' "$ENVFILE" > "$SANITIZED"
    # mktemp creates mode 0600; the container user (mambauser) must be able to
    # read the bind-mounted file, otherwise micromamba reports "bad file".
    chmod 0644 "$SANITIZED"
}

# Print the sorted lockfile: the 4 header lines are kept in place, the
# package URLs after them are sorted alphabetically.
sort_pkglist() {
    # `conda list --explicit` adds a 5th "# created-by" header line; drop it
    # so the file has the same header as `micromamba env export --explicit`.
    grep -v '^# created-by' | (
        sed -u 4q
        sort
    )
}

install_packages_docker() {
    local ENVFILE
    ENVFILE=$(realpath "$1")
    shift

    sanitize_envfile "$ENVFILE"

    # Prepare arguments for the command. The extra packages are joined into
    # one string for the container's shell (an unindexed array expansion
    # passed only the first package).
    local FILE_ARG="--file /tmp/environment.yml"
    local PKGS_ARGS="$*"

    # Get concretized package list.
    local PKGLIST
    PKGLIST=$(docker run --rm --network=host \
        -v "$SANITIZED:/tmp/environment.yml:ro" \
        mambaorg/micromamba:1.1.0 bash -c "\
            micromamba install -y -n base $FILE_ARG $PKGS_ARGS > /dev/null && \
            micromamba env export --explicit")

    echo "$PKGLIST" | sort_pkglist
}

install_packages_local() {
    local ENVFILE
    ENVFILE=$(realpath "$1")
    shift
    local PACKAGES=("$@")

    if [[ "$(uname -s)-$(uname -m)" != "Linux-x86_64" ]]; then
        echo "--no-docker needs a linux-64 machine (the lock is for linux-64)" >&2
        exit 1
    fi

    local SOLVER
    for SOLVER in mamba micromamba conda; do
        command -v "$SOLVER" > /dev/null 2>&1 && break
        SOLVER=""
    done
    if [[ -z "$SOLVER" ]]; then
        echo "--no-docker needs mamba, micromamba or conda on PATH" >&2
        exit 1
    fi

    sanitize_envfile "$ENVFILE"

    # Solve into a scratch prefix; TMPDIR is honoured (mktemp uses it).
    LOCAL_PREFIX=$(mktemp -d --suffix=.cal-disp-lock)
    rm -rf "$LOCAL_PREFIX"
    echo "Solving $ENVFILE with $SOLVER into $LOCAL_PREFIX ..." >&2
    "$SOLVER" env create -y -p "$LOCAL_PREFIX" --file "$SANITIZED" > /dev/null
    if [[ "${#PACKAGES[@]}" -gt 0 ]]; then
        "$SOLVER" install -y -p "$LOCAL_PREFIX" -c conda-forge "${PACKAGES[@]}" > /dev/null
    fi

    # Get concretized package list (micromamba has no --md5; conda/mamba do).
    if [[ "$SOLVER" == "micromamba" ]]; then
        "$SOLVER" env export -p "$LOCAL_PREFIX" --explicit | sort_pkglist
    else
        "$SOLVER" list -p "$LOCAL_PREFIX" --explicit --md5 | sort_pkglist
    fi
}

main() {
    local ENVFILE=""
    local PACKAGES=()
    local USE_DOCKER=1

    while [[ "$#" -gt 0 ]]; do
        case $1 in
        --file)
            shift
            if [[ -z "${1-}" ]]; then
                echo "No file provided after --file" >&2
                exit 1
            fi
            ENVFILE="$1"
            shift
            ;;
        --pkgs)
            shift
            while [[ "$#" -gt 0 && ! "$1" =~ ^-- ]]; do
                PACKAGES+=("$1")
                shift
            done
            ;;
        --no-docker)
            USE_DOCKER=0
            shift
            ;;
        -h | --help)
            echo "$HELP"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "$HELP"
            exit 1
            ;;
        esac
    done

    if [[ -z "$ENVFILE" ]]; then
        echo 'No environment file provided' >&2
        echo "$HELP"
        exit 1
    fi

    if [[ "$USE_DOCKER" -eq 0 ]]; then
        # ${arr[@]+"${arr[@]}"}: an empty array is not "unbound" (bash 3.2)
        install_packages_local "$ENVFILE" ${PACKAGES[@]+"${PACKAGES[@]}"}
    else
        # If no packages were passed, only the environment file is installed.
        install_packages_docker "$ENVFILE" ${PACKAGES[@]+"${PACKAGES[@]}"}
    fi
}

main "$@"
