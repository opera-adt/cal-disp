#!/usr/bin/env bash

# Enable common error handling options.
set -o errexit
set -o nounset
set -o pipefail

readonly HELP='usage: ./create-lockfile.sh --file ENVFILE [--pkgs PACKAGE ...] > specfile.txt

Create a conda lockfile from an environment YAML file and additional packages for reproducible environments.

options:
--file ENVFILE    Specify a YAML file containing package specifications.
--pkgs PACKAGE    Specify additional packages separated by spaces. Example: --pkgs numpy scipy
-h, --help        Show this help message and exit
'

install_packages() {
    local ENVFILE=$(realpath "$1")
    shift
    local PACKAGES="$@"

    # The lockfile captures only the conda layer. Strip any pip editable
    # install of the local project (e.g. a `- pip:` block with `- -e .`)
    # before locking: this runs in a container that has no project source
    # mounted, so pip would fail with "does not appear to be a Python project".
    # The cal-disp package itself is installed separately (see docker/Dockerfile),
    # so environment.yml can keep `-e .` for local dev without breaking locking.
    # Global (not local) so the EXIT trap can still see it after this function
    # returns. A RETURN trap would re-fire on main's return where it is unset.
    SANITIZED=$(mktemp --suffix=.yml)
    trap 'rm -f "${SANITIZED:-}"' EXIT
    grep -vE '^[[:space:]]*(- pip:[[:space:]]*|- -e[[:space:]].*)$' "$ENVFILE" > "$SANITIZED"
    # mktemp creates mode 0600; the container user (mambauser) must be able to
    # read the bind-mounted file, otherwise micromamba reports "bad file".
    chmod 0644 "$SANITIZED"

    # Prepare arguments for the command
    local FILE_ARG="--file /tmp/environment.yml"
    if [[ -n "$PACKAGES" ]]; then
        PKGS_ARGS=(${PACKAGES[@]})
    else
        PKGS_ARGS=""
    fi

    # Get concretized package list.
    local PKGLIST
    PKGLIST=$(docker run --rm --network=host \
        -v "$SANITIZED:/tmp/environment.yml:ro" \
        mambaorg/micromamba:1.1.0 bash -c "\
            micromamba install -y -n base $FILE_ARG $PKGS_ARGS > /dev/null && \
            micromamba env export --explicit")

    # Sort packages alphabetically.
    # (The first 4 lines are assumed to be header lines and ignored.)
    echo "$PKGLIST" | (
        sed -u 4q
        sort
    )
}

main() {
    local ENVFILE=""
    local PACKAGES=()

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

    # If no packages were passed, install only the packages in the environment file.
    if [[ "${#PACKAGES[@]}" -eq 0 ]]; then
        install_packages "$ENVFILE"
    else
        install_packages "$ENVFILE" "${PACKAGES[@]}"
    fi
}

main "$@"
