#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

readonly USAGE="usage: $0 [-t TAG] [-u USER_ID] [-b BASE] [-v VERSION]"
readonly HELP="$USAGE

Build the Docker image for cal-disp. Run from the repository root.

options:
  -t, --tag TAG        Docker image name/tag (default: opera-adt/cal-disp:latest)
  -u, --user-id ID     User ID for docker image (default: 1000)
  -b, --base BASE      Base image (default: ubuntu:22.04)
  -v, --version VER    PEP 440 version baked into the image (default: from
                       git describe, see docker/pep440-version.sh)
  -h, --help           Show this help and exit
"

# Defaults
tag="opera-adt/cal-disp:latest"
base=""
user_id=""
version=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--tag)
            tag="$2"
            shift 2
            ;;
        -u|--user-id)
            user_id="$2"
            shift 2
            ;;
        -b|--base)
            base="$2"
            shift 2
            ;;
        -v|--version)
            version="$2"
            shift 2
            ;;
        -h|--help)
            echo "$HELP"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "$USAGE"
            exit 1
            ;;
    esac
done

if [[ ! -f docker/Dockerfile ]]; then
    echo "Run this script from the repository root (docker/Dockerfile not found)" >&2
    exit 1
fi

# The image has no .git, so the version setuptools_scm would derive from
# it is computed here and passed in (otherwise every image reports the
# same placeholder version in product metadata).
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -z "$version" ]]; then
    version="$("$script_dir/pep440-version.sh")"
fi
echo "cal-disp version: $version"

# Build docker command
build_args=(
    "docker" "build"
    "--network=host"
    "--tag" "$tag"
    "--file" "docker/Dockerfile"
    "--build-arg" "CAL_DISP_VERSION=$version"
)

[[ -n "$base" ]] && build_args+=("--build-arg" "BASE=$base")
[[ -n "$user_id" ]] && build_args+=("--build-arg" "MAMBA_USER_ID=$user_id")

build_args+=(".")

# Execute build
echo "${build_args[@]}"
"${build_args[@]}"

# Usage examples
cat << EOF

To run the image:
  docker run --rm -it $tag cal-disp --help

To run on a PGE runconfig (the working directory in the image is /home/work):
  docker run --user \$(id -u):\$(id -g) -v \$PWD:/home/work --rm -it $tag cal-disp run runconfig.yaml
EOF

# where...
#     --user $(id -u):$(id -g)  # Needed to avoid permission issues when writing to the mounted volume.
#     -v $PWD:/home/work  # Mounts the current directory to /home/work, the WORKDIR of the image.
#     --rm  # Removes the container after it exits.
#     -it  # Needed to keep the container running after the command exits.
#     opera-adt/cal-disp:latest  # The name of the image to run.
#     cal-disp run runconfig.yaml # The `cal-disp` command that is run in the container
