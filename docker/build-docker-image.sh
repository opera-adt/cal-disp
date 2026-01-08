#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

readonly USAGE="usage: $0 [-t TAG] [-u USER_ID] [-b BASE]"
readonly HELP="$USAGE

Build the docker image for opera_tropo.

options:
  -t, --tag TAG        Docker image name/tag (default: opera-adt/tropo:latest)
  -u, --user-id ID     User ID for docker image (default: 1000)
  -b, --base BASE      Base image (default: ubuntu:22.04)
  -h, --help           Show this help and exit
"

# Defaults
tag="opera-adt/opera_adt:latest"
base=""
user_id=""

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

# Build docker command
build_args=(
    "docker" "build"
    "--network=host"
    "--tag" "$tag"
    "--file" "docker/Dockerfile"
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

To run on a PGE runconfig:
  docker run --user \$(id -u):\$(id -g) -v \$PWD:/work --rm -it $tag cal-disp run runconfig.yaml
EOF

# where...
#     --user $(id -u):$(id -g)  # Needed to avoid permission issues when writing to the mounted volume.
#     -v $PWD:/work  # Mounts the current directory to the /work directory in the container.
#     --rm  # Removes the container after it exits.
#     -it  # Needed to keep the container running after the command exits.
#     opera-adt/disp-s1:latest  # The name of the image to run.
#     disp-s1 runconfig.yaml # The `disp-s1` command that is run in the container
