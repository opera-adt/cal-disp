# cal-disp Docker

Docker container for OPERA DISP-S1 calibration tools. All commands below are
run from the repository root.

## Quick Start
```bash
# Build
docker/build-docker-image.sh -t opera-adt/cal-disp:latest

# Run
docker run --rm -it opera-adt/cal-disp:latest cal-disp --help
docker run --rm -it opera-adt/cal-disp:latest cal-disp --version
```

## Building

### Basic Build
```bash
docker/build-docker-image.sh
```

Builds `opera-adt/cal-disp:latest` using default settings (Ubuntu 22.04, user
ID 1000). The script derives the package version from `git describe`
(`docker/pep440-version.sh`, e.g. `v0.2-27-g88c0c46` becomes
`0.2.post1.dev27+g88c0c46`, the same scheme setuptools_scm uses) and passes it
as `--build-arg CAL_DISP_VERSION=...`. The build context has no `.git`, so
without that argument the image reports the placeholder `0.0.0.dev0+unknown`
in `cal-disp --version` and in the `cal_disp_software_version` product
metadata. Check it after a build:

```bash
docker run --rm opera-adt/cal-disp:latest cal-disp --version
```

### Custom Build
```bash
# Custom tag
docker/build-docker-image.sh --tag myorg/cal-disp:v1.0

# Custom user ID (match your host user)
docker/build-docker-image.sh --user-id $(id -u)

# Custom base image
docker/build-docker-image.sh --base ubuntu:24.04

# Explicit version (a release build from a tagged checkout gives e.g. 1.0.0 anyway)
docker/build-docker-image.sh --version 1.0.0

# Combined
docker/build-docker-image.sh --tag myorg/cal-disp:v1.0 --user-id $(id -u) --base ubuntu:24.04
```

### Build Options
```
-t, --tag TAG        Image name/tag (default: opera-adt/cal-disp:latest)
-u, --user-id ID     User ID in container (default: 1000)
-b, --base BASE      Base image (default: ubuntu:22.04)
-v, --version VER    PEP 440 version baked into the image (default: from git describe)
-h, --help           Show help
```

Plain `docker build` works too, but pass the version yourself:

```bash
docker build --network=host -f docker/Dockerfile \
    --build-arg CAL_DISP_VERSION=$(docker/pep440-version.sh) -t opera-adt/cal-disp:latest .
```

## Running

The image's working directory is `/home/work`; mount the folder that holds
your inputs and runconfig there so relative paths in the runconfig resolve.

### Show Help
```bash
docker run --rm -it opera-adt/cal-disp:latest cal-disp --help
```

### Process with Runconfig
```bash
docker run --rm -it \
    --user $(id -u):$(id -g) \
    -v $PWD:/home/work \
    opera-adt/cal-disp:latest \
    cal-disp run runconfig.yaml
```

### Interactive Session
```bash
docker run --rm -it \
    --user $(id -u):$(id -g) \
    -v $PWD:/home/work \
    opera-adt/cal-disp:latest \
    bash
```

### Download Products
The `download` group has the sub-commands `disp-s1`, `unr`, `tropo` and
`burst-bounds` (`cal-disp download --help`):
```bash
docker run --rm -it \
    --user $(id -u):$(id -g) \
    -v $PWD:/home/work \
    opera-adt/cal-disp:latest \
    cal-disp download unr --frame-id 8882 -o unr_data

docker run --rm -it \
    --user $(id -u):$(id -g) \
    -v $PWD:/home/work -v $HOME/.netrc:/home/conda/.netrc:ro \
    opera-adt/cal-disp:latest \
    cal-disp download disp-s1 --frame-id 8882 -o disp_data --start 2022-07-01 --end 2022-07-31
```

## Docker Run Flags

- `--rm` - Remove container after exit
- `-it` - Interactive terminal
- `--user $(id -u):$(id -g)` - Run as your user (avoids permission issues)
- `-v $PWD:/home/work` - Mount current directory to `/home/work`, the image's working directory

Matplotlib's config/cache directory is `/home/conda/.config/matplotlib`
(world-writable), so browse-image generation works for any `--user` even when
`/tmp` is read-only.

## Reproducible Environments

The Docker image is built in two layers:

1. the conda layer from `docker/conda-lock.txt`, an `@EXPLICIT` lockfile with
   exact package URLs and md5 checksums, generated from `environment.yml`;
2. `pip install cal-disp[download]` on top, which only adds the pure-Python
   dependencies (and the pinned Venti commit from `pyproject.toml`).

Compiled packages must come from the lockfile, never from PyPI wheels: the
`netCDF4` wheel bundles its own HDF5, which does not match the HDF5 that
`h5py`/GDAL link against and can crash when several threads read files.
`environment.yml` lists `netcdf4` for that reason, and the lockfile must keep
it (a lock generated from an older `environment.yml` omitted it, so pip pulled
the wheel).

### Updating Dependencies
```bash
# 1. Edit environment.yml
vim environment.yml

# 2. Regenerate the lockfile (linux-64). Without Docker, use the local mamba:
docker/create-lockfile.sh --file environment.yml > docker/conda-lock.txt
docker/create-lockfile.sh --file environment.yml --no-docker > docker/conda-lock.txt

# 3. Check it: python 3.13.x, and netcdf4/libnetcdf/hdf5/h5py all present
grep -E '/(python|netcdf4|libnetcdf|hdf5|h5py)-[0-9]' docker/conda-lock.txt

# 4. Rebuild image
docker/build-docker-image.sh

# 5. Commit both files
git add environment.yml docker/conda-lock.txt
git commit -m "Update dependencies"
```

`--no-docker` solves with `mamba` (or `micromamba`/`conda`) on the local
machine into a scratch prefix (under `$TMPDIR`) and removes it afterwards; it
must run on linux-64, since that is the platform the image is built for. The
Docker mode solves in a clean `mambaorg/micromamba` container. Both write the
same format (4 header lines, then the sorted package URLs); the committed
lock was generated with `--no-docker`.

To check a lock before building the image, create an environment from it and
confirm the HDF5 versions agree:

```bash
mamba create -p /tmp/lockcheck --file docker/conda-lock.txt
/tmp/lockcheck/bin/python -c "import netCDF4, h5py; print(netCDF4.__hdf5libversion__, h5py.version.hdf5_version)"
```

### Why Lockfiles?

**Problem**: `environment.yml` with version ranges installs different packages over time.

**Solution**: Lockfile pins exact package URLs with hashes.
```
# docker/conda-lock.txt
@EXPLICIT
https://conda.anaconda.org/conda-forge/linux-64/numpy-2.5.3-py313hb5f73ae_0.conda#efea38e4c681d5c33aade39330bea01c
```

**Benefits**:
- Identical builds months apart
- Fast, deterministic installs
- Verifiable with checksums

### Creating Lockfiles
```bash
# From environment.yml
docker/create-lockfile.sh --file environment.yml > docker/conda-lock.txt

# With extra packages
docker/create-lockfile.sh --file environment.yml --pkgs pytest black > docker/conda-lock.txt
```

## Troubleshooting

### Permission Denied on Output Files

**Problem**: Files created by container are owned by root.

**Solution**: Use `--user $(id -u):$(id -g)` flag.
```bash
docker run --user $(id -u):$(id -g) -v $PWD:/home/work ...
```

### Runconfig Not Found

**Problem**: `cal-disp run runconfig.yaml` reports the file does not exist.

**Solution**: The image's working directory is `/home/work`; mount your folder
there (`-v $PWD:/home/work`), not at `/work`.

### Network Issues During Build

**Problem**: Cannot reach package repositories.

**Solution**: Check network settings, try `--network=host`:
```bash
docker build --network=host ...
```

Already included in `docker/build-docker-image.sh`.

### Build Cache Issues

**Problem**: Old dependencies cached.

**Solution**: Force rebuild:
```bash
docker build --no-cache ...
```
