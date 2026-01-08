# cal-disp Docker

Docker container for OPERA DISP-S1 calibration tools.

## Quick Start
```bash
# Build
./build-docker-image.sh -t opera_adt/cal_disp:latest

# Run
docker run --rm -it opera-adt/cal-disp:latest cal-disp --help
```

## Building

### Basic Build
```bash
./build-docker-image.sh
```

Builds `opera-adt/cal-disp:latest` using default settings (Ubuntu 22.04, user ID 1000).

### Custom Build
```bash
# Custom tag
./build-docker-image.sh --tag myorg/cal-disp:v1.0

# Custom user ID (match your host user)
./build-docker-image.sh --user-id $(id -u)

# Custom base image
./build-docker-image.sh --base ubuntu:24.04

# Combined
./build-docker-image.sh --tag myorg/cal-disp:v1.0 --user-id $(id -u) --base ubuntu:24.04
```

### Build Options
```
-t, --tag TAG        Image name/tag (default: opera-adt/cal-disp:latest)
-u, --user-id ID     User ID in container (default: 1000)
-b, --base BASE      Base image (default: ubuntu:22.04)
-h, --help           Show help
```

## Running

### Show Help
```bash
docker run --rm -it opera-adt/cal-disp:latest cal-disp --help
```

### Process with Runconfig
```bash
docker run --rm -it \
    --user $(id -u):$(id -g) \
    -v $PWD:/work \
    opera-adt/cal-disp:latest \
    cal-disp run runconfig.yaml
```

### Interactive Session
```bash
docker run --rm -it \
    --user $(id -u):$(id -g) \
    -v $PWD:/work \
    opera-adt/cal-disp:latest \
    bash
```

### Download Products
```bash
docker run --rm -it \
    --user $(id -u):$(id -g) \
    -v $PWD:/work \
    opera-adt/cal-disp:latest \
    cal-disp download --frame 12345 --start 2024-01-01 --end 2024-12-31
```

## Docker Run Flags

- `--rm` - Remove container after exit
- `-it` - Interactive terminal
- `--user $(id -u):$(id -g)` - Run as your user (avoids permission issues)
- `-v $PWD:/work` - Mount current directory to `/work` in container

## Reproducible Environments

The Docker image is built from a conda lockfile for reproducibility.

### Updating Dependencies
```bash
# 1. Edit environment.yml
vim environment.yml

# 2. Regenerate lockfile
./create-lockfile.sh --file environment.yml > conda-lock.txt

# 3. Rebuild image
./build-docker-image.sh

# 4. Commit both files
git add environment.yml conda-lock.txt
git commit -m "Update dependencies"
```

### Why Lockfiles?

**Problem**: `environment.yml` with version ranges installs different packages over time.

**Solution**: Lockfile pins exact package URLs with hashes.
```
# conda-lock.txt
@EXPLICIT
https://conda.anaconda.org/.../numpy-1.24.3-...-abc123.tar.bz2
```

**Benefits**:
- Identical builds months apart
- Fast, deterministic installs
- Verifiable with checksums

### Creating Lockfiles
```bash
# From environment.yml
./create-lockfile.sh --file environment.yml > conda-lock.txt

# With extra packages
./create-lockfile.sh --file environment.yml --pkgs pytest black > conda-lock.txt
```

Uses Docker + micromamba to resolve dependencies in a clean environment.

### Platform-Specific

Create separate lockfiles per platform:
```bash
# Linux (default)
./create-lockfile.sh --file environment.yml > conda-lock-linux-64.txt

# macOS (run on Mac)
./create-lockfile.sh --file environment.yml > conda-lock-osx-64.txt
```

## Troubleshooting

### Permission Denied on Output Files

**Problem**: Files created by container are owned by root.

**Solution**: Use `--user $(id -u):$(id -g)` flag.
```bash
docker run --user $(id -u):$(id -g) -v $PWD:/work ...
```

### Network Issues During Build

**Problem**: Cannot reach package repositories.

**Solution**: Check network settings, try `--network=host`:
```bash
docker build --network=host ...
```

Already included in `build-docker-image.sh`.

### Build Cache Issues

**Problem**: Old dependencies cached.

**Solution**: Force rebuild:
```bash
docker build --no-cache ...
```
