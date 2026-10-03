[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/opera-adt/cal-disp/main.svg)](https://results.pre-commit.ci/latest/github/opera-adt/cal-disp/main)

# CAL-DISP

OPERA Calibration for DISP, a SAS repository.
Calibration workflows for OPERA DISP products.

Creates the science application software (SAS) using the [Venti](https://github.com/opera-adt/Venti) library.

## Installation

### Prerequisites

- Python 3.11–3.13 (`environment.yml` bounds it: the conda-forge `netcdf4`
  build for 3.14 warns about the numpy ABI)
- mamba/conda

### Setup

1. **Clone the repository:**
```bash
git clone https://github.com/opera-adt/cal-disp.git
```

2. **Create environment:**
```bash
mamba env create --name my-cal-env --file cal-disp/environment.yml
conda activate my-cal-env
```

3. **Install the package:**
```bash
# Install cal-disp with download capabilities
python -m pip install -e "cal-disp[download]"

# Or basic install only
python -m pip install -e cal-disp/
```

This installs [Venti](https://github.com/opera-adt/Venti) at the commit
pinned in `pyproject.toml` (the one the golden dataset and the regression
tests were produced with). To develop against a Venti checkout instead,
`pip install -e venti/` afterwards. Keep the compiled packages (`netcdf4`,
`h5py`, `gdal`, `rasterio`, ...) from conda-forge as `environment.yml` lists
them: the PyPI `netCDF4` wheel bundles a different HDF5 than `h5py`/GDAL link
against, which can crash multi-threaded reads.

### Setup with pixi

[pixi](https://pixi.sh) creates the same conda + pip environment from the
`[tool.pixi]` tables in `pyproject.toml`, locked in `pixi.lock`, with no
manual activation step:

```bash
git clone https://github.com/opera-adt/cal-disp.git
cd cal-disp
pixi install            # runtime environment with the download extra
pixi run cal-disp --help

pixi run -e dev test    # run the test suite in the dev environment
pixi run -e dev lint    # pre-commit (ruff, black, mypy) on all files
pixi shell -e dev       # open a shell with the dev environment activated
```

Other tasks: `pixi run build-golden` and `pixi run validate --golden-dir DIR`
wrap the scripts in `scripts/`. Re-run `pixi install` after pulling changes
to `pyproject.toml` or `pixi.lock`.

**Docker:** See [docker/README.md](docker/README.md). The image installs the
conda layer from `docker/conda-lock.txt`; regenerate that lock with
`docker/create-lockfile.sh --file environment.yml [--no-docker] > docker/conda-lock.txt`
whenever `environment.yml` changes.

---

## Quick Start

### 1. Prepare static layers (LOS and DEM)

Use the [Venti staging scripts](https://github.com/opera-adt/Venti/tree/main/scripts/staging/) to download the static layers for your frame:

```bash
python los_cli.py --frame-id 8882
python dem_cli.py --frame-id 8882
```

### 2. Download DISP-S1 products

```bash
cal-disp download disp-s1 \
    --frame-id 8882 \
    --start 2022-07-01 \
    --end 2022-08-01 \
    -o ./data/disp
```

Use `-n` to parallelize downloads:
```bash
cal-disp download disp-s1 --frame-id 8882 -o ./data/disp -n 8
```

### 3. Download ancillary data

**UNR GNSS timeseries** (required):
```bash
cal-disp download unr --frame-id 8882 -o ./data/unr
```

This stages UNR's `constant` grid (precomputed linear rates, IGS20). For the
time-variable positions use `--grid-type variable`. Files are named
`<id>_IGS20_<grid_type>.tenv8`, and the grid type must match both
`unr_grid_type` in the runconfig and `grid_type` in the algorithm parameters.

**Tropospheric corrections** (optional, one file per DISP acquisition):
```bash
cal-disp download tropo \
    -i ./data/disp/OPERA_L3_DISP-S1_IW_F08882_VV_*.nc \
    -o ./data/tropo

# Add --interp to download 2 scenes per date for temporal interpolation
cal-disp download tropo -i ./data/disp/OPERA_L3_DISP-S1_*.nc -o ./data/tropo --interp
```

**Burst boundary tiles** (optional):
```bash
cal-disp download burst-bounds \
    -i ./data/disp/OPERA_L3_DISP-S1_IW_F08882_VV_*.nc \
    -o ./data/burst_bounds
```

### 4. Configure the workflow

Copy [configs/algorithm_parameters.yaml](configs/algorithm_parameters.yaml) and edit as needed, then generate the run config:

```bash
cal-disp config \
    --disp-file ./data/disp/OPERA_L3_DISP-S1_IW_F08882_VV_*.nc \
    --frame-id 8882 \
    --unr-grid-latlon ./data/unr/grid_latlon_lookup_v0.2.txt \
    --unr-grid-dir ./data/unr \
    --unr-grid-version 0.2 \
    --unr-grid-type constant \
    --algorithm-params configs/algorithm_parameters.yaml \
    --los-file line_of_sight_enu.tif \
    --dem-file dem.tif \
    --output-dir outputs/ \
    --work-dir scratch/
```

This writes `runconfig.yaml` into the work directory (`--work-dir`). Use
`-c PATH` to write it elsewhere: the path is used as given, directories included
(e.g. `-c configs/runconfig.yaml`).

**Optional flags:**

| Flag | Description |
|---|---|
| `--ref-tropo-files` | TROPO files for reference date (repeat for multiple) |
| `--sec-tropo-files` | TROPO files for secondary date (repeat for multiple) |
| `--mask-file` | Byte mask (0=invalid, 1=good); recorded but **not applied** in this release (the run logs a warning) |
| `--algorithm-overrides` | Frame-specific parameter overrides (JSON), applied on top of `--algorithm-params`, e.g. `{"8882": {"downsample_factor": 3}}` |
| `--defo-area-db` | Deforming areas database (GeoJSON) |
| `--event-db` | Events database (GeoJSON) |
| `-w` / `--n-workers` | Number of parallel workers (default: 4) |
| `-t` / `--threads-per-worker` | Threads per worker (default: 1) |
| `--block-shape` | Processing block size in pixels (default: 512 512) |

### 5. Run the workflow

```bash
cal-disp run scratch/runconfig.yaml
```

Enable debug logging with:
```bash
cal-disp --debug run scratch/runconfig.yaml
```

### 6. Validate output (optional)

Compare a new output against a reference product:

```bash
cal-disp validate reference.nc output.nc

# Custom tolerance (default 1e-6): --tolerance sets both --rtol and --atol
cal-disp validate reference.nc output.nc --tolerance 1e-5
cal-disp validate reference.nc output.nc --rtol 1e-5 --atol 1e-7

# Only the root, identification and metadata groups (skip the auxiliary group)
cal-disp validate reference.nc output.nc --group main
```

The comparison covers every group (a missing group is a failure), each
variable's dtype, shape, attributes and values (reference zeros included), the
CRS and grid transform, the identification and metadata values, and the browse
PNG next to the output (at most 2048 px per side). Only build versions, the
processing start time and the embedded runconfig are ignored. All differences
are listed and the command exits with status 1 if there is any.

### Golden dataset

Build the golden dataset (frame F08882) and validate it. Requires the
`download` extra and Earthdata credentials in `~/.netrc`:

```bash
scripts/build_golden_output.sh --output-dir golden_data_disp_cal
scripts/run_validation.sh --golden-dir golden_data_disp_cal
```

Or, with the Docker image, from the folder that holds `golden_data_disp_cal/`:

```bash
docker run --rm --user $(id -u):$(id -g) -v $PWD:/home/work -w /home/work <image> \
    cal-disp run golden_data_disp_cal/configs/runconfig.yaml
docker run --rm -v $PWD:/home/work -w /home/work <image> \
    opera_cal-disp validate golden_data_disp_cal/golden_output/<golden>.nc \
    golden_data_disp_cal/output/<test>.nc
```

A golden product validates only in the software environment that made it: for a
delivery, build it inside the delivered image (see the script header).

---

## Algorithm Parameters

Key parameters in `configs/algorithm_parameters.yaml`:

| Parameter | Default | Description |
|---|---|---|
| `grid_type` | `constant` | GNSS model type: `constant` (velocity field) or `variable` (epoch-specific) |
| `reference_frame` | `IGS20` | GNSS reference frame: `IGS20` or `IGS14` |
| `window_size_meters` | `30000.0` | Side length of the moving window for plane fitting (metres) |
| `posting_meters` | `30.0` | DISP pixel spacing (30 m for DISP-S1) |
| `downsample_factor` | `6` | Integer downsampling before surface fitting (1 = disabled) |
| `calibration_surface_smoothing_method` | `gaussian` | Smoothing filter: `gaussian`, `gaussian_fft`, `hanning_fft`, `savitzky_golay` |
| `unwrap_error_correction` | `false` | Venti region-offset unwrap-error correction (off: its mask-island segmentation quantises real signal into λ/2 steps) |

---

## Development

```bash
# Install with dev dependencies
python -m pip install -e ".[download,test]"

# Set up pre-commit hooks
pre-commit install

# Run tests
pytest
```

We use:
- [black](https://black.readthedocs.io/) for formatting
- [ruff](https://docs.astral.sh/ruff/) for linting
- [mypy](https://mypy-lang.org/) for type checking
- [numpydoc](https://numpydoc.readthedocs.io/) for docstrings

Pre-commit runs these automatically. If checks fail, fix the issues and re-add files before committing.

## Contributing

1. Fork the repo
2. Create a branch (`git checkout -b feature/my-feature`)
3. Make changes and add tests
4. Run `pytest` and ensure pre-commit passes
5. Open a PR

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

BSD-3-Clause or Apache-2.0 - see `LICENSE.txt`

---

Developed at JPL/Caltech under contract with NASA.
Questions? [Open an issue](https://github.com/opera-adt/cal-disp/issues)
