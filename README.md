[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/opera-adt/cal-disp/main.svg)](https://results.pre-commit.ci/latest/github/opera-adt/cal-disp/main)

# CAL-DISP

OPERA Calibration for DISP, a SAS repository.
Calibration workflows for OPERA DISP products.

Creates the science application software (SAS) using the [Venti](https://github.com/opera-adt/Venti) library.

## Installation

### Prerequisites

- Python ≥3.11
- mamba/conda

### Setup

1. **Clone repositories:**
```bash
git clone https://github.com/opera-adt/venti.git
git clone https://github.com/opera-adt/cal-disp.git
```

2. **Create environment:**
```bash
mamba env create --name my-cal-env --file cal-disp/environment.yml
conda activate my-cal-env
```

3. **Install packages:**
```bash
# Install venti
python -m pip install -e venti/

# Install cal-disp with download capabilities
python -m pip install -e "cal-disp[download]"

# Or basic install only
python -m pip install -e cal-disp/
```

**Docker:** See [Docker_README](docker/README.md)

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

This writes `runconfig.yaml` to the work directory (`scratch/` by default). Use `-c` to set a custom output path.

**Optional flags:**

| Flag | Description |
|---|---|
| `--ref-tropo-files` | TROPO files for reference date (repeat for multiple) |
| `--sec-tropo-files` | TROPO files for secondary date (repeat for multiple) |
| `--iono-files` | Ionospheric correction files |
| `--tiles-files` | Calibration tile bounds files |
| `--mask-file` | Byte mask (0=invalid, 1=good) |
| `--algorithm-overrides` | Frame-specific parameter overrides (JSON) |
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

# Custom tolerance
cal-disp validate reference.nc output.nc --tolerance 1e-5

# Validate only the main data group
cal-disp validate reference.nc output.nc --group main
```

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
| `unwrap_error_correction` | `true` | Apply watershed-based unwrap-error correction |

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
