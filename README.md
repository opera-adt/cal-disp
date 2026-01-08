[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/opera-adt/cal-disp/main.svg)](https://results.pre-commit.ci/latest/github/opera-adt/cal-disp/main)

# CAL-DISP

OPERA Calibration for DISP, a SAS repository.
Calibration workflows for OPERA DISP products.

Creates the science application software (SAS) using the [Venti]() library.

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
   # Install venti (add at later stage)
   python -m pip install -e venti/

   # Install cal-disp with download capabilities
   python -m pip install -e "cal-disp[download]"

   # Or basic install only
   python -m pip install -e cal-disp/
```

**Docker:** See [Docker_README](docker/README.md)

## Quick Start
```bash
# Download data
cal-disp download disp-s1 --frame-id 8882 -s 2022-07-01 -e 2022-08-01 -o ./data
cal-disp download tropo -i ./data/OPERA_L3_DISP-S1_*.nc -o ./data/tropo
cal-disp download unr --frame-id 8882 -o ./data/unr

# Create config and run
cal-disp config \
    -d ./data/OPERA_L3_DISP-S1_*.nc \
    -cl ./data/unr/grid_latlon_lookup_v0.2.txt \
    -cd ./data/unr \
    --los-file ./data/static/line_of_sight_enu.tif \
    --dem-file ./data/static/dem.tif \
    -c runconfig.yaml

cal-disp run runconfig.yaml
```

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
