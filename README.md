[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/opera-adt/cal-disp/main.svg)](https://results.pre-commit.ci/latest/github/opera-adt/cal-disp/main)

# CAL-DISP
OPERA Calibration for DISP, a SAS repository.
Calibration workflows for OPERA DISP products.

Creates the science application software (SAS) using the [Venti]() library.

## Development setup


### Prerequisite installs
1. Download source code:
```bash
git clone https://github.com/opera-adt/venti.git
git clone https://github.com/opera-adt/cal-disp.git
```
2. Install dependencies, either to a new environment:
```bash
mamba env create --name my-cal-env --file cal-disp/environment.yml
conda activate my-cal-env
```
or install within your existing env with mamba.

3. Install `venti` and `cal-disp` via pip in editable mode
```bash
# Note: add venti at later stage
python -m pip install -e  cal-disp/
```
Install with module needed for download
```bash
pip install -e "cal-disp[download]"
```

#### Option 2: Docker

```bash
# Clone the repository
git clone https://github.com/opera-adt/cal-disp.git
cd cal-disp


# Build Docker image
./scripts/create-lockfile.sh --file environment.yml > docker/conda-lock.txt
./scripts/build-docker.sh
```

See [Docker Setup Guide](docker/README.md) for detailed instructions.

### Basic Usage

```bash
# Download DISP-S1 data for a specific frame
cal-disp download disp-s1 --frame-id 8882 -s 2022-07-01 -e 2022-08-01 -o ./data

# Download tropospheric corrections
cal-disp download tropo -i ./data/OPERA_L3_DISP-S1_*.nc -o ./data/tropo

# Download UNR GPS data
cal-disp download unr --frame-id 8882 -o ./data/unr

# Create calibration configuration
cal-disp config \
    -d ./data/OPERA_L3_DISP-S1_*.nc \
    -cl ./data/unr/grid_latlon_lookup_v0.2.txt \
    -cd ./data/unr \
    --los-file ./data/static/line_of_sight_enu.tif \
    --dem-file ./data/static/dem.tif \
    -c runconfig.yaml

# Run calibration
cal-disp run runconfig.yaml
```

## Installation Details

### Prerequisites

- Python ≥ 3.11
- conda/mamba (recommended) or pip
- For download functionality: internet access to ASF and UNR data sources

### Environment Setup

#### Creating the Environment

**From environment.yml:**
```bash
# Create new environment
mamba env create -f environment.yml

# Or update existing environment
mamba env update -f environment.yml --prune
```

**From scratch:**
```bash
# Create environment with basic dependencies
mamba create -n cal-disp python=3.11 gdal geopandas rasterio -c conda-forge
conda activate cal-disp
```

#### Installing cal-disp

**Basic installation (no download capabilities):**
```bash
pip install -e .
```

**With download extras (recommended):**
```bash
pip install -e ".[download]"
```

This installs additional dependencies:
- `asf_search` - For downloading DISP-S1 products
- `opera-utils` - OPERA utility functions

**With all development tools:**
```bash
pip install -e ".[download,test,docs]"
```

### Verifying Installation

```bash
# Check version
cal-disp --version

# View available commands
cal-disp --help

# Test download command
cal-disp download --help
```

## Development Setup

### Setting Up for Development

1. **Clone and install in development mode:**
   ```bash
   git clone https://github.com/opera-adt/cal-disp.git
   cd cal-disp
   mamba env create -f environment.yml
   conda activate my-cal-disp
   pip install -e ".[download,test]"
   ```

2. **Install pre-commit hooks:**
   ```bash
   pre-commit install
   ```

   This automatically runs:
   - Code formatting with [black](https://black.readthedocs.io/)
   - Linting with [ruff](https://docs.astral.sh/ruff/)
   - Type checking with [mypy](https://mypy-lang.org/)
   - Import sorting with [isort](https://pycqa.github.io/isort/)

3. **Configure your editor:**
   - Install black and ruff plugins
   - Enable format-on-save for automatic formatting
   - Configure mypy for type hints

### Code Style

We follow these conventions:

- **Docstrings**: [numpydoc format](https://numpydoc.readthedocs.io/en/latest/format.html)
- **Formatting**: [black](https://black.readthedocs.io/) (line length: 88)
- **Type hints**: Required for public APIs
- **Imports**: Sorted with [isort](https://pycqa.github.io/isort/)

### Running Tests

Install test dependencies:
```bash
pip install -e ".[test]"
```

Run the test suite:
```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=cal_disp --cov-report=html

# Run specific test file
pytest tests/test_download.py

# Run with verbose output
pytest -v
```

### Setup for contributing

We use [pre-commit](https://pre-commit.com/) to automatically run linting, formatting, and [mypy type checking](https://www.mypy-lang.org/).
Additionally, we follow [`numpydoc` conventions for docstrings](https://numpydoc.readthedocs.io/en/latest/format.html).
To install pre-commit locally, run:

```bash
pre-commit install
```
This adds a pre-commit hooks so that linting/formatting is done automatically. If code does not pass the checks, you will be prompted to fix it before committing.
Remember to re-add any files you want to commit which have been altered by `pre-commit`. You can do this by re-running `git add` on the files.

Since we use [black](https://black.readthedocs.io/en/stable/) for formatting and [flake8](https://flake8.pycqa.org/en/latest/) for linting, it can be helpful to install these plugins into your editor so that code gets formatted and linted as you save.

### Running the unit tests

After making functional changes and/or have added new tests, you should run pytest to check that everything is working as expected.

First, install the extra test dependencies:
```bash
python -m pip install --no-deps -e .[test]
```

Then run the tests:

```bash
pytest
```

## Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Run tests (`pytest`)
5. Ensure pre-commit checks pass
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to the branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## License

This project is licensed under either the BSD-3-Clause or Apache-2.0 license - see the [LICENSE.txt](LICENSE.txt) file for details.

## Acknowledgments

This software was developed at the Jet Propulsion Laboratory, California Institute of Technology, under contract with the National Aeronautics and Space Administration.

## Contact

For questions or support, please [open an issue](https://github.com/opera-adt/cal-disp/issues) on GitHub.
