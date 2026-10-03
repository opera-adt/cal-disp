# Quickstart

This page takes one DISP-S1 frame from a fresh checkout to a validated
DISP-CAL-S1 product. Frame **F08882** is used throughout; it is the frame of
the golden dataset, so every command here can be checked against a known
result.

## 1. Install

=== "pixi (recommended)"

    [pixi](https://pixi.sh) creates the locked conda + pip environment from
    `pyproject.toml` and `pixi.lock`:

    ```bash
    git clone https://github.com/opera-adt/cal-disp.git
    cd cal-disp
    pixi install                 # runtime environment with the download extra
    pixi run cal-disp --help
    pixi shell                   # or prefix every command with `pixi run`
    ```

=== "conda / mamba"

    ```bash
    git clone https://github.com/opera-adt/cal-disp.git
    mamba env create --name cal-disp --file cal-disp/environment.yml
    conda activate cal-disp
    python -m pip install -e "cal-disp[download]"
    ```

    The `download` extra adds `asf_search` and friends for `cal-disp download`.

=== "Docker"

    ```bash
    cd cal-disp
    docker/build-docker-image.sh          # tags the image with the git version
    docker run --rm --user $(id -u):$(id -g) -v $PWD:/home/work -w /home/work \
        cal-disp:latest cal-disp --help
    ```

    The image installs the conda layer from `docker/conda-lock.txt`; see
    [docker/README.md](https://github.com/opera-adt/cal-disp/blob/main/docker/README.md).

Downloads from ASF/Earthdata need credentials in `~/.netrc`:

```text
machine urs.earthdata.nasa.gov login <user> password <password>
```

## 2. Stage the inputs

```bash
mkdir -p data/{disp,static,unr,tropo}

# DISP-S1 product(s) for the frame and date range
cal-disp download disp-s1 --frame-id 8882 --start 2022-01-01 --end 2022-08-01 -o data/disp

# UNR gridded GNSS (constant-velocity grid, IGS20)
cal-disp download unr --frame-id 8882 -o data/unr

# TROPO zenith delay for each DISP acquisition (optional but recommended)
cal-disp download tropo -i data/disp/OPERA_L3_DISP-S1_IW_F08882_VV_*.nc -o data/tropo
```

The DISP-S1-STATIC line-of-sight and DEM layers come from ASF as well; the
golden-dataset script shows the URLs
(`scripts/build_golden_output.sh`), or use the
[Venti staging scripts](https://github.com/opera-adt/Venti/tree/main/scripts/staging/).
Put `*_line_of_sight_enu.tif` and `*_dem.tif` in `data/static`.

!!! note "Grid type"
    `cal-disp download unr` stages the `constant` grid by default. The grid
    type must match both `unr_grid_type` in the runconfig and `grid_type` in
    the algorithm parameters; the workflow refuses to run otherwise.

## 3. Configure

Copy `configs/algorithm_parameters.yaml`, edit it if needed, and generate the
runconfig:

```bash
cal-disp config \
    --disp-file data/disp/OPERA_L3_DISP-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_*.nc \
    --frame-id 8882 \
    --unr-grid-latlon data/unr/grid_latlon_lookup_v0.3.txt \
    --unr-grid-dir data/unr \
    --unr-grid-version 0.3 \
    --unr-grid-type constant \
    --los-file data/static/*_line_of_sight_enu.tif \
    --dem-file data/static/*_dem.tif \
    --ref-tropo-files data/tropo/OPERA_L4_TROPO-ZENITH_20220111T*.nc \
    --sec-tropo-files data/tropo/OPERA_L4_TROPO-ZENITH_20220722T*.nc \
    --algorithm-params configs/algorithm_parameters.yaml \
    --output-dir output \
    --work-dir output/scratch \
    -c runconfig.yaml
```

`cal-disp config` validates every path and writes the runconfig with
explanatory comments. Fields the current release does not support
(`iono_files`, `tiles_files`, non-NetCDF output) are rejected at load time
rather than silently ignored.

## 4. Run

```bash
cal-disp run runconfig.yaml            # add --debug before `run` for DEBUG logs
```

A full frame takes a few minutes and up to ~18 GB of memory. The product,
its browse image and `cal_disp.log` land in `output/`.

## 5. Validate

Against a reference product (e.g. the golden product for the same frame):

```bash
cal-disp validate reference.nc output/OPERA_L4_DISP-CAL-S1_*.nc
```

`validate` compares every group: dimensions, dtypes, data values within
tolerance (`--rtol`/`--atol`, or `--tolerance` for both), variable and global
attributes, the CRS and geotransform, and the identification/metadata
values, ignoring only fields that legitimately differ between runs
(production time, software versions, embedded runconfig).

## 6. The golden dataset

The golden dataset is frame F08882, pair 2022-01-11 / 2022-07-22, in the
delivery layout (`configs/`, `input_data/`, `golden_output/`, `output/`).
Build it and check an installation against it:

```bash
scripts/build_golden_output.sh --output-dir golden_data_disp_cal
scripts/run_validation.sh --golden-dir golden_data_disp_cal
```

or, with pixi, `pixi run build-golden --output-dir DIR` and
`pixi run validate --golden-dir DIR`. The
[workflow walkthrough notebook](notebooks/00_calibration_walkthrough.ipynb)
runs the same computation step by step and ends with the same comparison.

!!! warning
    A golden product is reproducible to the validation tolerance only in the
    software environment that produced it. For a delivery, build it inside
    the delivered Docker image.
