# Development

## Environment

```bash
pixi install -e dev        # runtime + test + lint + notebook tooling
pixi shell -e dev
```

or with conda: `mamba env create -f environment.yml`, then
`pip install -e ".[download,test]"` and `pre-commit install`.

Python is bounded to `>=3.11,<3.14`. Venti is pinned to a commit in
`pyproject.toml`; bumping it means regenerating both lock files and the
golden dataset.

## Tests and linting

```bash
pixi run test              # pytest, from the repository root
pixi run lint              # pre-commit: ruff, black, mypy, hygiene hooks
```

There is a single pytest configuration, in `pyproject.toml`. Warnings are
errors; the few upstream warnings that are known to be benign are listed
there with a reason. Tests never touch the network: download tests run
against a local HTTP server, and workflow tests block sockets.

## Lock files

| File | Used by | Regenerate with |
|---|---|---|
| `pixi.lock` | `pixi install` | `pixi install` after editing `pyproject.toml` |
| `docker/conda-lock.txt` | the Docker image (conda layer only) | `docker/create-lockfile.sh --file environment.yml --no-docker > docker/conda-lock.txt` |

`netcdf4`, `h5py` and GDAL must all come from conda-forge so that they share
one HDF5 library; check with the
[environment provenance notebook](notebooks/12_environment_provenance.ipynb).

## Docker

`docker/build-docker-image.sh` passes the git version into the image
(`--build-arg CAL_DISP_VERSION`), so products report the real
`cal_disp_software_version`. Mount the working folder at `/home/work`.

## Golden dataset

`scripts/build_golden_output.sh` downloads the inputs for frame F08882,
writes the configs and runs the calibration once; `scripts/run_validation.sh`
re-runs it and compares with `cal-disp validate`. Rebuild the golden product
whenever a change alters the output on purpose (algorithm parameters,
product encoding or metadata) and say so in the commit message.

## Notebooks

`docs/notebooks/` holds executed notebooks that double as regression
evidence:

- `00_calibration_walkthrough.ipynb` — the workflow step by step, ending in
  the golden comparison.
- `01`–`13` — one notebook per defect fixed in the gamma-release review,
  each showing the behaviour before and after the fix.

They read the golden dataset from `test_golden/` next to the repository root,
or from `$CAL_DISP_GOLDEN_DIR`, and write to `$TMPDIR/cal_disp_notebooks/`.
Re-execute one with:

```bash
pixi run -e dev jupyter nbconvert --to notebook --execute --inplace \
    docs/notebooks/02_float32_tropo_delay.ipynb
```

On hosts with fragmented transparent huge pages, set
`NUMPY_MADVISE_HUGEPAGE=0`; otherwise full-frame runs can be many times
slower.

## Planned work

Deferred items, including those that would change the golden dataset, are
tracked in [TODO.md](https://github.com/opera-adt/cal-disp/blob/main/TODO.md).

## Release checklist

1. `pixi run test` and `pixi run lint` pass.
2. Lock files are current and `docker/conda-lock.txt` contains `netcdf4`.
3. Golden dataset rebuilt with the release code; `scripts/run_validation.sh`
   passes; the walkthrough notebook's final comparison passes.
4. Tag `vX.Y`; build the Docker image from the tag so the product metadata
   carries that version.
