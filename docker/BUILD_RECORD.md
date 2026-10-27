# Image build record

| Tag | Built | Host | Source | Size | Image ID |
|---|---|---|---|---|---|
| `cal-disp:0.5.4` | 2026-10-22 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.4`, `--network=host`) | tag `v0.5.4` (`c3e96bb`), version `0.5.4` | 1.28 GB | `sha256:ab025a6ffbacb6389f427f301f7c2eb7b89b8936dc49f9b235ebec69e4539161` |
| `cal-disp:0.5.6` | 2026-10-25 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.6`, `--network=host`) | tag `v0.5.6` (`8d96246`), version `0.5.6` | 1.27 GB | `sha256:36e18d977296c87a0dc0cca48f8b9bba9cee7d7f9b5044e090baf08ebc0ca5ee` |
| `cal-disp:0.5.6-rc` | 2026-10-25 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.6-rc -v 0.5.6rc1`, `--network=host`) | branch `release/0.5.6-fit-mask-file` (`bd64e58`), version `0.5.6rc1` | 1.27 GB | `sha256:6caaa8c02569029e8db1b45ccd6ac6a8102633aabd584c5648f03f36fc7f9fd3` |
| `cal-disp:0.5.5` | 2026-10-22 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.5`, `--network=host`) | tag `v0.5.5` (`c1f3579`), version `0.5.5` | 1.27 GB | `sha256:de6dc5d25efdf86b1e6bdd16e36e0c943231ca0625253f8474151a170cb57f39` |
| `cal-disp:0.5.5-rc` | 2026-10-22 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.5-rc -v 0.5.5rc1`, `--network=host`) | branch `release/0.5.5-plate-motion` (`d35b0ab`, #33 + #34), version `0.5.5rc1` | 1.27 GB | `sha256:4d1dc8bb3ae74153a8405d04699ff59056b6af42882ffe223bcb2cac6efb022f` |
| `cal-disp:0.5.4-rc` | 2026-10-22 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.4-rc -v 0.5.4rc1`, `--network=host`) | branch `release/0.5.4-hawaii-stratified`, version `0.5.4rc1` | 1.28 GB | `sha256:3a6081ed2274475ee1fc3d04024413b4224fffa027d543911f98f4b676dd7c55` |
| `cal-disp:0.5.3` | 2026-10-20 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.3`, `--network=host`) | tag `v0.5.3` (`3f2c241`), version `0.5.3` | 1.28 GB | `sha256:37b1b94432b12044d77efc8fb00c69fc83a67439de5fce5046c4765424edcad8` |
| `cal-disp:0.5.3-rc` | 2026-10-20 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.3-rc -v 0.5.3rc1`, `--network=host`) | branch `release/0.5.3-stratified-tropo`, version `0.5.3rc1` | 1.28 GB | `sha256:6baadd83b52bb2d8de823ef86ec70aa1b42e601070f78027f0ab22f572d28724` |
| `cal-disp:0.5.2` | 2026-10-20 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.2`, `--network=host`) | tag `v0.5.2` (`1d64f2f`), version `0.5.2` | 1.28 GB | `sha256:b4b0d7047299cdf8e5484dd5a65f607bf33f43823e90616e84246d712e52d2a3` |
| `cal-disp:0.5.2-rc` | 2026-10-20 | aurora (`docker/build-docker-image.sh -t cal-disp:0.5.2-rc -v 0.5.2rc1`, `--network=host`) | branch `release/0.5.2-frame-f18904` `9f05f5c`, version `0.5.2rc1` | 1.28 GB | `sha256:4225dde676de01a372a97365d08a4e8983063592020c59783241e01137bdab34` |
| `cal-disp:0.5.1` | 2026-10-16 | aurora (`docker/build-docker-image.sh`, `--network=host`) | tag `v0.5.1` (`54e6382`), version `0.5.1` | 1.28 GB | `sha256:f333f7bfc62e09322317e2bf6c940ac82f5267adf41385e7717d76fe3a8866c9` |
| `cal-disp:0.5.1-rc` | 2026-10-16 | aurora (`docker/build-docker-image.sh -v 0.5.1rc1`, `--network=host`) | branch `feature/t37-3-metadata` `df23416`, version `0.5.1rc1` | 1.28 GB | `sha256:a7d61bc66b1b3d612d43a52cd753eac62d546237eb4fcaaffe44daa08052ae2b` |
| `cal-disp:0.5.0` | 2026-10-15 | aurora (`docker/build-docker-image.sh`, `--network=host`) | tag `v0.5.0` (`438b84e`), version `0.5.0` | 1.28 GB | `sha256:f29cd9131c4000f719207759645111ba497ab1d077aca029c3e605161a04ffe1` |
| `cal-disp:0.5.0-rc` | 2026-10-15 | aurora (`docker/build-docker-image.sh -v 0.5.0rc1`, `--network=host`) | T49 branch `54fe893` + the Dockerfile change dropping the `download` extra (same commit as this record), version `0.5.0rc1` | 1.28 GB | `sha256:d3f0fafcabaf6d2d0898e6c5f4402a0423b42f4097daeca4dc88d95493fadc18` |
| `cal-disp:0.3.0-rc` | 2026-10-05 | aurora (`docker/build-docker-image.sh`, `--network=host`) | `gamma-release` `7085ec9`, version `0.2.post1.dev79+g7085ec96e` | 1.5 GB | `sha256:fc115560ebf97907c112d7096a7fb7796506e794122a3c8e369c87d969691333` |

Key versions inside: Python 3.13.15, GDAL 3.13.3, PROJ 9.8.1, rasterio 1.5.1,
numpy 2.5.3, scipy 1.18.1, opera-utils 0.25.8, venti `2a7e61f`. Full conda
list: `installed_packages_0.3.0rc.csv` (214 packages).

## Validation inside the image

```bash
docker run --rm --user $(id -u):$(id -g) -e NUMPY_MADVISE_HUGEPAGE=0 -e HOME=/tmp \
  -v $PWD:/home/work -w /home/work cal-disp:0.3.0-rc \
  bash -lc 'eval "$(micromamba shell hook --shell bash)" && micromamba activate cal-disp \
            && scripts/run_validation.sh --golden-dir test_golden'
```

2026-10-05: **passed at 1e-6** against the golden built on aurora with the
pixi `dev` env (GDAL 3.11 line); peak memory 6.31 GB. The golden is therefore
reproducible across the two environments for the gamma configuration.

## 0.5.0-rc (2026-10-15)

Venti `2d529ae` (tree identical to venti-dev `main` `f1da8aa`). Dependency
budget: OK (`scripts/check_env_budget.py`; matplotlib is the recorded
exception). The v0.5 golden was built inside this image
(`scripts/build_golden_output.sh --skip-download`) and validated inside it:
passed at 1e-6, peak 5.76 GB (`test_golden/validation_0.5.0.log`). The gamma
file reproduces the archived 0.3 golden's data layers at 1e-6.

## 0.5.0 (2026-10-15)

Built from the `v0.5.0` tag; reports `cal-disp 0.5.0`, Venti `2d529ae`. Budget
check OK. The v0.5 golden (built in `0.5.0-rc`) validates in this image at
1e-6, peak 5.75 GB (`test_golden/validation_0.5.0_release_image.log`).

## 0.5.1 (2026-10-16)

Built from the `v0.5.1` tag; reports `cal-disp 0.5.1`, Venti `35ec67e`. Budget
check OK. The 0.5.1 golden (built in `0.5.1-rc`) validates in this image at
1e-6, peak 5.76 GB (`test_golden/validation_0.5.1_release_image.log`).

Note (2026-10-20): the local images listed above for 0.3.0-rc through 0.5.2-rc were
removed from aurora by a Docker prune outside this build; only `cal-disp:0.5.2` is
local. None had been pushed; each is reproducible from its tag or branch.
