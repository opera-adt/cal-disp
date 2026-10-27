# Reproducible base image
ARG BASE=ubuntu:22.04
FROM ${BASE}

# Metadata
LABEL org.opencontainers.image.title="cal-disp"
LABEL org.opencontainers.image.description="OPERA DISP Calibration"
LABEL org.opencontainers.image.authors="OPERA ADT Team"
LABEL org.opencontainers.image.licenses="Apache-2.0"
LABEL org.opencontainers.image.url="https://github.com/opera-adt/cal-disp"
LABEL org.opencontainers.image.source="https://github.com/opera-adt/cal-disp"

ARG DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=true

# Install system dependencies (git: pip clones the pinned Venti commit)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    vim \
    bzip2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install a pinned micromamba release (a moving "latest" makes the base
# layer differ between builds). Release binaries:
# https://github.com/mamba-org/micromamba-releases/releases
ARG MICROMAMBA_VERSION=2.1.1-0
RUN set -ex && \
    curl -fsSL "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" \
        -o /usr/local/bin/micromamba && \
    chmod +x /usr/local/bin/micromamba && \
    /usr/local/bin/micromamba --version

# Create non-root user. The image is also run with `--user $(id -u):$(id -g)`
# (any uid), which must be able to reach the environment under /home/conda
# and to write the matplotlib config/cache dir: /tmp, matplotlib's fallback,
# may be read-only in a PGE container.
ARG MAMBA_USER_ID=1000
ARG MAMBA_USER_GID=1000
RUN groupadd -g "${MAMBA_USER_GID}" --system conda && \
    useradd -l -u "${MAMBA_USER_ID}" -g "${MAMBA_USER_GID}" --system \
        -d /home/conda -m -s /bin/bash conda && \
    mkdir -p /home/work /home/conda/.config/matplotlib && \
    chown -R conda:conda /home/work /home/conda && \
    chmod 755 /home/conda && \
    chmod -R a+rwX /home/conda/.config

USER ${MAMBA_USER_ID}
WORKDIR /home/work

# run commands in a bash login shell
SHELL ["/bin/bash", "-l", "-c"]

# Environment setup
ENV HOME=/home/conda
ENV MAMBA_ROOT_PREFIX=/home/conda/micromamba
ENV PATH="${MAMBA_ROOT_PREFIX}/envs/cal-disp/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=cal-disp
ENV MPLCONFIGDIR=/home/conda/.config/matplotlib

# Install the conda layer from the lockfile (docker/create-lockfile.sh)
COPY --chown=conda:conda docker/conda-lock.txt /tmp/conda-lock.txt
RUN micromamba create -y -n cal-disp -f /tmp/conda-lock.txt && \
    micromamba clean -afy && \
    rm -f /tmp/conda-lock.txt

# Copy only what `pip install` needs (keeps the layer small: no data, docs,
# tests or notebooks; .dockerignore drops .git, so the version cannot be
# read from git and must be passed in).
RUN mkdir -p /tmp/cal-disp
COPY --chown=conda:conda pyproject.toml README.md LICENSE /tmp/cal-disp/
COPY --chown=conda:conda requirements /tmp/cal-disp/requirements
COPY --chown=conda:conda src /tmp/cal-disp/src

# Version reported by `cal-disp --version` and written into product metadata
# (cal_disp_software_version). build-docker-image.sh passes the PEP 440
# version of the checkout (docker/pep440-version.sh); a plain `docker build`
# without --build-arg gets a placeholder that is visibly not a release.
ARG CAL_DISP_VERSION=0.0.0.dev0+unknown
ENV SETUPTOOLS_SCM_PRETEND_VERSION_FOR_CAL_DISP=${CAL_DISP_VERSION}
RUN eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate cal-disp && \
    python -m pip install --no-cache-dir "/tmp/cal-disp[download]" && \
    rm -rf /tmp/cal-disp && \
    cal-disp --version

CMD ["cal-disp", "--help"]
