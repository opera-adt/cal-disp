# Reproducible base image
ARG BASE=ubuntu:22.04
FROM ${BASE}

# Metadata
LABEL org.opencontainers.image.title="cal-disp"
LABEL org.opencontainers.image.description="OPERA DISP calibration"
LABEL org.opencontainers.image.authors="OPERA ADT Team"
LABEL org.opencontainers.image.licenses="Apache-2.0"
LABEL org.opencontainers.image.url="https://github.com/opera-adt/cal-disp"
LABEL org.opencontainers.image.source="https://github.com/opera-adt/cal-disp"

ARG DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=true

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    git \
    vim \
    bzip2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install micromamba - CRITICAL: Extract to /tmp to avoid corrupting /bin/
RUN set -ex && \
    cd /tmp && \
    curl -L https://micro.mamba.pm/api/micromamba/linux-64/latest -o micromamba.tar.bz2 && \
    tar -xjf micromamba.tar.bz2 && \
    mv bin/micromamba /usr/local/bin/micromamba && \
    rm -rf /tmp/micromamba.tar.bz2 /tmp/bin && \
    /usr/local/bin/micromamba --version

# Verify system is intact
RUN ls -la /bin/bash /bin/sh && echo "Filesystem intact"

# Create non-root user
ARG MAMBA_USER_ID=1000
ARG MAMBA_USER_GID=1000
RUN groupadd -g "${MAMBA_USER_GID}" --system conda && \
    useradd -l -u "${MAMBA_USER_ID}" -g "${MAMBA_USER_GID}" --system \
        -d /home/conda -m -s /bin/bash conda && \
    mkdir -p /home/work && \
    chown -R conda:conda /home/work /home/conda

USER ${MAMBA_USER_ID}
WORKDIR /home/work

# run commands in a bash login shell
SHELL ["/bin/bash", "-l", "-c"]

# Environment setup
ENV HOME=/home/conda
ENV MAMBA_ROOT_PREFIX=/home/conda/micromamba
ENV PATH="${MAMBA_ROOT_PREFIX}/envs/cal-disp/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=cal-disp
ENV MPLCONFIGDIR=/tmp/matplotlib

# Copy lockfile and source code
COPY --chown=conda:conda docker/conda-lock.txt /tmp/conda-lock.txt
COPY --chown=conda:conda . /tmp/cal-disp

# Install from lockfile
RUN micromamba create -y -n cal-disp -f /tmp/conda-lock.txt && \
    micromamba clean -afy

# Install cal-disp package with download extras
ENV SETUPTOOLS_SCM_PRETEND_VERSION_FOR_CAL_DISP=1.0.0
RUN eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate cal-disp && \
    python -m pip install --no-cache-dir /tmp/cal-disp[download] && \
    rm -rf /tmp/cal-disp /tmp/conda-lock.txt

CMD ["cal-disp", "--help"]
