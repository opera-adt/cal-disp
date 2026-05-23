#!/bin/bash
# build_golden_output.sh
#
# Generate a golden dataset from REAL downloaded OPERA data.
# Use this for manual validation that the pipeline processes actual products
# correctly.  For lightweight automated CI testing use create_golden_dataset.py
# instead, which creates small synthetic inputs without any network access.
#
# Prerequisites
# -------------
#   - cal-disp installed in the active environment
#   - Earthdata credentials configured (for cal-disp download commands)
#   - Static layer GeoTIFFs (LOS ENU + DEM) for the chosen frame
#
# Usage
# -----
#   # Minimal — reads output root from $CAL_DISP_TEST_DATA
#   ./scripts/build_golden_output.sh \
#       --static-dir /path/to/static_layers
#
#   # Fully explicit
#   ./scripts/build_golden_output.sh \
#       --output-dir /path/to/test_data \
#       --frame-id   8882 \
#       --start      2016-07-01 \
#       --end        2016-08-31 \
#       --static-dir /path/to/static_layers \
#       --algo-params configs/algorithm_parameters.yaml \
#       --unr-version 0.2 \
#       --unr-type constant \
#       --skip-tropo
#
# Output layout (under <output-dir>/)
# ------------------------------------
#   golden/
#     disp/           downloaded OPERA_L3_DISP-S1_*.nc
#     gnss/           grid_latlon_lookup.txt + *.tenv8
#     los.tif         (symlinked from --static-dir)
#     dem.tif         (symlinked from --static-dir)
#     algorithm_parameters.yaml
#   golden_output/
#     OPERA_L4_DISP-CAL-S1_*.nc   expected CalProduct
#
# After running, enable integration tests with:
#   export CAL_DISP_TEST_DATA=<output-dir>
#   pytest -m integration

set -euo pipefail

# Defaults
FRAME_ID="8882"
START_DATE="2016-07-01"
END_DATE="2016-08-31"
UNR_VERSION="0.2"
UNR_TYPE="constant"
SKIP_TROPO=false
OUTPUT_DIR="${CAL_DISP_TEST_DATA:-}"
STATIC_DIR=""
ALGO_PARAMS="$(dirname "$0")/../configs/algorithm_parameters.yaml"

# Argument parsing
usage() {
    grep "^#" "$0" | grep -v "^#!/" | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)   OUTPUT_DIR="$2";   shift 2 ;;
        --frame-id)     FRAME_ID="$2";     shift 2 ;;
        --start)        START_DATE="$2";   shift 2 ;;
        --end)          END_DATE="$2";     shift 2 ;;
        --static-dir)   STATIC_DIR="$2";   shift 2 ;;
        --algo-params)  ALGO_PARAMS="$2";  shift 2 ;;
        --unr-version)  UNR_VERSION="$2";  shift 2 ;;
        --unr-type)     UNR_TYPE="$2";     shift 2 ;;
        --skip-tropo)   SKIP_TROPO=true;   shift ;;
        -h|--help)      usage ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate inputs
if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "ERROR: Provide --output-dir or set \$CAL_DISP_TEST_DATA." >&2
    exit 1
fi

if [[ -z "${STATIC_DIR}" ]]; then
    echo "ERROR: --static-dir is required (path to directory with LOS and DEM GeoTIFFs)." >&2
    exit 1
fi

LOS_FILE=$(find "${STATIC_DIR}" -name "*line_of_sight_enu.tif" -type f | head -n 1)
DEM_FILE=$(find "${STATIC_DIR}" -name "*dem.tif" -type f | head -n 1)

if [[ -z "${LOS_FILE}" ]]; then
    echo "ERROR: No *line_of_sight_enu.tif found in ${STATIC_DIR}." >&2
    exit 1
fi
if [[ -z "${DEM_FILE}" ]]; then
    echo "ERROR: No *dem.tif found in ${STATIC_DIR}." >&2
    exit 1
fi

if [[ ! -f "${ALGO_PARAMS}" ]]; then
    echo "ERROR: Algorithm parameters file not found: ${ALGO_PARAMS}" >&2
    exit 1
fi

ALGO_PARAMS="$(realpath "${ALGO_PARAMS}")"
OUTPUT_DIR="$(realpath "${OUTPUT_DIR}")"

GOLDEN_DIR="${OUTPUT_DIR}/golden"
GOLDEN_OUTPUT_DIR="${OUTPUT_DIR}/golden_output"

# Directory structure
echo "=== Building real-data golden dataset ==="
echo "  frame        : ${FRAME_ID}"
echo "  date range   : ${START_DATE} → ${END_DATE}"
echo "  golden inputs: ${GOLDEN_DIR}"
echo "  golden output: ${GOLDEN_OUTPUT_DIR}"
echo ""

mkdir -p \
    "${GOLDEN_DIR}/disp" \
    "${GOLDEN_DIR}/gnss" \
    "${GOLDEN_OUTPUT_DIR}"

WORK_DIR="${OUTPUT_DIR}/_work"
mkdir -p "${WORK_DIR}"

# Step 1: Download DISP-S1 data
echo "[1/5] Downloading DISP-S1 data..."
cal-disp download disp-s1 \
    --frame-id "${FRAME_ID}" \
    --start    "${START_DATE}" \
    --end      "${END_DATE}" \
    -o         "${GOLDEN_DIR}/disp"

DISP_FILE=$(find "${GOLDEN_DIR}/disp" -name "OPERA_L3_DISP-S1_*.nc" -type f | head -n 1)
if [[ -z "${DISP_FILE}" ]]; then
    echo "ERROR: No DISP file downloaded." >&2; exit 1
fi
echo "  using: $(basename "${DISP_FILE}")"

# Step 2: Download UNR GNSS data
echo "[2/5] Downloading UNR GNSS data..."
cal-disp download unr \
    --frame-id "${FRAME_ID}" \
    -o         "${GOLDEN_DIR}/gnss"

UNR_LOOKUP=$(find "${GOLDEN_DIR}/gnss" -name "grid_latlon_lookup*.txt" -type f | head -n 1)
if [[ -z "${UNR_LOOKUP}" ]]; then
    echo "ERROR: UNR lookup table not found after download." >&2; exit 1
fi

# Step 3: Download tropospheric data (optional)
if [[ "${SKIP_TROPO}" == false ]]; then
    echo "[3/5] Downloading tropospheric data..."
    cal-disp download tropo \
        --input-file "${DISP_FILE}" \
        -o           "${GOLDEN_DIR}/tropo"
else
    echo "[3/5] Skipping tropospheric download (--skip-tropo)."
fi

# Step 4: Link static layers into golden dir
echo "[4/5] Linking static layers..."
ln -sf "$(realpath "${LOS_FILE}")" "${GOLDEN_DIR}/los.tif"
ln -sf "$(realpath "${DEM_FILE}")" "${GOLDEN_DIR}/dem.tif"
cp "${ALGO_PARAMS}" "${GOLDEN_DIR}/algorithm_parameters.yaml"
echo "  los.tif  → $(realpath "${LOS_FILE}")"
echo "  dem.tif  → $(realpath "${DEM_FILE}")"

# Step 5: Generate config and run calibration
echo "[5/5] Running calibration workflow..."
CONFIG_FILE="${WORK_DIR}/runconfig.yaml"

cal-disp config \
    -d  "${DISP_FILE}" \
    -ul "${UNR_LOOKUP}" \
    -ud "${GOLDEN_DIR}/gnss" \
    -uv "${UNR_VERSION}" \
    -ut "${UNR_TYPE}" \
    --los-file  "${GOLDEN_DIR}/los.tif" \
    --dem-file  "${GOLDEN_DIR}/dem.tif" \
    -a  "${GOLDEN_DIR}/algorithm_parameters.yaml" \
    -c  "${CONFIG_FILE}" \
    --frame-id  "${FRAME_ID}" \
    -o  "${GOLDEN_OUTPUT_DIR}" \
    --work-dir  "${WORK_DIR}"

cal-disp run "${CONFIG_FILE}"

# Summary
echo ""
echo "Done."
echo ""
echo "Golden inputs : ${GOLDEN_DIR}"
echo "Golden output : ${GOLDEN_OUTPUT_DIR}"
ls -lh "${GOLDEN_OUTPUT_DIR}"/*.nc 2>/dev/null || true
echo ""
echo "To validate a new run against this reference:"
echo "  cal-disp validate <new_output.nc> ${GOLDEN_OUTPUT_DIR}/<reference.nc>"
echo ""
echo "To run the automated workflow integration tests (no golden data needed):"
echo "  pytest -m integration"
