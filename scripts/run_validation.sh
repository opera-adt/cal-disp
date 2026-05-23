#!/bin/bash
# run_validation.sh
#
# Validate a cal-disp installation against a golden dataset.
#
# Runs the calibration workflow on the golden inputs, then compares the
# output against the committed golden reference using `cal-disp validate`.
#
# Usage
# -----
#   ./scripts/run_validation.sh --golden-dir /path/to/golden_dataset
#
# Options
#   --golden-dir DIR    Root directory of the golden dataset (required).
#                       Must contain golden/ and golden_output/ subdirectories
#                       as produced by create_golden_dataset.py or
#                       build_golden_output.sh.
#   --tolerance FLOAT   Floating-point tolerance passed to cal-disp validate.
#                       Default: 1e-6.
#   --group GROUP       Which product group to validate: main, auxiliary, all.
#                       Default: all.
#   -h, --help          Show this help message.
#
# Exit codes
#   0  Validation passed.
#   1  Validation failed or workflow error.

set -euo pipefail

# Defaults
GOLDEN_DIR=""
TOLERANCE="1e-6"
GROUP="all"

# Argument parsing
usage() {
    grep "^#" "$0" | grep -v "^#!/" | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --golden-dir)  GOLDEN_DIR="$2"; shift 2 ;;
        --tolerance)   TOLERANCE="$2";  shift 2 ;;
        --group)       GROUP="$2";      shift 2 ;;
        -h|--help)     usage ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate inputs
if [[ -z "${GOLDEN_DIR}" ]]; then
    echo "ERROR: --golden-dir is required." >&2
    echo "Usage: $0 --golden-dir /path/to/golden_dataset" >&2
    exit 1
fi

GOLDEN_DIR="$(realpath "${GOLDEN_DIR}")"
INPUT_DIR="${GOLDEN_DIR}/input_data"
CONFIGS_DIR="${GOLDEN_DIR}/configs"
REFERENCE_DIR="${GOLDEN_DIR}/golden_output"
OUTPUT_DIR="${GOLDEN_DIR}/output"

for d in "${INPUT_DIR}" "${CONFIGS_DIR}" "${REFERENCE_DIR}"; do
    if [[ ! -d "${d}" ]]; then
        echo "ERROR: Expected directory not found: ${d}" >&2
        echo "Make sure --golden-dir points to the root produced by" >&2
        echo "create_golden_dataset.py or build_golden_output.sh." >&2
        exit 1
    fi
done

DISP_FILE=$(find "${INPUT_DIR}/disp" -name "OPERA_L3_DISP-S1_*.nc" -type f | head -n 1)
LOS_FILE=$(find "${INPUT_DIR}/static_input" -name "*line_of_sight_enu.tif" -type f | head -n 1)
DEM_FILE=$(find "${INPUT_DIR}/static_input" -name "*_dem.tif" -type f | head -n 1)
LOOKUP_FILE="${INPUT_DIR}/gnss/grid_latlon_lookup.txt"
GNSS_DIR="${INPUT_DIR}/gnss"
ALGO_FILE="${CONFIGS_DIR}/algorithm_parameters.yaml"
REFERENCE_NC=$(find "${REFERENCE_DIR}" -name "OPERA_L4_DISP-CAL-S1_*.nc" -type f | head -n 1)

# Tropo files — sorted by filename gives chronological order (ref then sec)
mapfile -t TROPO_FILES < <(find "${INPUT_DIR}/tropo" -name "OPERA_L4_TROPO-ZENITH_*.nc" -type f | sort)
REF_TROPO_FILE="${TROPO_FILES[0]:-}"
SEC_TROPO_FILE="${TROPO_FILES[1]:-}"

for f in "${DISP_FILE}" "${LOS_FILE}" "${DEM_FILE}" "${LOOKUP_FILE}" "${ALGO_FILE}" \
         "${REFERENCE_NC}" "${REF_TROPO_FILE}" "${SEC_TROPO_FILE}"; do
    if [[ -z "${f}" || ! -f "${f}" ]]; then
        echo "ERROR: Required file not found: ${f:-<no match>}" >&2
        exit 1
    fi
done

# Parse frame-id and unr-type from the algorithm parameters file
# Read unr grid type (constant/variable) if present, else default to constant
UNR_TYPE="constant"
if grep -q "grid_type" "${ALGO_FILE}" 2>/dev/null; then
    UNR_TYPE=$(grep "grid_type" "${ALGO_FILE}" | awk '{print $2}' | tr -d '"' | head -n 1)
fi

# Parse frame-id from DISP filename  (e.g. F36540 -> 36540)
FRAME_ID=$(basename "${DISP_FILE}" | grep -oP '(?<=_F)\d+')

# Working and output directories (use the delivered layout)
WORK_DIR="${GOLDEN_DIR}/output/_work"
TEST_OUTPUT_DIR="${GOLDEN_DIR}/output"
CONFIG_FILE="${WORK_DIR}/runconfig.yaml"

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}" "${TEST_OUTPUT_DIR}"

# Generate config
echo "=== Running validation ==="
echo "  golden inputs : ${INPUT_DIR}"
echo "  reference     : $(basename "${REFERENCE_NC}")"
echo "  frame         : ${FRAME_ID}"
echo "  unr-type      : ${UNR_TYPE}"
echo ""

echo "[1/3] Generating configuration..."
cal-disp config \
    -d  "${DISP_FILE}" \
    -ul "${LOOKUP_FILE}" \
    -ud "${GNSS_DIR}" \
    -uv "0.3" \
    -ut "${UNR_TYPE}" \
    --los-file  "${LOS_FILE}" \
    --dem-file  "${DEM_FILE}" \
    --ref-tropo-files "${REF_TROPO_FILE}" \
    --sec-tropo-files "${SEC_TROPO_FILE}" \
    -a  "${ALGO_FILE}" \
    --frame-id "${FRAME_ID}" \
    -o  "${TEST_OUTPUT_DIR}" \
    --work-dir "${WORK_DIR}" \
    -c  "${CONFIG_FILE}"

# Run calibration
echo "[2/3] Running calibration..."
cal-disp run "${CONFIG_FILE}"

TEST_NC=$(find "${TEST_OUTPUT_DIR}" -name "OPERA_L4_DISP-CAL-S1_*.nc" -type f | head -n 1)
if [[ -z "${TEST_NC}" ]]; then
    echo "ERROR: Calibration produced no output NetCDF." >&2
    exit 1
fi

# Validate
echo "[3/3] Validating output..."
if cal-disp validate \
    --tolerance "${TOLERANCE}" \
    --group     "${GROUP}" \
    "${REFERENCE_NC}" \
    "${TEST_NC}"; then
    echo ""
    echo "✓ Validation passed."
    EXIT_CODE=0
else
    echo ""
    echo "✗ Validation FAILED. See above for details." >&2
    EXIT_CODE=1
fi

# Clean up scratch directory only; leave output/ intact for inspection
rm -rf "${WORK_DIR}"

exit "${EXIT_CODE}"
