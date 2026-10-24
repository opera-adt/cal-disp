#!/bin/bash
# build_golden_output.sh
#
# Build the DISP-CAL golden dataset (frame F08882, pair 2022-01-11 / 2022-07-22)
# in the delivery layout, ready for run_validation.sh and the delivered commands.
#
# Requirements
# ------------
#   - cal-disp installed with the download extra (pip install "cal-disp[download]")
#   - curl, and Earthdata login credentials in ~/.netrc
#
# Usage
# -----
#   scripts/build_golden_output.sh [--output-dir DIR] [--skip-download]
#
#   --output-dir DIR   Dataset folder to create (default: ./golden_data_disp_cal).
#   --skip-download    Use the inputs already in DIR/input_data.
#
# Output layout (DIR/)
# --------------------
#   configs/        algorithm_parameters.yaml, runconfig.yaml (paths relative
#                   to the parent of DIR, e.g. the folder mounted at /home/work)
#   input_data/     disp/, static_input/, tropo/, gnss/,
#                   static_ancillary/ (Venti frame table, area and event databases)
#   golden_output/  golden product (.nc, .png)
#   output/         empty; test products from run_validation.sh or cal-disp run
#
# The golden product is only reproducible to the validation tolerance in the
# same software environment. For a delivery, run this script inside the
# delivered Docker image, e.g. from the folder that holds cal-disp/:
#   docker run --rm --user $(id -u):$(id -g) -v $PWD:/home/work -w /home/work \
#       -v ~/.netrc:/home/conda/.netrc:ro <image> \
#       bash cal-disp/scripts/build_golden_output.sh
#
# Verify
# ------
#   scripts/run_validation.sh --golden-dir DIR
# or, as delivered (from the parent of DIR):
#   cal-disp run DIR/configs/runconfig.yaml
#   cal-disp validate DIR/golden_output/<golden>.nc DIR/output/<test>.nc

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Golden case
FRAME_ID=8882
DISP_NAME="OPERA_L3_DISP-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_20251027T005420Z"
STATIC_NAME="OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0"
REF_DATE="20220111"
SEC_DATE="20220722"
UNR_VERSION="0.3"
UNR_TYPE="constant"
ASF_URL="https://cumulus.asf.earthdatacloud.nasa.gov/OPERA"

OUTPUT_DIR="golden_data_disp_cal"
SKIP_DOWNLOAD=false

usage() {
    grep "^#" "$0" | grep -v "^#!/" | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=true; shift ;;
        -h|--help)       usage ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Work from the parent of the dataset folder, so runconfig paths are relative
mkdir -p "${OUTPUT_DIR}"
cd "$(dirname "$(realpath "${OUTPUT_DIR}")")"
D="$(basename "${OUTPUT_DIR}")"

if compgen -G "${D}/golden_output/OPERA_L4_DISP-CAL-S1_*.nc" > /dev/null; then
    echo "ERROR: ${D}/golden_output already has a golden product; remove it first." >&2
    exit 1
fi
mkdir -p "${D}"/configs "${D}"/input_data/{disp,static_input,tropo,gnss,static_ancillary} \
    "${D}"/golden_output "${D}"/output

DISP_FILE="${D}/input_data/disp/${DISP_NAME}.nc"
LOS_FILE="${D}/input_data/static_input/${STATIC_NAME}_line_of_sight_enu.tif"
DEM_FILE="${D}/input_data/static_input/${STATIC_NAME}_dem.tif"
GNSS_DIR="${D}/input_data/gnss"
TROPO_DIR="${D}/input_data/tropo"

echo "=== Building golden dataset in $(pwd)/${D}"

# 1. Inputs
if [[ "${SKIP_DOWNLOAD}" == false ]]; then
    fetch() {  # Earthdata download (credentials from ~/.netrc)
        [[ -s "$2" ]] && return 0
        curl --fail --silent --show-error --netrc --location \
            --cookie ~/.edl_cookies --cookie-jar ~/.edl_cookies -o "$2" "$1"
    }
    echo "[1/4] Downloading inputs..."
    fetch "${ASF_URL}/OPERA_L3_DISP-S1_V1/${DISP_NAME}/${DISP_NAME}.nc" "${DISP_FILE}"
    fetch "${ASF_URL}/OPERA_L3_DISP-S1-STATIC_V1/${STATIC_NAME}/$(basename "${LOS_FILE}")" "${LOS_FILE}"
    fetch "${ASF_URL}/OPERA_L3_DISP-S1-STATIC_V1/${STATIC_NAME}/$(basename "${DEM_FILE}")" "${DEM_FILE}"
    if ! compgen -G "${TROPO_DIR}/OPERA_L4_TROPO-ZENITH_*.nc" > /dev/null; then
        cal-disp download tropo -i "${DISP_FILE}" -o "${TROPO_DIR}"
    fi
    cal-disp download unr --frame-id "${FRAME_ID}" --grid-type "${UNR_TYPE}" -o "${GNSS_DIR}"
else
    echo "[1/4] Using existing inputs (--skip-download)."
fi

REF_TROPO=$(find "${TROPO_DIR}" -name "OPERA_L4_TROPO-ZENITH_${REF_DATE}T*.nc" | sort | head -n 1)
SEC_TROPO=$(find "${TROPO_DIR}" -name "OPERA_L4_TROPO-ZENITH_${SEC_DATE}T*.nc" | sort | head -n 1)
LOOKUP="${GNSS_DIR}/grid_latlon_lookup_v${UNR_VERSION}.txt"
for f in "${DISP_FILE}" "${LOS_FILE}" "${DEM_FILE}" "${LOOKUP}" "${REF_TROPO}" "${SEC_TROPO}"; do
    if [[ -z "${f}" || ! -s "${f}" ]]; then
        echo "ERROR: missing input: ${f:-TROPO file for ${REF_DATE}/${SEC_DATE}}" >&2
        exit 1
    fi
done

# 2. Configs
echo "[2/4] Writing configs..."
# v0.5 configuration (Venti plan T49; VnV sign-off 2026-10-15): the repo's
# configs/algorithm_parameters_v05.yaml, and the static ancillaries production
# uses: the materialized Venti frame table and the area and event databases
# of the installed Venti. Gamma 0.3 (decision 0001) is reproduced with
# configs/algorithm_parameters_gamma.yaml.
cp "${SCRIPT_DIR}/../configs/algorithm_parameters_v05.yaml" "${D}/configs/algorithm_parameters.yaml"
ANC="${D}/input_data/static_ancillary"
python - "${ANC}" <<'PY'
import shutil, sys
from pathlib import Path

import venti
from venti.frames import load_frame_table

out = Path(sys.argv[1])
load_frame_table().write(out / "frame_parameters.json", materialize=True)
for name in ("defo_area_db_v1.geojson", "event_db_v1.geojson"):
    shutil.copy(Path(venti.__file__).parent / "data" / name, out / name)
PY

cal-disp config \
    -d  "${DISP_FILE}" \
    -f  "${FRAME_ID}" \
    -ul "${LOOKUP}" \
    -ud "${GNSS_DIR}" \
    -uv "${UNR_VERSION}" \
    -ut "${UNR_TYPE}" \
    --los-file "${LOS_FILE}" \
    --dem-file "${DEM_FILE}" \
    --ref-tropo-files "${REF_TROPO}" \
    --sec-tropo-files "${SEC_TROPO}" \
    -a  "${D}/configs/algorithm_parameters.yaml" \
    --algorithm-overrides "${ANC}/frame_parameters.json" \
    --defo-area-db "${ANC}/defo_area_db_v1.geojson" \
    --event-db "${ANC}/event_db_v1.geojson" \
    -o  "${D}/output" \
    --work-dir "${D}/output/_work" \
    --keep-relative \
    -c  "${D}/configs/runconfig.yaml" > /dev/null
rm -rf "${D}/output/_work"

# 3. Golden run
echo "[3/4] Running calibration..."
cal-disp run "${D}/configs/runconfig.yaml"
mv "${D}"/output/OPERA_L4_DISP-CAL-S1_* "${D}/golden_output/"
rm -rf "${D}/output/_work" "${D}/output/cal_disp.log"

echo "[4/4] Done."
echo "  Golden product: $(ls "${D}"/golden_output/*.nc)"
echo ""
echo "Verify with:  scripts/run_validation.sh --golden-dir $(pwd)/${D}"
echo "or, from $(pwd):"
echo "  cal-disp run ${D}/configs/runconfig.yaml"
echo "  cal-disp validate ${D}/golden_output/<golden>.nc ${D}/output/<test>.nc"
