#!/bin/bash
set -e  # Exit on error

# Variables
FRAME_ID="8882"
START_DATE="2022-07-01"
END_DATE="2022-08-01"

# Step 1: Create directory structure
echo "Creating directories..."
mkdir -p delivery_data_cal_disp/configs
mkdir -p delivery_data_cal_disp/input_data/{disp,static_input,tropo,unr}
mkdir -p delivery_data_cal_disp/golden_output

# Step 2: Navigate to working directory
cd delivery_data_cal_disp || exit 1

# Get absolute path for working directory
WORK_DIR=$(pwd)
echo "Working directory: ${WORK_DIR}"

# Step 3: Download DISP-S1 data
echo "Downloading DISP-S1 data..."
cal-disp download disp-s1 \
    --frame-id "${FRAME_ID}" \
    -s "${START_DATE}" \
    -e "${END_DATE}" \
    -o ./input_data/disp

# Find the downloaded DISP file
DISP_FILE=$(find ./input_data/disp -name "OPERA_L3_DISP-S1_*.nc" -type f | head -n 1)
echo "Using DISP file: ${DISP_FILE}"

# Convert to absolute path
DISP_FILE_ABS="${WORK_DIR}/${DISP_FILE#./}"
echo "Absolute DISP path: ${DISP_FILE_ABS}"

# Step 4: Download tropospheric correction data
echo "Downloading tropospheric data..."
cal-disp download tropo \
    -i "${DISP_FILE}" \
    -o ./input_data/tropo

# Step 5: Download UNR GPS data
echo "Downloading UNR data..."
cal-disp download unr \
    --frame-id "${FRAME_ID}" \
    -o ./input_data/unr

# Step 5.5: Copy static layers to input directory
echo "Copying static layers..."
cp ../../data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_line_of_sight_enu.tif \
   ./input_data/static_input/

cp ../../data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_dem.tif \
   ./input_data/static_input/

cp ../../configs/algorithm_parameters.yaml \
   ./configs/

# Step 6: Create configuration file
echo "Creating configuration..."
CONFIG_FILE="runconfig.yaml"
echo "Config file will be created at: ${CONFIG_FILE}"

cal-disp config \
    -d "${DISP_FILE_ABS}" \
    -ul "${WORK_DIR}/input_data/unr/grid_latlon_lookup_v0.2.txt" \
    -ud "${WORK_DIR}/input_data/unr" \
    -uv "0.2" \
    -ut "variable" \
    --los-file "${WORK_DIR}/input_data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_line_of_sight_enu.tif" \
    --dem-file "${WORK_DIR}/input_data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_dem.tif" \
    -a "${WORK_DIR}/configs/algorithm_parameters.yaml" \
    -c "${CONFIG_FILE}" \
    --frame-id "${FRAME_ID}" \
    -o "${WORK_DIR}/golden_output" \
    --work-dir "${WORK_DIR}"

cp "${WORK_DIR}"/"${CONFIG_FILE}" \
   ./configs/

# Verify config file was created
if [ ! -f "${CONFIG_FILE}" ]; then
    echo "ERROR: Config file was not created at ${WORK_DIR}/${CONFIG_FILE}"
    echo "Checking for config file in current directory..."
    find . -name "runconfig.yaml" -type f
    exit 1
fi

echo "✓ Config file created: ${CONFIG_FILE}"

# Step 7: Run calibration
echo "Running calibration..."
cal-disp run "${CONFIG_FILE}"

# Step 8: Move outputs to golden_output directory (if needed)
echo "Organizing outputs..."
if [ -d "${WORK_DIR}/golden_output" ]; then
    mv ./outputs/*.nc ./golden_output/ 2>/dev/null || true
    mv ./outputs/*.png ./golden_output/ 2>/dev/null || true
    echo "Moved outputs from ./outputs/ to ./golden_output/"
fi

echo ""
echo "Done! Results are in golden_output/"
echo "Location: ${WORK_DIR}/golden_output/"
ls -lh "${WORK_DIR}/golden_output/"
