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
cp ../data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_line_of_sight_enu.tif \
   ./input_data/static_input/

cp ../data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_dem.tif \
   ./input_data/static_input/

# Step 6: Create configuration file
echo "Creating configuration..."
CONFIG_FILE="configs/runconfig.yaml"

cal-disp config \
    -d "${DISP_FILE}" \
    -cl ./input_data/unr/grid_latlon_lookup_v0.2.txt \
    -cd ./input_data/unr \
    --los-file ./input_data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_line_of_sight_enu.tif \
    --dem-file ./input_data/static_input/OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_dem.tif \
    -a configs/algorithm_parameters.yaml \
    -c "${CONFIG_FILE}"

# Step 7: Run calibration
echo "Running calibration..."
cal-disp run "${CONFIG_FILE}"

# Step 8: Move outputs to golden_output directory
echo "Organizing outputs..."
if [ -d "./outputs" ]; then
    mv ./outputs/*.nc ./golden_output/ 2>/dev/null || true
    mv ./outputs/*.png ./golden_output/ 2>/dev/null || true
fi

echo "Done! Results are in golden_output/"
