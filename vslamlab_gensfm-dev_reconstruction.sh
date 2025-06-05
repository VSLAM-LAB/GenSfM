#!/bin/bash

# Default values
matcher_type="exhaustive"
use_gpu="1"
verbose="0"
settings_yaml=""
sequence_path=""
exp_folder=""
exp_id=""
calibration_yaml=""
rgb_txt=""

# Function to split key-value pairs and assign them to variables
split_and_assign() {
  local input=$1
  local key=$(echo $input | cut -d':' -f1)
  local value=$(echo $input | cut -d':' -f2-)
  eval $key=$value
}

# Read Inputs
for ((i=1; i<=$#; i++)); do
    split_and_assign "${!i}"
done

exp_id=$(printf "%05d" ${exp_id})

echo -e "\n================= Experiment Configuration ================="
echo "  Sequence Path     : $sequence_path"
echo "  Experiment Folder : $exp_folder"
echo "  Experiment ID     : $exp_id"
echo "  Verbose           : $verbose"
echo "  Matcher Type      : $matcher_type"
echo "  Use GPU           : $use_gpu"
echo "  Settings YAML     : $settings_yaml"
echo "  Calibration YAML  : $calibration_yaml"
echo "  RGB TXT           : $rgb_txt"
echo "============================================================"

# Create folder to save colmap files
exp_folder_colmap="${exp_folder}/colmap_${exp_id}"
rm -rf "$exp_folder_colmap"
mkdir "$exp_folder_colmap"

# Run COLMAP scripts for matching and mapping
export QT_QPA_PLATFORM_PLUGIN_PATH="$CONDA_PREFIX/plugins/platforms"
gensfm_args="$sequence_path $exp_folder $exp_id $settings_yaml $calibration_yaml $rgb_txt"

./vslamlab_gensfm-dev_matcher.sh $gensfm_args $matcher_type $use_gpu
./vslamlab_gensfm-dev_mapper.sh $gensfm_args 

# Visualization with colmap gui
if [ "$verbose" -eq 1 ]; then
  exp_folder_colmap="${exp_folder}/colmap_${exp_id}"
  rgb_path="${sequence_path}/$(awk '{print $2}' "${rgb_txt}" | awk -F'/' 'NR==1 {print $1}')"
  database="${exp_folder_colmap}/colmap_database.db"
  ./build/src/exe/gen_colmap gui --import_path "${exp_folder_colmap}/0" --database_path ${database} --image_path ${rgb_path}
fi

# # Remove colmap data
# rm -rf ${exp_folder_colmap}


