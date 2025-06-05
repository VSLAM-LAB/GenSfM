#!/bin/bash
echo "Executing colmap_mapper.sh ..."

sequence_path="$1" 
exp_folder="$2" 
exp_id="$3" 
settings_yaml="$4"
calibration_yaml="$5"
rgb_txt="$6"

exp_folder_colmap="${exp_folder}/colmap_${exp_id}"
rgb_path="${sequence_path}/$(awk '{print $2}' "${rgb_txt}" | awk -F'/' 'NR==1 {print $1}')"

calibration_model=$(grep -oP '(?<=Camera\.model:\s)[\w]+' "$calibration_yaml")
echo "        camera model : $calibration_model"

# Reading settings from yaml file
mapper_snapshot_images_freq=$(yq '.mapper.Mapper_snapshot_images_freq // 1' $settings_yaml)

echo "    colmap mapper ..."
database="${exp_folder_colmap}/colmap_database.db"

./build/src/exe/gen_colmap mapper \
  --database_path ${database} \
  --image_path ${rgb_path} \
  --output_path ${exp_folder_colmap} #\
  #--Mapper.snapshot_path "${exp_folder_colmap}/snapshot" \
  #--Mapper.snapshot_images_freq ${mapper_snapshot_images_freq}

echo "    colmap model_converter ..."
./build/src/exe/gen_colmap model_converter \
	--input_path ${exp_folder_colmap}/0 --output_path ${exp_folder_colmap} --output_type TXT



