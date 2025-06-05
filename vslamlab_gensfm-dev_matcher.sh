#!/bin/bash
echo ""
echo "Executing colmap_matcher.sh ..."

sequence_path="$1"
exp_folder="$2" 
exp_id="$3" 
settings_yaml="$4"
calibration_yaml="$5"
rgb_txt="$6"
matcher_type="$7"
use_gpu="$8"

exp_folder_colmap="${exp_folder}/colmap_${exp_id}"
rgb_path="${sequence_path}/$(awk '{print $2}' "${rgb_txt}" | awk -F'/' 'NR==1 {print $1}')"

# Reading settings from yaml file
# feature_extractor_SiftExtraction_num_octaves=$(yq '.feature_extractor.SiftExtraction_num_octaves // 4.0' $settings_yaml)
# feature_extractor_SiftExtraction_octave_resolution=$(yq '.feature_extractor.SiftExtraction_octave_resolution // 3.0' $settings_yaml)
# feature_extractor_SiftExtraction_peak_threshold=$(yq '.feature_extractor.SiftExtraction_peak_threshold // 0.0066666666666666671' $settings_yaml)
# feature_extractor_SiftExtraction_edge_threshold=$(yq '.feature_extractor.SiftExtraction_edge_threshold // 10.0' $settings_yaml)
# feature_extractor_SiftExtraction_dsp_min_scale=$(yq '.feature_extractor.SiftExtraction_dsp_min_scale // 0.1666666666666666' $settings_yaml)
# feature_extractor_SiftExtraction_dsp_max_scale=$(yq '.feature_extractor.SiftExtraction_dsp_max_scale // 3.0' $settings_yaml)
# feature_extractor_SiftExtraction_dsp_num_scales=$(yq '.feature_extractor.SiftExtraction_dsp_num_scales // 10.0' $settings_yaml)

# matcher_SiftMatching_max_ratio=$(yq '.matcher.SiftMatching_max_ratio // 0.80000000000000004' $settings_yaml)
# matcher_SiftMatching_max_distance=$(yq '.matcher.SiftMatching_max_distance // 0.69999999999999996' $settings_yaml)
# matcher_TwoViewGeometry_min_num_inliers=$(yq '.matcher.TwoViewGeometry_min_num_inliers // 15.0' $settings_yaml)
# matcher_TwoViewGeometry_max_error=$(yq '.matcher.TwoViewGeometry_max_error // 4.0' $settings_yaml)
# matcher_TwoViewGeometry_confidence=$(yq '.matcher.TwoViewGeometry_confidence // 0.999' $settings_yaml)
# matcher_TwoViewGeometry_min_inlier_ratio=$(yq '.matcher.TwoViewGeometry_min_inlier_ratio // 0.25' $settings_yaml)
# matcher_SequentialMatching_overlap=$(yq '.matcher.SequentialMatching_overlap // 10.0' $settings_yaml)
# matcher_SequentialMatching_quadratic_overlap=$(yq '.matcher.SequentialMatching_quadratic_overlap // 1.0' $settings_yaml)
# matcher_ExhaustiveMatching_block_size=$(yq '.matcher.ExhaustiveMatching_block_size // 50.0' $settings_yaml)

# Create colmap image list
colmap_image_list="${exp_folder_colmap}/colmap_image_list.txt"
awk '{split($2, arr, "/"); print arr[2]}' "$rgb_txt" > "$colmap_image_list"

# Create Colmap Database
pwd
database="${exp_folder_colmap}/colmap_database.db"
rm -rf ${database}
./build/src/exe/gen_colmap database_creator --database_path ${database}

# Feature extractor
echo "    colmap feature_extractor ..."
echo "        camera model : IMPLICIT_DISTORTION"
./build/src/exe/gen_colmap feature_extractor \
  --database_path ${database} \
  --image_path ${rgb_path} \
  --image_list_path ${colmap_image_list} \
  --ImageReader.camera_model IMPLICIT_DISTORTION \
  --ImageReader.single_camera_per_folder 1

# Exhaustive Feature Matcher
if [ "${matcher_type}" == "exhaustive" ]
then
	echo "    colmap exhaustive_matcher ..."
  ./build/src/exe/gen_colmap exhaustive_matcher \
  --database_path ${database} \
  --SiftMatching.multiple_models=1
fi

# # Sequential Feature Matcher
# if [ "${matcher_type}" == "sequential" ]
# then
#   num_rgb=$(wc -l < ${rgb_txt})

#   # Pick vocabulary tree based on the number of images
#   vocabulary_tree="Baselines/colmap/vocab_tree_flickr100K_words32K.bin"
#   if [ "$num_rgb" -gt 1000 ]; then
#     vocabulary_tree="Baselines/colmap/vocab_tree_flickr100K_words256K.bin"
#   fi
#   if [ "$num_rgb" -gt 10000 ]; then
#     vocabulary_tree="Baselines/colmap/vocab_tree_flickr100K_words1M.bin"
#   fi

#   echo "    colmap sequential_matcher ..."
#   echo "        Vocabulary Tree: $vocabulary_tree"
#   ./Baselines/GenSfM-DEV/build/src/exe/gen_colmap sequential_matcher \
#     --database_path "${database}" \
#     --SequentialMatching.loop_detection 1 \
#     --SequentialMatching.vocab_tree_path ${vocabulary_tree} \
#     --SiftMatching.use_gpu "${use_gpu}" \
#     --SiftMatching.max_ratio "${matcher_SiftMatching_max_ratio}" \
#     --SiftMatching.max_distance "${matcher_SiftMatching_max_distance}" \
#     --TwoViewGeometry.min_num_inliers "${matcher_TwoViewGeometry_min_num_inliers}" \
#     --TwoViewGeometry.max_error "${matcher_TwoViewGeometry_max_error}" \
#     --TwoViewGeometry.confidence "${matcher_TwoViewGeometry_confidence}" \
#     --TwoViewGeometry.min_inlier_ratio "${matcher_TwoViewGeometry_min_inlier_ratio}" \
#     --SequentialMatching.overlap "${matcher_SequentialMatching_overlap}" \
#     --SequentialMatching.quadratic_overlap "${matcher_SequentialMatching_quadratic_overlap}"
# fi