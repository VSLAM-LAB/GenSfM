#!/bin/bash

delete_if_exists() {
  local folder=$1
  build_folder="${folder}/build"
  if [ -d "$build_folder" ]; then
    rm -rf "$build_folder"
  fi
}

# Check inputs
force_build=false
verbose=false
for input in "$@"
do
    echo "Processing input: $input"
    if [ "$input" = "-f" ]; then
  	force_build=true   
    fi
    if [ "$input" = "-v" ]; then
  	verbose=true   
    fi
done

# Baseline Dir
gensfm_PATH=$(realpath "$0")
gensfm_DIR=$(dirname "$gensfm_PATH")

## Compile colmap
source_folder="${gensfm_DIR}"     
build_folder="$source_folder/build"

if [ "$force_build" = true ]; then
	delete_if_exists ${source_folder}
fi

if [ "$verbose" = true ]; then
        echo "[gensfm-dev][build.sh] Compile gensfm ... "  
	cmake -G Ninja -B $build_folder -S $source_folder\
  -DCMAKE_PREFIX_PATH=$source_folder \
  -DCMAKE_INSTALL_PREFIX=$source_folder \
  -DCMAKE_BUILD_TYPE=Release \
  -DPython_EXECUTABLE=$(which python) \
  -DPython_INCLUDE_DIR=$(python -c "from sysconfig import get_paths; print(get_paths()['include'])") \
  -DPython_LIBRARY=$(find $CONDA_PREFIX -name "libpython3.10.so" | head -n 1) \
  -DBOOST_STATIC=OFF
	ninja install -C $build_folder
else    
	echo "[gensfm-dev][build.sh] Compile gensfm (output disabled) ..."  
	cmake -G Ninja -B $build_folder -S $source_folder\
  -DCMAKE_PREFIX_PATH=$source_folder \
  -DCMAKE_INSTALL_PREFIX=$source_folder \
  -DCMAKE_BUILD_TYPE=Release \
  -DPython_EXECUTABLE=$(which python) \
  -DPython_INCLUDE_DIR=$(python -c "from sysconfig import get_paths; print(get_paths()['include'])") \
  -DPython_LIBRARY=$(find $CONDA_PREFIX -name "libpython3.10.so" | head -n 1) \
  -DBOOST_STATIC=OFF> /dev/null 2>&1
	ninja install -C $build_folder > /dev/null 2>&1
fi


