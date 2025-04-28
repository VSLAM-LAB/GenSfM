
Structure-from-Motion with a Non-Parametric Camera Model (CVPR 2025 Highlight)
======
This is an extension of the incremental Structure-from-Motion framework [COLMAP](https://github.com/colmap/colmap), allowing for the usage of Implicit Distortion camera model, which features a non-parametric generic representation of the calibration map and an adaptive partial calibration procedure. This is a cleaned up re-implementation of the original code used in the paper

```
Structure-from-Motion with a Non-Parametric Camera Model
Yihan Wang*, Linfei Pan*, Marc Pollefeys, Viktor Larsson
IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR) 2025
```
and if you have any problems running the code or find any bugs, please create an issue.

## Getting Started
Please build GenSfM using the following commands:
```shell
mkdir build
cd build
cmake .. -GNinja
ninja
```
## End-to-End Example
We provide two example datasets that could be downloaded from [here](https://drive.google.com/drive/folders/1ut-anuoDhFHkG3e54uN05fR_g_l20NFb?usp=sharing). Please download them and put them under `data` folder. Navigate to `./build/src/exe` to run GenSfM.
### Run from database
If a COLMAP database already exists, GenSfM can directly use it to perform incremental mapping by: 

```shell
./gen_colmap mapper \
            --database_path  ./data/CAB_facade_mixed/full_mix_centered.db \
            --image_path  ./data/CAB_facade_mixed/images \
            --output_path ./data/CAB_facade_mixed/sparse \
# Add the following commands if you wish to monitor the middle stages of the reconstruction
            --Mapper.snapshot_path ./data/CAB_facade_mixed/sparse/snapshot \
            --Mapper.snapshot_images_freq 1
```
### Run from images
If there is not a COLMAP database yet, you need to establish it first. 

```shell
./gen_colmap feature_extractor \
            --image_path ./data/Fisheye_grossmunster/images \
            --database_path ./data/Fisheye_grossmunster/database.db \ 
# GenSFM could be ran by choosing IMPLICIT_DISTORTION camera model            
            --ImageReader.camera_model IMPLICIT_DISTORTION \
# Set the following argument if images are captured by the same camera per folder
            --ImageReader.single_camera_per_folder 1


./gen_colmap exhaustive_matcher \
            --database_path ./data/Fisheye_grossmunster/database.db \ 
            --SiftMatching.multiple_models=1

./gen_colmap mapper \
            --database_path  ./data/Fisheye_grossmunster/database.db \
            --image_path  ./data/Fisheye_grossmunster/images \
            --output_path ./data/Fisheye_grossmunster/sparse \
            --Mapper.snapshot_path ./data/Fisheye_grossmunster/sparse/snapshot \
            --Mapper.snapshot_images_freq 1
```
### Catadioptric feature matching
We provide experimental scripts under experimental_scripts/cata_feature_matching to facilitate the feature matching for catadioptric images. To perform this, first run `split_images.m` to psuedo-undistort the catadioptric images into 8 rectangular views, and then use the following steps to perform feature merging and matching:
```bash
./gen_colmap feature_extractor --database_path database_split.db --image_path images_split/

./gen_colmap database_creator --database_path database.db

python merge_features.py database_split.db database.db 4288,2848,2190,1463,350,1250,796 # Change to the parameters obtained from `split_images.m`

./gen_colmap exhaustive_matcher --database_path database.db --SiftMatching.multiple_models=1

```
### Visualize the reconstructions
The results are written out in the COLMAP sparse reconstruction format. Please refer to [COLMAP](https://github.com/colmap/colmap) for more details. The reconstruction can be visualized using the GUI, for example:
```bash
./gen_colmap gui \
        --database_path ./data/Fisheye_grossmunster/database.db    \
        --image_path ./data/Fisheye_grossmunster/images     \
        --import_path ./data/Fisheye_grossmunster/sparse/0
```
### Undistorting images using estimated calibration
We also provide scripts for undistorting fisheye or catadioptric images using the estimated calibration. Please first follow [Implicit Distortion](https://github.com/cvg/implicit_dist.git) to install the dependencies. 

Then, store the calibration spline control points {$\theta_i$, $r_i$} estimated during reconstruction to a txt file as the examples in `./experimental_scripts/undistortion/test_data`. Then see the notebook `./experimental_scripts/undistortion/test_undistortion.ipynb` to run the undistortion.

## Acknowledgement
This pipeline builds on COLMAP, RadialSfM, and Implicit_dst, for which you should cite:

    @inproceedings{schoenberger2016sfm,
        author={Sch\"{o}nberger, Johannes Lutz and Frahm, Jan-Michael},
        title={Structure-from-Motion Revisited},
        booktitle={Conference on Computer Vision and Pattern Recognition (CVPR)},
        year={2016},
    }
```
@inproceedings{schoenberger2016mvs,
        author={Sch\"{o}nberger, Johannes Lutz and Zheng, Enliang and Pollefeys, Marc and Frahm, Jan-Michael},
        title={Pixelwise View Selection for Unstructured Multi-View Stereo},
        booktitle={European Conference on Computer Vision (ECCV)},
        year={2016},
}
```
```
@inproceedings{larsson2020calibration,
        title={Calibration-free structure-from-motion with calibrated radial trifocal tensors},
        author={Larsson, Viktor and Zobernig, Nicolas and Taskin, Kasim and Pollefeys, Marc},
        booktitle={European Conference on Computer Vision (ECCV)},
        year={2020},
}
```
```
@inproceedings{pan2022camera,
        title={Camera pose estimation using implicit distortion models},
        author={Pan, Linfei and Pollefeys, Marc and Larsson, Viktor},
        booktitle={Conference on Computer Vision and Pattern Recognition (CVPR)}, 
        year={2022}
}