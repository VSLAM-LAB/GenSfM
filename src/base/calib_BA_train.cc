#ifndef CALIB_BA_TAIN_CC_
#define CALIB_BA_TAIN_CC_
#include <iostream>
#include <vector>
#include <cassert>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>

#include <filesystem>
// #include "base/projection.h"
// #include "estimators/triangulation.h"
// #include "estimators/implicit_camera_pose.h"
// #include "estimators/implicit_cost_matrix.h"
// #include "estimators/implicit_pose_refinement.h"
#include "util/math.h"
#include "util/types.h"
// #include "util/misc.h"
#include "spline.h"
#include "controllers/bundle_adjustment.h"
#include "optim/bundle_adjustment.h"
#include "base/camera_models.h"
#include "base/correspondence_graph.h"
#include "base/projection.h"
#include "base/reconstruction.h"



int main(){
    

    // the base path to the files
    std::string base_path = "/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/reconstructions/";
    std::string output_base = "/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/spline/";
    // list of files to read
    for (const auto& entry : std::filesystem::directory_iterator(base_path)){
        // list the subdirectories
        std::string subdirectory = entry.path();
        std::cout << "BA refining for: " << subdirectory << std::endl;
        // read reconstruction from subdirectory
        colmap::Reconstruction reconstruction;
        reconstruction.ReadText(subdirectory);
        std::cout << "reconstructions reading successful" << std::endl;
        // run bundle adjustment
        std::cout << "Number of registered images: " << reconstruction.NumRegImages() << std::endl;
        std::cout << "Number of 3D points: " << reconstruction.NumPoints3D() << std::endl;
        
        colmap::BundleAdjustmentOptions ba_options;
        ba_options.refine_focal_length = false;
        ba_options.refine_principal_point = false;
        ba_options.refine_extra_params = true;
        ba_options.refine_extrinsics = false;
        // Configure bundle adjustment.
        colmap::BundleAdjustmentConfig ba_config;
        // configure bundle adjustment config

        for (const colmap::image_t image_id : reconstruction.RegImageIds()) {
            ba_config.AddImage(image_id);
        }
        std::cout << "images added" << std::endl;
        // set all images constant in the bundle adjustment config
        for (const colmap::image_t image_id : reconstruction.RegImageIds()) {
            ba_config.SetConstantPose(image_id);
        }
        // set all points constant in the bundle adjustment config
        for (const colmap::point3D_t point3D_id : reconstruction.Point3DIds()) {
            ba_config.AddConstantPoint(point3D_id);
        }
        std::cout << "points added" << std::endl;
        // run bundle adjustment
        colmap::BundleAdjuster bundle_adjuster(ba_options, ba_config);
        bundle_adjuster.Solve(&reconstruction);
        // write the spline to the output directory
        std::string output_path = output_base + subdirectory.substr(subdirectory.find_last_of("/")+1) + "_spline.txt";


    }
}

#endif