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

#include "util/misc.h"


int main(int argc, char*argv[]){
    bool opt_extra_params = true;
    if (argc > 1) {
        opt_extra_params = std::stoi(argv[1]);
    }

    std::cout << "opt_extra_params: " << opt_extra_params << std::endl; 

    

    // the base path to the files
    std::string base_path = "/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/reconstructions/";
    std::string output_base = "/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/spline/";
    std::string output_base_pose = "/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/poses/";

    
    
    if (opt_extra_params == false) {
        base_path = "/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/reconstructions_test/";
    }
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
        ba_options.refine_extra_params = opt_extra_params;
        ba_options.refine_extrinsics = true;
        // Configure bundle adjustment.
        colmap::BundleAdjustmentConfig ba_config;
        // configure bundle adjustment config

        for (const colmap::image_t image_id : reconstruction.RegImageIds()) {
            ba_config.AddImage(image_id);
        }
        std::cout << "images added" << std::endl;
        // set all images constant in the bundle adjustment config
        // for (const colmap::image_t image_id : reconstruction.RegImageIds()) {
        //     ba_config.SetConstantPose(image_id);
        // }
        // ba_config.SetConstantPose(reg)
        // set all points constant in the bundle adjustment config
        for (const colmap::point3D_t point3D_id : reconstruction.Point3DIds()) {
            ba_config.AddConstantPoint(point3D_id);
        }
        for (auto &[camera_id, const_camera] : reconstruction.Cameras()) {
            colmap::Camera &camera = reconstruction.Camera(camera_id);
            camera.SetCalibrated(true);
        }
        ba_options.solver_options.minimizer_progress_to_stdout=true;
        ba_options.print_summary=true;

        std::cout << "points added" << std::endl;
        // run bundle adjustment
        colmap::BundleAdjuster bundle_adjuster(ba_options, ba_config);
        std::cout<<"ba_config.NumImages: ," <<
        // CHECK_NOTNULL(reconstruction);
        // std::cout<<"CHECK NOT NULL SUCCEED"<<std::endl;
        bundle_adjuster.Solve(&reconstruction);
        // write the spline to the output directory
        std::string output_path = output_base + subdirectory.substr(subdirectory.find_last_of("/")+1) + "_spline.txt";
        std::vector<double> theta;
        std::vector<double> radii;
        colmap::Camera &camera = reconstruction.Camera(1);
        int num_control = (camera.Params().size() - 2) / 2 ;
        for ( int i =0; i< num_control; i++){
            theta.push_back(camera.Params()[i+2]);
            radii.push_back(camera.Params()[i+2+num_control]);
        }
        if (opt_extra_params == true){
        std::ofstream file(output_path);
         for ( int i =0; i< num_control; i++){
            file << camera.Params()[i+2] <<  " " << camera.Params()[i+2+num_control]<< std::endl;
        }}
        if (opt_extra_params == false){
            // std::string out_pose_path = std::string("/home/yihan/cvg/implicit_radial_sfm/experimental_scripts/babelcalib/reconstructions/refined_test_pose/")+ subdirectory.substr(subdirectory.find_last_of("/")+1) + "/pose.txt";

            colmap::CreateDirIfNotExists(output_base_pose);
            std::string out_pose_path = output_base_pose + subdirectory.substr(subdirectory.find_last_of("/")+1) + "_pose.txt";

            std::ofstream file(out_pose_path);
            for (const colmap::image_t image_id : reconstruction.RegImageIds()) {
                const colmap::Image& image = reconstruction.Image(image_id);
                file << image_id;
                for (size_t i = 0; i < 4; i++) {
                file << " "<< image.Qvec()[i];

                }
                for (size_t i = 0; i < 3; i++) {
                file << " "<< image.Tvec()[i];

                }
                file << std::endl;
            }
        }

    }
}

#endif