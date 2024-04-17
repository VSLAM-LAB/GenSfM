#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <ceres/ceres.h>

#include "implicit_local_bundle_adjustment.h"
#include "implicit_cost_functions.h"
#include "implicit_cost_matrix.h"
#include "implicit_intrinsic.h"
#include "implicit_utils.h"


namespace colmap {

ceres::LossFunction* setup_loss_function_ba_local(ImplicitBundleAdjustmentOptions::LossFunction loss_func, double loss_scale) {
    switch(loss_func) {
        case ImplicitBundleAdjustmentOptions::HUBER:
            return new ceres::HuberLoss(loss_scale);                
        case ImplicitBundleAdjustmentOptions::CAUCHY:
            return new ceres::CauchyLoss(loss_scale);                
        case ImplicitBundleAdjustmentOptions::TRIVIAL:
        default:
            return nullptr;                
    }
}

void local_bundle_adjustment(std::vector<std::vector<Eigen::Vector2d>> &points2D,
                                std::vector<Eigen::Vector3d> &points3D,
                                std::vector<std::vector<int>> &pointsInd,
                                std::vector<image_t>& image_ids,
                                CostMatrix &cost_matrix, const CostMatrixOptions& cm_opt,
                                std::vector<CameraPose> &poses, Eigen::Vector2d &pp, 
                                BundleAdjustmentConfig ba_config,
                                ImplicitBundleAdjustmentOptions ba_opt,
                                std::unordered_map<point3D_t, size_t> pointID_to_globalIndex) {

    if (ba_opt.upgrade_result) {
        std::cout << "pose upgrade for " << points2D.size() << " images" << std::endl;
        // First, upgrade the pose estimation
        int num_cams = points2D.size();
        // convert to suitable data structure for pose refinement
        std::vector<std::vector<Eigen::Vector3d>> points3D_sep(num_cams);
        for (int i = 0; i < num_cams; i++) {
            points3D_sep[i].resize(points2D[i].size());
            for (int j = 0; j < points2D[i].size(); j++) {
                points3D_sep[i][j] = points3D[pointsInd[i][j]];
            }
        }
        // print points3D dimensions before the pose refinement
        // std::cout << "before pose refinement points3D_sep size: " << points3D_sep.size() << std::endl;

        // poses = pose_refinement_multi(points2D, points3D_sep, cost_matrix, pp, poses, ba_opt);
        // std::cout << "after pose refinement points3D_sep size: " << points3D_sep.size() << std::endl;
    }
    // print points2D dimensions after the pose refinement

    // print points3D dimensions after the pose refinement
    std::cout << "points3D size: " << points3D.size() << std::endl;
  


    // if (ba_opt.filter_result) {
    //     std::cout << "ba_opt.filter_result" << std::endl;
    //     filter_result_ba_local(points2D, points3D, pointsInd, poses, pp, ba_opt);
    // }

    // exlucde points with too little occurrence in images
    std::vector<int> occur_count(points3D.size(), 0);
    for (size_t i = 0; i < pointsInd.size(); i++) {
        for (size_t j = 0; j < pointsInd[i].size(); j++) {
            occur_count[pointsInd[i][j]] += 1;
        }
    }

    int counter = 0;
    for (size_t i = 0; i < pointsInd.size(); i++) {
        auto ite_2D = points2D[i].begin();
        auto ite_ind = pointsInd[i].begin();
        while (ite_ind != pointsInd[i].end()) {
            // if (occur_count[*ite_ind] < ba_opt.min_curr_num) {
            if (occur_count[*ite_ind] < 1) {
                points2D[i].erase(ite_2D);
                pointsInd[i].erase(ite_ind);
                counter++;
            } else {
                ite_2D++;
                ite_ind++;
            }
        }
        for (int j = 0; j < pointsInd[i].size(); j++) {
            // if (occur_count[pointsInd[i][j]] < ba_opt.min_curr_num)
            if (occur_count[pointsInd[i][j]] < 1)
                std::cout << "exlusion error!" << std::endl;
        }
    }
    
    cost_matrix = build_cost_matrix_multi(points2D, cm_opt, pp);
    
    std::vector<Eigen::Quaterniond> qs;
    std::vector<Eigen::Vector3d> ts;

    int n_img = points2D.size();
    for (size_t k = 0; k < n_img; ++k) {
        qs.emplace_back(poses[k].q());
        ts.emplace_back(poses[k].t);
    }
    
    std::vector<std::vector<Eigen::Vector2d>> points2D_center = points2D;
    for (size_t cam_k = 0; cam_k < n_img; ++cam_k) {
        for (size_t i = 0; i < points2D[cam_k].size(); ++i) {
            points2D_center[cam_k][i] -= pp;
        }
    }

    // std::cout << "excluded point number: " << counter << std::endl;

    std::cout << "start bundle adjustment for " << points2D.size() << " images" << std::endl;
    // Then, do bundle adjustment
    int ite = 0;
    while (ite < ba_opt.max_ite_num) {
        
        std::cout << "ite number: " << ite << std::endl;
        double ratio = local_bundle_adjustment_inner(points2D_center, points3D, pointsInd, image_ids, cost_matrix, qs, ts, ba_config, ba_opt, pointID_to_globalIndex);
        ite++;

        std::cout << "optimize_projection done, decrease ratio: " << ratio << std::endl;
        if (ratio < ba_opt.stop_ratio) {
            std::cout << "ratio < " << ba_opt.stop_ratio << ", BA terminated" << std::endl;
            break;
        }
    }
    std::cout << "poses size: " << poses.size()<< std::endl;
    std::cout << "qs size: " << qs.size()<< std::endl;
    std::cout << "ts size: " << ts.size()<< std::endl;
    std::cout << "n_img: " << n_img << std::endl;
    for (size_t k = 0; k < n_img; ++k) {
        // std::cout << "k: " << k << std::endl;
        // std::cout << "poses[k].q: " << poses[k].q()<< std::endl;
        // std::cout << "qs[k]: " << qs[k]<< std::endl;
        // std::cout << "ts[k]: " << ts[k]<< std::endl;
        // std::cout << "poses[k].t: " << poses[k].t<< std::endl;
        poses[k].q(qs[k]);
        poses[k].t = ts[k];
    }
   std::cout << "BA terminated successfully "<< std::endl;
}


double local_bundle_adjustment_inner(const std::vector<std::vector<Eigen::Vector2d>> &points2D_center,
                                std::vector<Eigen::Vector3d> &points3D,
                                const std::vector<std::vector<int>> &pointsInd,
                                std::vector<image_t>& image_ids,
                                const CostMatrix &cost_matrix,
                                std::vector<Eigen::Quaterniond> &qs, std::vector<Eigen::Vector3d> &ts,
                                BundleAdjustmentConfig ba_config,
                                ImplicitBundleAdjustmentOptions ba_opt,
                                std::unordered_map<point3D_t, size_t> pointID_to_globalIndex) {

    size_t n_img = points2D_center.size();
    std::cout << "n_img: " << n_img << std::endl;
    std::vector<Eigen::Vector3d> points3D_new = points3D;


    ceres::Problem problem;

    // Radial reprojection error

    ceres::LossFunction* loss_function_radial = setup_loss_function_ba_local(ba_opt.loss_radial, ba_opt.loss_scale_radial); 
    // std::cout << "image_ids: " << image_ids.size() << std::endl;   

    // Adding all 3D points to the problem
    std::vector<double*> parameter_blocks;
    for (auto& point : points3D_new) {
        parameter_blocks.push_back(point.data());
    }

    // Adding all 3D points to the problem
    std::cout << "==========  Adding 3D points ==========" << std::endl;


    for (size_t cam_k = 0; cam_k < n_img; ++cam_k) {
        for (size_t i = 0; i < points2D_center[cam_k].size(); ++i) {
            int global_index = pointsInd[cam_k][i];
            ceres::CostFunction* reg_cost = BARadialReprojError::CreateCost(points2D_center[cam_k][i]);
    //         if(global_index == 1400){
    //             std::cout << "global_index is 0" << std::endl;
    //             // find the point3D_id as the index of pointsInd for which the global_index is 0
    
    //             point3D_t point3D_id = std::find_if(pointID_to_globalIndex.begin(), pointID_to_globalIndex.end(),
    //                                     [](const std::pair<point3D_t, size_t>& pair) {
    //                                         return pair.second == 1400;
    //                                     })->first; // Access the point3D_id from the iterator returned by find_if

    // std::cout << "point3D_id for global_index 0: " << point3D_id << std::endl;
    //         }   
            // std::cout << "Adding global_index: " << global_index << std::endl;
            
            problem.AddResidualBlock(reg_cost, loss_function_radial, qs[cam_k].coeffs().data(), ts[cam_k].data(), points3D_new[pointsInd[cam_k][i]].data());
        }
    }
    std::cout << "Total point in point3D_new: " << points3D_new.size() << std::endl;
    std::cout << "Total point in pointsInd: " << pointsInd.size() << std::endl;

    // // Set constant points based on ba_config
    // for (auto& [point3D_id, index] : pointID_to_globalIndex) {
    //     if (!ba_config.HasVariablePoint(point3D_id)) {
    //         // if (problem.HasParameterBlock(points3D[index].data())) {
    //         //     problem.SetParameterBlockConstant(points3D[index].data());
    //         // }
    //         std::cout << "index: " << index << std::endl;
    //         std::cout << "point3D_id not variable: " << point3D_id << std::endl;
    //         // check if has the parameter block
    //         if(problem.HasParameterBlock(points3D_new[index].data())){
    //             std::cout << "has the parameter block"  << std::endl;
    //             problem.SetParameterBlockConstant(points3D_new[index].data());            
    //         }else{
    //         std::cout << "doesn't have the parameter block"  << std::endl;
    //         }
    //     }
    // }


    // // ////////////// original code //////////////

    // size_t n_img = points2D_center.size();
    // std::cout << "n_img: " << n_img << std::endl;
    // std::vector<Eigen::Vector3d> points3D_new = points3D;

    // ceres::Problem problem;
    // //Radial reprojection error
    // ceres::LossFunction* loss_function_radial = setup_loss_function_ba_local(ba_opt.loss_radial, ba_opt.loss_scale_radial);
    // std::cout << "image_ids: " << image_ids.size() << std::endl;
    // for (size_t cam_k = 0; cam_k < n_img; ++cam_k) {
    //     for (size_t i = 0; i < points2D_center[cam_k].size(); ++i) {
    //         int global_index = pointID_to_globalIndex[pointsInd[cam_k][i]];
    //         ceres::CostFunction* reg_cost = BARadialReprojError::CreateCost(points2D_center[cam_k][i]);
    //         problem.AddResidualBlock(reg_cost, loss_function_radial, qs[cam_k].coeffs().data(), ts[cam_k].data(), points3D_new[pointsInd[cam_k][i]].data());
    //         // problem.AddResidualBlock(reg_cost, loss_function_radial, qs[cam_k].coeffs().data(), ts[cam_k].data(), points3D_new[global_index].data());
    //     }
    // }

    // /////////////////////////////////////////////////////////



    std::vector<std::vector<double*>> params(cost_matrix.pt_index.size());
    
    ceres::LossFunction* loss_function_dist = setup_loss_function_ba_local(ba_opt.loss_dist, ba_opt.loss_scale_dist);        

    // Implicit distortion cost (the cost matrix, regularization)
    for (size_t i = 0; i < cost_matrix.pt_index.size(); ++i) {
        ceres::CostFunction* reg_cost = BACostMatrixRowCost::CreateCost(
            points2D_center, pointsInd, points3D, points3D_new, cost_matrix.pt_index[i], cost_matrix.cam_index[i], cost_matrix.values[i], qs, ts, params[i]);

        problem.AddResidualBlock(reg_cost, loss_function_dist, params[i]);
    }

    // Setup parameterizations and constant parameter blocks
    for (size_t k = 0; k < n_img; ++k) {
        image_t img_id = image_ids[k];
        bool is_constant_pose = ba_config.HasConstantPose(img_id);
        double *q = qs[k].coeffs().data();
        if (!is_constant_pose) {
            // check if has parameter block
            // std::cout<<"has parameter block q?"<<problem.HasParameterBlock(q)<<std::endl;
            // std::cout<<"has parameter block ts?"<<problem.HasParameterBlock(ts[k].data())<<std::endl;
        
            problem.SetParameterization(q, new ceres::EigenQuaternionParameterization());
        } else {
            // std::cout << "has parameter block q? "<<problem.HasParameterBlock(q) << std::endl;
            // std::cout << "has parameter block ts? "<<problem.HasParameterBlock(ts[k].data()) << std::endl;
            problem.SetParameterBlockConstant(q);
            problem.SetParameterBlockConstant(ts[k].data());
        }

        problem.SetParameterization(q, new ceres::EigenQuaternionParameterization());

        if (pointsInd[k].size() != points2D_center[k].size()) {
            std::cout << "size error" << std::endl;
        }
    }

    ceres::Solver::Options options;
    options.max_num_iterations = 5;
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    options.minimizer_progress_to_stdout = ba_opt.verbose; // true if you want more debug output
    options.num_threads = 8;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    points3D = points3D_new;

    // return the decrease ratio
    return (summary.initial_cost - summary.final_cost) / double(summary.initial_cost);
}

void filter_result_ba_local(std::vector<std::vector<Eigen::Vector2d>> &points2D,
                                std::vector<Eigen::Vector3d> &points3D,
                                std::vector<std::vector<int>> &pointsInd,
                                std::vector<CameraPose>& poses, Eigen::Vector2d &pp, 
                                ImplicitBundleAdjustmentOptions ba_opt) {

    std::cout << "points2D size" << points2D.size()<< std::endl;
    
   
    std::vector<std::vector<Eigen::Vector3d>> points3D_sep(points2D.size());
    for (int i = 0; i < points2D.size(); i++) {
        points3D_sep[i].resize(points2D[i].size());
        for (int j = 0; j < points2D[i].size(); j++) {
            points3D_sep[i][j] = points3D[pointsInd[i][j]];
        }
    }
    std::cout << "points3D_sep size" << points3D_sep.size()<< std::endl;
    
    std::vector<std::vector<double>> fs_diff;
    calculate_fmed_diff(points2D, points3D_sep, poses, pp, fs_diff);
    // calculate the total number of points in points3D_sep
    int total_num = 0;
    for (int i = 0; i < points3D_sep.size(); i++) {
        total_num += points3D_sep[i].size();
    }
    std::cout << "points3D_sep number before calculating fmed_diff: " << total_num << std::endl;


    int counter = 0;    
    for (int i = 0; i < points2D.size(); i++) {
        int ori_size = points2D[i].size();
        exclude_outliers(fs_diff[i], points2D[i], pointsInd[i], true, ba_opt.filter_thres);
        counter += ori_size - points2D[i].size();
    }
    std::cout << "filtering number: " << counter << std::endl;
}


}