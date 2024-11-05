#include "implicit_pose_refinement.h"
#include "implicit_cost_functions.h"
#include "implicit_utils.h"
#include "manifold.h"
#include <ceres/ceres.h>

namespace colmap{
ceres::LossFunction* setup_loss_function(PoseRefinementOptions::LossFunction loss_func, double loss_scale) {
    switch(loss_func) {
        case PoseRefinementOptions::HUBER:
            return new ceres::HuberLoss(loss_scale);                
        case PoseRefinementOptions::CAUCHY:
            return new ceres::CauchyLoss(loss_scale);                
        case PoseRefinementOptions::TRIVIAL:
        default:
            return nullptr;                
    }
}

CameraPose pose_refinement(const std::vector<Eigen::Vector2d> &points2D, 
                                const std::vector<Eigen::Vector3d> &points3D,
                                const CostMatrix &cost_matrix, const Eigen::Vector2d &pp,
                                const CameraPose &initial_pose, PoseRefinementOptions refinement_opt){

    std::vector<CameraPose> output;
    const std::vector<CameraPose> initial_poses = {initial_pose};
    
    output = pose_refinement_multi({points2D}, {points3D}, cost_matrix, pp, initial_poses, refinement_opt);
    return output[0];
}

std::vector<CameraPose> pose_refinement_multi(
                                const std::vector<std::vector<Eigen::Vector2d>> &points2D, 
                                const std::vector<std::vector<Eigen::Vector3d>> &points3D,
                                const CostMatrix &cost_matrix, const Eigen::Vector2d &pp,
                                const std::vector<CameraPose> &initial_poses, PoseRefinementOptions refinement_opt) {

    size_t n_img = points2D.size();    

    std::vector<std::vector<Eigen::Vector2d>> points2D_center = points2D;
     for (size_t cam_k = 0; cam_k < n_img; ++cam_k) {
        for (size_t i = 0; i < points2D[cam_k].size(); ++i) {
            points2D_center[cam_k][i] -= pp;
        }
     }

    std::vector<Eigen::Quaterniond> qs;
    std::vector<Eigen::Vector3d> ts;

    for (size_t k = 0; k < n_img; ++k) {
        qs.emplace_back(initial_poses[k].q());
        ts.emplace_back(initial_poses[k].t);
    }

    ceres::Problem problem;
    
    // Radial projection error
    ceres::LossFunction* loss_function_radial = setup_loss_function(refinement_opt.loss_radial, refinement_opt.loss_scale_radial);        

    for (size_t cam_k = 0; cam_k < n_img; ++cam_k) {
        for (size_t i = 0; i < points2D[cam_k].size(); ++i) {
            ceres::CostFunction* reg_cost = RadialReprojError::CreateCost(points2D_center[cam_k][i], points3D[cam_k][i]);
            problem.AddResidualBlock(reg_cost, loss_function_radial, qs[cam_k].coeffs().data(), ts[cam_k].data());
        }
        }

    // This is used for setting up the dynamic parameter vectors
    // (workaround for having duplicate parameters in the blocks)
    std::vector<std::vector<double*>> params(cost_matrix.pt_index.size());
    
    // Implicit distortion cost (the cost matrix, regularization)
    ceres::LossFunction* loss_function_dist = setup_loss_function(refinement_opt.loss_dist, refinement_opt.loss_scale_dist);        


    for (size_t i = 0; i < cost_matrix.pt_index.size(); ++i) {
        ceres::CostFunction* reg_cost = CostMatrixRowCost::CreateCost(
                points2D_center, points3D, cost_matrix.pt_index[i], cost_matrix.cam_index[i],
                cost_matrix.values[i], qs, ts, params[i]);

        problem.AddResidualBlock(reg_cost, loss_function_dist, params[i]);
    }

    // Setup parameterizations and constant parameter blocks
    for (size_t k = 0; k < n_img; ++k) {
        if (!problem.HasParameterBlock(qs[k].coeffs().data()))
            continue;
        double *q = qs[k].coeffs().data();
        double *t = ts[k].data();

        SetEigenQuaternionManifold(&problem, q);
    }

    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    // options.minimizer_progress_to_stdout = refinement_opt.verbose; // true if you want more debug output
    options.minimizer_progress_to_stdout = false; // true if you want more debug output
    options.num_threads = 8;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::vector<CameraPose> output;
    for (size_t k = 0; k < n_img; ++k) {
        CameraPose pose;
        pose.q(qs[k]);
        pose.t = ts[k];
        output.push_back(pose);
    }
    return output;
}

CameraPose pose_refinement_multi_point2d_single_image(
                                const std::vector<std::vector<Eigen::Vector2d>> &points2D, 
                                const std::vector<std::vector<Eigen::Vector3d>> &points3D,
                                const CostMatrix &cost_matrix, const Eigen::Vector2d &pp,
                                const std::vector<CameraPose> &initial_poses,
                                int image_id, PoseRefinementOptions refinement_opt) {

    size_t n_img = points2D.size();
    if (image_id < 0 || image_id >= n_img) {
        std::cerr << "Invalid image_id" << std::endl;
    }

    std::vector<std::vector<Eigen::Vector2d>> points2D_center = points2D;
     for (size_t cam_k = 0; cam_k < n_img; ++cam_k) {
        for (size_t i = 0; i < points2D[cam_k].size(); ++i) {
            points2D_center[cam_k][i] -= pp;
        }
     }

    std::vector<Eigen::Quaterniond> qs;
    std::vector<Eigen::Vector3d> ts;

    for (size_t k = 0; k < n_img; ++k) {
        qs.emplace_back(initial_poses[k].q());
        ts.emplace_back(initial_poses[k].t);
    }

    ceres::Problem problem;
    
    // Radial projection error
    ceres::LossFunction* loss_function_radial = setup_loss_function(refinement_opt.loss_radial, refinement_opt.loss_scale_radial);        

    for (size_t cam_k = 0; cam_k < n_img; ++cam_k) {
        for (size_t i = 0; i < points2D[cam_k].size(); ++i) {
            ceres::CostFunction* reg_cost = RadialReprojError::CreateCost(points2D_center[cam_k][i], points3D[cam_k][i]);
            problem.AddResidualBlock(reg_cost, loss_function_radial, qs[cam_k].coeffs().data(), ts[cam_k].data());
        }
        }

    // This is used for setting up the dynamic parameter vectors
    // (workaround for having duplicate parameters in the blocks)
    std::vector<std::vector<double*>> params(cost_matrix.pt_index.size());
    
    // Implicit distortion cost (the cost matrix, regularization)
    ceres::LossFunction* loss_function_dist = setup_loss_function(refinement_opt.loss_dist, refinement_opt.loss_scale_dist);        


    for (size_t i = 0; i < cost_matrix.pt_index.size(); ++i) {
        ceres::CostFunction* reg_cost = CostMatrixRowCost::CreateCost(
                points2D_center, points3D, cost_matrix.pt_index[i], cost_matrix.cam_index[i],
                cost_matrix.values[i], qs, ts, params[i]);

        problem.AddResidualBlock(reg_cost, loss_function_dist, params[i]);
    }

    // Setup parameterizations and constant parameter blocks
    for (size_t k = 0; k < n_img; ++k) {
        double *q = qs[k].coeffs().data();
        double *t = ts[k].data();
        
        if (!problem.HasParameterBlock(q))
            continue;
        SetEigenQuaternionManifold(&problem, q);

        if (k != image_id) {
            problem.SetParameterBlockConstant(q);
            problem.SetParameterBlockConstant(t);
        }
    }

    ceres::Solver::Options options;
    options.max_num_iterations = 100;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    // options.minimizer_progress_to_stdout = refinement_opt.verbose; // true if you want more debug output
    options.minimizer_progress_to_stdout = false; // true if you want more debug output
    options.num_threads = 8;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    CameraPose pose;
    pose.q(qs[image_id]);
    pose.t = ts[image_id];
    // output.push_back(pose);
    return pose;
}

void filter_result_pose_refinement(std::vector<Eigen::Vector2d> &points2D,
                                std::vector<Eigen::Vector3d> &points3D,
                                const CameraPose& pose, const Eigen::Vector2d &pp, 
                                PoseRefinementOptions refinement_opt) {
    
    std::vector<std::vector<Eigen::Vector2d>> points2D_vec = {points2D};
    std::vector<std::vector<Eigen::Vector3d>> points3D_vec = {points3D};
    filter_result_pose_refinement_multi(points2D_vec, points3D_vec, {pose}, pp, refinement_opt);

    points2D = points2D_vec[0];
    points3D = points3D_vec[0];
}

void filter_result_pose_refinement_multi(std::vector<std::vector<Eigen::Vector2d>> &points2D,
                                std::vector<std::vector<Eigen::Vector3d>> &points3D,
                                const std::vector<CameraPose>& poses, const Eigen::Vector2d &pp, 
                                PoseRefinementOptions refinement_opt) {

    std::vector<std::vector<double>> fs_diff;
    calculate_fmed_diff(points2D, points3D, poses, pp, fs_diff);
    
    int counter = 0;
    for (int i = 0; i < points2D.size(); i++) {
        int ori_size = points2D[i].size();
        exclude_outliers(fs_diff[i], points2D[i], points3D[i], true, refinement_opt.filter_thres);
        counter += ori_size - points2D[i].size();
    }
    // std::cout << "filtered number of entries: " << counter << std::endl;
}

void filter_result_pose_refinement_multi(std::vector<std::vector<Eigen::Vector2d>> &points2D,
                                std::vector<std::vector<Eigen::Vector3d>> &points3D,
                                const std::vector<CameraPose>& poses, const Eigen::Vector2d &pp, 
                                std::vector<std::vector<bool>>& is_outlier_multi,
                                PoseRefinementOptions refinement_opt) {

    std::vector<std::vector<double>> fs_diff;
    calculate_fmed_diff(points2D, points3D, poses, pp, fs_diff);

    is_outlier_multi.clear();
    is_outlier_multi.resize(points2D.size());
    
    int counter = 0;
    for (int i = 0; i < points2D.size(); i++) {
        int ori_size = points2D[i].size();
        exclude_outliers(fs_diff[i], points2D[i], points3D[i], true, refinement_opt.filter_thres, &(is_outlier_multi[i]));
        counter += ori_size - points2D[i].size();
    }
    // std::cout << "filtered number of entries: " << counter << std::endl;
}

void pose_refinement_1D_radial(const std::vector<Eigen::Vector2d> &points2D,
                            const std::vector<Eigen::Vector3d> &points3D,
                            CameraPose *pose, Eigen::Vector2d *pp, 
                            PoseRefinement1DRadialOptions opt) {

    std::vector<CameraPose> poses;
    poses.push_back(*pose);
    joint_pose_refinement_1D_radial({points2D}, {points3D}, &poses, pp, opt);
    *pose = poses[0];
}



void joint_pose_refinement_1D_radial(const std::vector<std::vector<Eigen::Vector2d>> &points2D,
                            const std::vector<std::vector<Eigen::Vector3d>> &points3D,
                            std::vector<CameraPose> *poses, Eigen::Vector2d *pp, 
                            PoseRefinement1DRadialOptions opt) {

    const size_t num_cams = points2D.size();
    std::vector<Eigen::Quaterniond> qs;
    std::vector<Eigen::Vector2d> ts;

    for (size_t k = 0; k < num_cams; ++k) {
        qs.emplace_back(poses->at(k).q());
        ts.push_back(poses->at(k).t.topRows<2>());
    }

    ceres::Problem problem;
    ceres::LossFunction* loss_function = new ceres::HuberLoss(6.0);

    for (size_t cam_k = 0; cam_k < num_cams; ++cam_k) {
        for (size_t i = 0; i < points2D[cam_k].size(); ++i) {
            ceres::CostFunction* reg_cost = RadialReprojError_Implicit::CreateCost(points2D[cam_k][i], points3D[cam_k][i]);
            problem.AddResidualBlock(reg_cost, loss_function, qs[cam_k].coeffs().data(), ts[cam_k].data(), pp->data());
        }
    }
    
    for (size_t cam_k = 0; cam_k < num_cams; ++cam_k) {
        SetEigenQuaternionManifold(&problem, qs[cam_k].coeffs().data());
    }


    ceres::Solver::Options options;
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    // options.minimizer_progress_to_stdout = opt.verbose;
    options.minimizer_progress_to_stdout = false;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    for (size_t cam_k = 0; cam_k < num_cams; ++cam_k) {
        poses->at(cam_k).q(qs[cam_k]);
        poses->at(cam_k).t.topRows<2>() = ts[cam_k];
    }
}
} // namespace colmap