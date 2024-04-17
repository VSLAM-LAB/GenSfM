#ifndef IMPLICIT_DIST_LOCAL_BUNDLE_AJDUSTMENT_H_
#define IMPLICIT_DIST_LOCAL_BUNDLE_AJDUSTMENT_H_

#include <vector>
#include <Eigen/Dense>

#include "optim/bundle_adjustment.h"
#include "implicit_camera_pose.h"
#include "implicit_pose_refinement.h"
#include "implicit_bundle_adjustment.h"



namespace colmap {

    // struct ImplicitBundleAdjustmentOptions : PoseRefinementOptions {
    //     int max_ite_num = 10;
    //     int min_curr_num = 5;
    //     double stop_ratio = 1e-3;
    //     double filter_thres = 10;
        
    //     bool upgrade_result = false;
    //     bool filter_result = true; // filter before BA starts

    //     ImplicitBundleAdjustmentOptions clone() const {
    //         ImplicitBundleAdjustmentOptions copy = *this;
    //         return copy;
    //     }
    // };
    
    

    void local_bundle_adjustment(std::vector<std::vector<Eigen::Vector2d>> &points2D,
                            std::vector<Eigen::Vector3d> &points3D,
                            std::vector<std::vector<int>> &pointsInd,
                            std::vector<image_t>& image_ids,
                            CostMatrix &cost_matrix, const CostMatrixOptions& cm_opt,
                            std::vector<CameraPose> &poses, Eigen::Vector2d &pp, 
                            BundleAdjustmentConfig ba_config,
                            ImplicitBundleAdjustmentOptions ba_opt,
                            std::unordered_map<point3D_t, size_t> pointID_to_globalIndex);

    double local_bundle_adjustment_inner(const std::vector<std::vector<Eigen::Vector2d>> &points2D_center,
                            std::vector<Eigen::Vector3d> &points3D,
                            const std::vector<std::vector<int>> &pointsInd,
                            std::vector<image_t>& image_ids,
                            const CostMatrix &cost_matrix,
                            std::vector<Eigen::Quaterniond> &qs, std::vector<Eigen::Vector3d> &ts, 
                            BundleAdjustmentConfig ba_config,
                            ImplicitBundleAdjustmentOptions ba_opt,
                            std::unordered_map<point3D_t, size_t> pointID_to_globalIndex);

    void filter_result_ba_local(std::vector<std::vector<Eigen::Vector2d>> &points2D,
                            std::vector<Eigen::Vector3d> &points3D,
                            std::vector<std::vector<int>> &pointsInd,
                            std::vector<CameraPose> &poses, Eigen::Vector2d &pp, 
                            ImplicitBundleAdjustmentOptions ba_opt);



}  // namespace colmap

#endif