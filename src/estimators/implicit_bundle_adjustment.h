#ifndef IMPLICIT_DIST_BUNDLE_AJDUSTMENT_H_
#define IMPLICIT_DIST_BUNDLE_AJDUSTMENT_H_

#include <vector>
#include <Eigen/Dense>

#include "implicit_camera_pose.h"
#include "implicit_pose_refinement.h"

namespace colmap {

    struct ImplicitBundleAdjustmentOptions : PoseRefinementOptions {
        int max_ite_num = 10;
        // int max_ite_num = 5;
        int min_curr_num = 2;
        double stop_ratio = 1e-3;
        double filter_thres = 10;

        bool upgrade_result = false;
        // bool filter_result = true; // filter before BA starts
        bool filter_result = false; // disable filtering for now

        ImplicitBundleAdjustmentOptions clone() const {
            ImplicitBundleAdjustmentOptions copy = *this;
            return copy;
        }
    };



    void bundle_adjustment(std::vector<std::vector<Eigen::Vector2d>>& points2D,
        std::vector<Eigen::Vector3d>& points3D,
        std::vector<std::vector<int>>& pointsInd,
        CostMatrix& cost_matrix, const CostMatrixOptions& cm_opt,
        std::vector<CameraPose>& poses, Eigen::Vector2d& pp,
        ImplicitBundleAdjustmentOptions ba_opt);

    double bundle_adjustment_inner(const std::vector<std::vector<Eigen::Vector2d>>& points2D_center,
        std::vector<Eigen::Vector3d>& points3D,
        const std::vector<std::vector<int>>& pointsInd,
        const CostMatrix& cost_matrix,
        std::vector<Eigen::Quaterniond>& qs, std::vector<Eigen::Vector3d>& ts,
        ImplicitBundleAdjustmentOptions ba_opt);


    void filter_result_ba(std::vector<std::vector<Eigen::Vector2d>>& points2D,
        std::vector<Eigen::Vector3d>& points3D,
        std::vector<std::vector<int>>& pointsInd,
        std::vector<CameraPose>& poses, Eigen::Vector2d& pp,
        ImplicitBundleAdjustmentOptions ba_opt);
}  // namespace colmap

#endif