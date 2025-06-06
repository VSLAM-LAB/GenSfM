// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Viktor Larsson

#ifndef COLMAP_SRC_ESTIMATORS_RADIAL_ABSOLUTE_POSE_H_
#define COLMAP_SRC_ESTIMATORS_RADIAL_ABSOLUTE_POSE_H_

#include <array>
#include <vector>
#include "estimators/implicit_cost_matrix.h"
#include "estimators/implicit_pose_refinement.h"
#include "estimators/implicit_intrinsic.h"
#include "estimators/implicit_utils.h"
#include "estimators/implicit_camera_pose.h"

#include <Eigen/Core>

#include "util/alignment.h"
#include "util/types.h"

namespace colmap {

// Analytic solver for the 1D radial absolute pose problem.
//
// The algorithm is based on the following paper:
//
//   Kukelova et al., Real-Time Solution to the Absolute Pose Problem with
//                    Unknown Radial Distortion and Focal Length, ICCV 2013
class RadialP5PEstimator {
 public:
  // The 2D image feature observations.
  typedef Eigen::Vector2d X_t;
  // The observed 3D features in the world frame.
  typedef Eigen::Vector3d Y_t;
  // The transformation from the world to the camera frame.
  typedef Eigen::Matrix3x4d M_t;

  // The minimum number of samples needed to estimate a model.
  static const int kMinNumSamples = 5;

  // Estimate radial camera pose from 5 correspondences
  //
  // @param points2D   2D images points
  // @param points3D   3D world points
  //
  // @return           Camera pose as a 3x4 matrix.
  static std::vector<M_t> Estimate(const std::vector<X_t>& points2D,
                                   const std::vector<Y_t>& points3D, bool initial = false);

  // Calculate the squared reprojection error given a set of 2D-3D point
  // correspondences and a projection matrix.
  //
  // @param points2D     2D image points as Nx2 matrix.
  // @param points3D     3D world points as Nx3 matrix.
  // @param proj_matrix  3x4 projection matrix.
  // @param residuals    Output vector of residuals.
  static void Residuals(const std::vector<X_t>& points2D,
                        const std::vector<Y_t>& points3D,
                        const M_t& proj_matrix, std::vector<double>* residuals);
};

// Estimates the foward offset t3 by assuming a single focal length
// This is only used for normalizing scale of the reconstruction and
// for visualization purposes.
double EstimateRadialCameraForwardOffset(
    const Eigen::Matrix3x4d proj_matrix,
    const std::vector<Eigen::Vector2d> points2D,
    const std::vector<Eigen::Vector3d> points3D,
    bool *negative_focal);

// Estimates the foward translation t3 using implicit distortion model
CameraPose EstimateCameraForwardOffsetImplictDistortion(
    const Eigen::Matrix3x4d proj_matrix,
    const std::vector<Eigen::Vector2d> &points2D,
    const std::vector<Eigen::Vector3d> &points3D,
    const Eigen::Vector2d pp = Eigen::Vector2d(0.0,0.0));
    
    }  // namespace colmap

#endif  // COLMAP_SRC_ESTIMATORS_RADIAL_ABSOLUTE_POSE_H_
