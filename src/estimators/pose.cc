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
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#include "estimators/manifold.h"
#include "estimators/pose.h"

#include "base/camera_models.h"
#include "base/cost_functions.h"
#include "base/essential_matrix.h"
#include "base/pose.h"
#include "estimators/absolute_pose.h"
#include "estimators/essential_matrix.h"
#include "estimators/radial_absolute_pose.h"
#include "optim/bundle_adjustment.h"
#include "util/matrix.h"
#include "util/misc.h"
#include "util/threading.h"

namespace colmap {
namespace {

typedef LORANSAC<P3PEstimator, EPNPEstimator> AbsolutePoseRANSAC;

void EstimateAbsolutePoseKernel(const Camera& camera,
                                const double focal_length_factor,
                                const std::vector<Eigen::Vector2d>& points2D,
                                const std::vector<Eigen::Vector3d>& points3D,
                                const RANSACOptions& options,
                                AbsolutePoseRANSAC::Report* report) {
  // Scale the focal length by the given factor.
  Camera scaled_camera = camera;
  const std::vector<size_t>& focal_length_idxs = camera.FocalLengthIdxs();
  for (const size_t idx : focal_length_idxs) {
    scaled_camera.Params(idx) *= focal_length_factor;
  }

  // Normalize image coordinates with current camera hypothesis.
  std::vector<Eigen::Vector2d> points2D_N(points2D.size());
  for (size_t i = 0; i < points2D.size(); ++i) {
    points2D_N[i] = scaled_camera.ImageToWorld(points2D[i]);
  }

  // Estimate pose for given focal length.
  auto custom_options = options;
  custom_options.max_error =
      scaled_camera.ImageToWorldThreshold(options.max_error);
  AbsolutePoseRANSAC ransac(custom_options);
  *report = ransac.Estimate(points2D_N, points3D);
}

}  // namespace

bool EstimateRadialAbsolutePose(const AbsolutePoseEstimationOptions& options,
                                const std::vector<Eigen::Vector2d>& points2D,
                                const std::vector<Eigen::Vector3d>& points3D,
                                Eigen::Vector4d* qvec, Eigen::Vector3d* tvec,
                                Camera* camera, size_t* num_inliers,
                                std::vector<char>* inlier_mask) {
  // CHECK_EQ(camera->ModelId(), Radial1DCameraModel::model_id);
  CHECK(camera->ModelId() == Radial1DCameraModel::model_id ||
      camera->ModelId() == ImplicitDistortionModel::model_id);
  
  // Subtract principal point
  std::vector<Eigen::Vector2d> points2D_N(points2D.size());
  for (size_t i = 0; i < points2D.size(); ++i) {
    points2D_N[i] = camera->ImageToWorld(points2D[i]); //normalized 2D points
  }

  LORANSAC<RadialP5PEstimator, RadialP5PEstimator, MEstimatorSupportMeasurer>
      ransac(options.ransac_options);
  auto report = ransac.Estimate(points2D_N, points3D); //// The transformation from the world to the camera frame, a 3x4 matrix is stored in report.model

  if (!report.success) {
    return false;
  }

  *num_inliers = report.support.num_inliers;
  *inlier_mask = report.inlier_mask;
  // output the inlier ratio to a file
  // std::ofstream file("~/cvg/implicit_radial_sfm/inliers.txt", std::ios_base::app);
  

  if (*num_inliers == 0) {
    return false;
  }

  // We guesstimate the forward translation as well
  // Note that this is not used during the reconstruction, but is only
  // for scale and visualizuation purposes
  std::vector<Eigen::Vector2d> inlier_points2D;
  std::vector<Eigen::Vector2d> inlier_points2D_normalized;
  std::vector<Eigen::Vector3d> inlier_points3D;
  inlier_points2D.reserve(*num_inliers);
  inlier_points3D.reserve(*num_inliers);
  for (size_t i = 0; i < inlier_mask->size(); ++i) {
    if ((*inlier_mask)[i]) {
      inlier_points2D_normalized.push_back(points2D_N[i]); // normalizing 2D points within the implicit distortion model
      inlier_points2D.push_back(points2D[i]);
      inlier_points3D.push_back(points3D[i]);
    }
  }
  ;
  double tz = EstimateRadialCameraForwardOffset(report.model, inlier_points2D_normalized,
                                                inlier_points3D, nullptr);
  // print out tz
  // std::cout << "tz:   " << tz << std::endl;
  //append tz into file "tz.txt"
  // std::ofstream file("/home/ivonne/radialsfm/tz.txt", std::ios_base::app);
  // file << tz << std::endl;
  // file.close();
  // print out the image id
  // // Extract pose parameters.
  *qvec = RotationMatrixToQuaternion(report.model.leftCols<3>());
  *tvec = report.model.rightCols<1>();
  (*tvec)(2) += tz;
  
  //principal point
  Eigen::Vector2d pp = Eigen::Vector2d(camera->PrincipalPointX(), camera->PrincipalPointY());

  // Estimate the forward translation using implicit distortion model.
  CameraPose implicit_pose = EstimateCameraForwardOffsetImplictDistortion(report.model, inlier_points2D, 
                                                inlier_points3D, pp);

  // //print out tz_imp
  // std::cout << "tz_imp: " << tz_imp << std::endl;
  // //append tz_imp into file "tz_imp.txt"
  // std::ofstream file2("/home/ivonne/radialsfm/tz_imp.txt", std::ios_base::app);
  // file2 << tz_imp << std::endl;
  // file2.close();

  

  // Extract pose parameters from implicit_pose
  if(inlier_points2D.size()>0){
  *qvec = implicit_pose.q_vec;
  *tvec = implicit_pose.t;
  }


  if (IsNaN(*qvec) || IsNaN(*tvec)) {
    return false;
  }

  return true;
}

bool EstimateAbsolutePose(const AbsolutePoseEstimationOptions& options,
                          const std::vector<Eigen::Vector2d>& points2D,
                          const std::vector<Eigen::Vector3d>& points3D,
                          Eigen::Vector4d* qvec, Eigen::Vector3d* tvec,
                          Camera* camera, size_t* num_inliers,
                          std::vector<char>* inlier_mask) {
  options.Check();

  if (camera->ModelId() == Radial1DCameraModel::model_id ||
      camera->ModelId() == ImplicitDistortionModel::model_id ){
    return EstimateRadialAbsolutePose(options, points2D, points3D, qvec, tvec,
                                      camera, num_inliers, inlier_mask);
  }

  std::vector<double> focal_length_factors;
  if (options.estimate_focal_length) {
    // Generate focal length factors using a quadratic function,
    // such that more samples are drawn for small focal lengths
    focal_length_factors.reserve(options.num_focal_length_samples + 1);
    const double fstep = 1.0 / options.num_focal_length_samples;
    const double fscale =
        options.max_focal_length_ratio - options.min_focal_length_ratio;
    for (double f = 0; f <= 1.0; f += fstep) {
      focal_length_factors.push_back(options.min_focal_length_ratio +
                                     fscale * f * f);
    }
  } else {
    focal_length_factors.reserve(1);
    focal_length_factors.push_back(1);
  }

  std::vector<std::future<void>> futures;
  futures.resize(focal_length_factors.size());
  std::vector<typename AbsolutePoseRANSAC::Report,
              Eigen::aligned_allocator<typename AbsolutePoseRANSAC::Report>>
      reports;
  reports.resize(focal_length_factors.size());

  ThreadPool thread_pool(std::min(
      options.num_threads, static_cast<int>(focal_length_factors.size())));

  for (size_t i = 0; i < focal_length_factors.size(); ++i) {
    futures[i] = thread_pool.AddTask(
        EstimateAbsolutePoseKernel, *camera, focal_length_factors[i], points2D,
        points3D, options.ransac_options, &reports[i]);
  }

  double focal_length_factor = 0;
  Eigen::Matrix3x4d proj_matrix;
  *num_inliers = 0;
  inlier_mask->clear();

  // Find best model among all focal lengths.
  for (size_t i = 0; i < focal_length_factors.size(); ++i) {
    futures[i].get();
    const auto report = reports[i];
    if (report.success && report.support.num_inliers > *num_inliers) {
      *num_inliers = report.support.num_inliers;
      proj_matrix = report.model;
      *inlier_mask = report.inlier_mask;
      focal_length_factor = focal_length_factors[i];
    }
  }

  if (*num_inliers == 0) {
    return false;
  }

  // Scale output camera with best estimated focal length.
  if (options.estimate_focal_length && *num_inliers > 0) {
    const std::vector<size_t>& focal_length_idxs = camera->FocalLengthIdxs();
    for (const size_t idx : focal_length_idxs) {
      camera->Params(idx) *= focal_length_factor;
    }
  }

  // Extract pose parameters.
  *qvec = RotationMatrixToQuaternion(proj_matrix.leftCols<3>());
  *tvec = proj_matrix.rightCols<1>();

  if (IsNaN(*qvec) || IsNaN(*tvec)) {
    return false;
  }

  return true;
}

size_t EstimateRelativePose(const RANSACOptions& ransac_options,
                            const std::vector<Eigen::Vector2d>& points1,
                            const std::vector<Eigen::Vector2d>& points2,
                            Eigen::Vector4d* qvec, Eigen::Vector3d* tvec) {
  RANSAC<EssentialMatrixFivePointEstimator> ransac(ransac_options);
  const auto report = ransac.Estimate(points1, points2);

  if (!report.success) {
    return 0;
  }

  std::vector<Eigen::Vector2d> inliers1(report.support.num_inliers);
  std::vector<Eigen::Vector2d> inliers2(report.support.num_inliers);

  size_t j = 0;
  for (size_t i = 0; i < points1.size(); ++i) {
    if (report.inlier_mask[i]) {
      inliers1[j] = points1[i];
      inliers2[j] = points2[i];
      j += 1;
    }
  }

  Eigen::Matrix3d R;

  std::vector<Eigen::Vector3d> points3D;
  PoseFromEssentialMatrix(report.model, inliers1, inliers2, &R, tvec,
                          &points3D);

  *qvec = RotationMatrixToQuaternion(R);

  if (IsNaN(*qvec) || IsNaN(*tvec)) {
    return 0;
  }

  return points3D.size();
}

bool RefineAbsolutePose(const AbsolutePoseRefinementOptions& options,
                        const std::vector<char>& inlier_mask,
                        const std::vector<Eigen::Vector2d>& points2D,
                        const std::vector<Eigen::Vector3d>& points3D,
                        Eigen::Vector4d* qvec, Eigen::Vector3d* tvec,
                        Camera* camera, bool optimize_tz) {
  CHECK_EQ(inlier_mask.size(), points2D.size());
  CHECK_EQ(points2D.size(), points3D.size());
  options.Check();

  ceres::LossFunction* loss_function =
      new ceres::CauchyLoss(options.loss_function_scale);

  double* camera_params_data = camera->ParamsData();
  double* qvec_data = qvec->data();
  double* tvec_data = tvec->data();

  std::vector<Eigen::Vector3d> points3D_copy = points3D;

  ceres::Problem problem;

  for (size_t i = 0; i < points2D.size(); ++i) {
    // Skip outlier observations
    if (!inlier_mask[i]) {
      continue;
    }

    ceres::CostFunction* cost_function = nullptr;
    if(!optimize_tz){
    cost_function = BundleAdjustmentCostFunction<Radial1DCameraModel>::Create(
        points2D[i]);
        }
    else{
    switch (camera->ModelId()) {
#define CAMERA_MODEL_CASE(CameraModel)                                  \
  case CameraModel::kModelId:                                           \
    cost_function =                                                     \
        BundleAdjustmentCostFunction<CameraModel>::Create(points2D[i]); \
    break;

      CAMERA_MODEL_SWITCH_CASES

#undef CAMERA_MODEL_CASE
    }
    }

    problem.AddResidualBlock(cost_function, loss_function, qvec_data, tvec_data,
                             points3D_copy[i].data(), camera_params_data);
    problem.SetParameterBlockConstant(points3D_copy[i].data());
  }
  
  if (problem.NumResiduals() > 0) {
    // Quaternion parameterization.
    *qvec = NormalizeQuaternion(*qvec);
    SetQuaternionManifold(&problem, qvec_data);

  
    if (camera->ModelId() == Radial1DCameraModel::model_id ||
         (camera->ModelId() == ImplicitDistortionModel::model_id && !optimize_tz)) {
    // if (camera->ModelId() == Radial1DCameraModel::model_id ||
        //  camera->ModelId() == ImplicitDistortionModel::model_id) {
    //  if (camera->ModelId() == Radial1DCameraModel::model_id) {
      // Only optimize over the first two elements of the translation vector for
      // radial cameras
      SetSubsetManifold(3, {2}, &problem, tvec_data);
    }

    // Camera parameterization.
    if (!options.refine_focal_length && !options.refine_extra_params) {
      problem.SetParameterBlockConstant(camera->ParamsData());
    } else {
      // Always set the principal point as fixed.
      std::vector<int> camera_params_const;
      const std::vector<size_t>& principal_point_idxs =
          camera->PrincipalPointIdxs();
      camera_params_const.insert(camera_params_const.end(),
                                 principal_point_idxs.begin(),
                                 principal_point_idxs.end());

      if (!options.refine_focal_length) {
        const std::vector<size_t>& focal_length_idxs =
            camera->FocalLengthIdxs();
        camera_params_const.insert(camera_params_const.end(),
                                   focal_length_idxs.begin(),
                                   focal_length_idxs.end());
      }

      if (!options.refine_extra_params) {
        const std::vector<size_t>& extra_params_idxs =
            camera->ExtraParamsIdxs();
        camera_params_const.insert(camera_params_const.end(),
                                   extra_params_idxs.begin(),
                                   extra_params_idxs.end());
      }

      if (camera_params_const.size() == camera->NumParams()) {
        problem.SetParameterBlockConstant(camera->ParamsData());
      } else {
        SetSubsetManifold(static_cast<int>(camera->NumParams()), camera_params_const, &problem, camera->ParamsData());
        
      }
    }
  }

  ceres::Solver::Options solver_options;
  solver_options.gradient_tolerance = options.gradient_tolerance;
  solver_options.max_num_iterations = options.max_num_iterations;
  solver_options.linear_solver_type = ceres::DENSE_QR;

  // The overhead of creating threads is too large.
  solver_options.num_threads = 1;
#if CERES_VERSION_MAJOR < 2
  solver_options.num_linear_solver_threads = 1;
#endif  // CERES_VERSION_MAJOR

  ceres::Solver::Summary summary;
  ceres::Solve(solver_options, &problem, &summary);

  if (solver_options.minimizer_progress_to_stdout) {
    std::cout << std::endl;
  }

  if (options.print_summary) {
    PrintHeading2("Pose refinement report");
    PrintSolverSummary(summary);
  }

  return summary.IsSolutionUsable();
}

bool RefineRelativePose(const ceres::Solver::Options& options,
                        const std::vector<Eigen::Vector2d>& points1,
                        const std::vector<Eigen::Vector2d>& points2,
                        Eigen::Vector4d* qvec, Eigen::Vector3d* tvec) {
  CHECK_EQ(points1.size(), points2.size());

  // CostFunction assumes unit quaternions.
  *qvec = NormalizeQuaternion(*qvec);

  const double kMaxL2Error = 1.0;
  ceres::LossFunction* loss_function = new ceres::CauchyLoss(kMaxL2Error);

  ceres::Problem problem;

  for (size_t i = 0; i < points1.size(); ++i) {
    ceres::CostFunction* cost_function =
        RelativePoseCostFunction::Create(points1[i], points2[i]);
    problem.AddResidualBlock(cost_function, loss_function, qvec->data(),
                             tvec->data());
  }

  SetQuaternionManifold(&problem, qvec->data());

  SetSphereManifold<3>(&problem, tvec->data());

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  return summary.IsSolutionUsable();
}

}  // namespace colmap
