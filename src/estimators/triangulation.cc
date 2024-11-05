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

#include "estimators/triangulation.h"

#include <Eigen/Geometry>

#include "base/projection.h"
#include "base/triangulation.h"
#include "base/camera_models.h"
#include "estimators/essential_matrix.h"
#include "optim/combination_sampler.h"
#include "optim/loransac.h"
#include "util/logging.h"
#include "util/math.h"

namespace colmap {

void TriangulationEstimator::SetMinTriAngle(const double min_tri_angle) {
  CHECK_GE(min_tri_angle, 0);
  min_tri_angle_ = min_tri_angle;
}

void TriangulationEstimator::SetResidualType(const ResidualType residual_type) {
  residual_type_ = residual_type;
}




std::vector<TriangulationEstimator::M_t> TriangulationEstimator::Estimate(
    const std::vector<X_t>& point_data,
    const std::vector<Y_t>& pose_data, bool initial) const {
  CHECK_GE(point_data.size(), 2);
  CHECK_EQ(point_data.size(), pose_data.size());

  // Check if any of the cameras in this sample is a 1D radial camera
  bool have_radial = false;
  std::vector<bool> is_radial;
  int radial_count = 0;
  for(const auto &p : pose_data) {
    if(p.camera->ModelId() == Radial1DCameraModel::model_id || p.camera->ModelId() == ImplicitDistortionModel::model_id){
      is_radial.push_back(true);
      radial_count++;
      have_radial = true;
      // full triangulation
      // have_radial = false;
      // break;
    } else {
      is_radial.push_back(false);
    }
    
  }
  // print the value of have_radial
  // std::cout << "Initial in Triangulation Estimator: " << initial << std::endl;
  // if (!initial) {
  //   have_radial = false;
  // }

  std::vector<Eigen::Vector3d> candidates;
  if(!have_radial) { // standard point-based triangulation

    if(point_data.size() == 3) {
      // This is a minimal sample and not LO-step in RANSAC
      // We run all possible combinations of three points

      candidates.push_back(TriangulatePoint(
        pose_data[0].proj_matrix, pose_data[1].proj_matrix,
        point_data[0].point_normalized, point_data[1].point_normalized));  
      candidates.push_back(TriangulatePoint(
        pose_data[0].proj_matrix, pose_data[2].proj_matrix,
        point_data[0].point_normalized, point_data[2].point_normalized));
      candidates.push_back(TriangulatePoint(
        pose_data[1].proj_matrix, pose_data[2].proj_matrix,
        point_data[1].point_normalized, point_data[2].point_normalized));
    
    } else {
      // This is a non-minimal sample (for refinement in RANSAC)
      std::vector<Eigen::Matrix3x4d> proj_matrices;
      proj_matrices.reserve(point_data.size());
      std::vector<Eigen::Vector2d> points;
      points.reserve(point_data.size());
      for (size_t i = 0; i < point_data.size(); ++i) {
        proj_matrices.push_back(pose_data[i].proj_matrix);
        points.push_back(point_data[i].point_normalized);
      }
      candidates.push_back(TriangulateMultiViewPoint(proj_matrices, points));
    }

  } else { // triangulation with radial cameras involved
    // If we have only two points and two radial cameras, we can not triangulate
    if (point_data.size() == 2 && radial_count == 2)
      return candidates;
    // For minimal mixed cases, we run all possible combinations of three points if mixed case is involved
    if (point_data.size() == 3 && radial_count < 3) {
      for (size_t i = 0; i < 3; i++) {
        std::vector<Eigen::Vector3d> lines;
        std::vector<Eigen::Matrix3x4d> proj_matrices;
        for (size_t k = 0; k < 2; k++) {
          proj_matrices.push_back(pose_data[i].proj_matrix);
          int curr_idx = (i + 1 + k) % 3;
          if (is_radial[curr_idx]) {
            // add line corresponding to the radial line
            lines.emplace_back(point_data[curr_idx].point_normalized(1), -point_data[curr_idx].point_normalized(0), 0.0);
          } else {        
            // add two lines corresponding to x/y coordinates
            lines.emplace_back(1.0, 0.0, -point_data[curr_idx].point_normalized(0));
            proj_matrices.push_back(pose_data[curr_idx].proj_matrix);        
            lines.emplace_back(0.0, 1.0, -point_data[curr_idx].point_normalized(1));
          }
        }
        if (lines.size() == 4) {
          candidates.push_back(TriangulateMultiViewPointFromLines(proj_matrices, lines));
        }
      }
    }
    std::vector<Eigen::Vector3d> lines;
    std::vector<Eigen::Matrix3x4d> proj_matrices;

    for (size_t i = 0; i < point_data.size(); ++i) {
      proj_matrices.push_back(pose_data[i].proj_matrix);
      if (is_radial[i]) {
        // add line corresponding to the radial line
        lines.emplace_back(point_data[i].point_normalized(1), -point_data[i].point_normalized(0), 0.0);
      } else {        
        // add two lines corresponding to x/y coordinates
        lines.emplace_back(1.0, 0.0, -point_data[i].point_normalized(0));
        proj_matrices.push_back(pose_data[i].proj_matrix);        
        lines.emplace_back(0.0, 1.0, -point_data[i].point_normalized(1));
      }        
    }

    candidates.push_back(TriangulateMultiViewPointFromLines(proj_matrices, lines));

  }

  std::vector<Eigen::Vector3d> output;

  for(Eigen::Vector3d &xyz : candidates) {
    // check cheirality for each camera (or half-plane constraint for radial cameras)
    bool cheiral_ok = true;
    for (size_t i = 0; i < pose_data.size(); ++i) {
      if (is_radial[i]) {
        Eigen::Vector2d n = pose_data[i].proj_matrix.topRows<2>() * xyz.homogeneous();
        cheiral_ok &= n.dot(point_data[i].point_normalized) > 0;
      } else {
        cheiral_ok &= HasPointPositiveDepth(pose_data[i].proj_matrix, xyz);
      }
    }

    if(!cheiral_ok)
      continue;

    bool tri_angle_ok = false;
    // check point-based triangulation angle
    for (size_t i = 0; i < pose_data.size(); ++i) {
      if (is_radial[i]) {
        continue;
      }
      for (size_t j = 0; j < i; ++j) {
        if(is_radial[j]){
          continue;
        }
        const double tri_angle = CalculateTriangulationAngle(
            pose_data[i].proj_center, pose_data[j].proj_center, xyz);
        tri_angle_ok |= (tri_angle >= min_tri_angle_);
      }
    }

    if(!have_radial && !tri_angle_ok) {
      // we have only pinhole-like cameras and poor triangulatiopn angle
      continue;
    }

    // TODO: triangulation angle equivalent for radial cameras

    output.push_back(xyz);
  }

  return output;
}

std::vector<TriangulationEstimator::M_t> TriangulationEstimator::EstimateStandard(
  
    const std::vector<X_t>& point_data,
    const std::vector<Y_t>& pose_data, bool initial) const {
  CHECK_GE(point_data.size(), 2);
  CHECK_EQ(point_data.size(), pose_data.size());
  // std::cout << "==============================Starting standard triangulation ===============================" << std::endl;

  

  std::vector<Eigen::Vector3d> candidates;
   // standard point-based triangulation

  if(point_data.size() == 3) {
    // This is a minimal sample and not LO-step in RANSAC
    // We run all possible combinations of three points

    candidates.push_back(TriangulatePoint(
      pose_data[0].proj_matrix, pose_data[1].proj_matrix,
      point_data[0].point_normalized_standard, point_data[1].point_normalized_standard));  
    candidates.push_back(TriangulatePoint(
      pose_data[0].proj_matrix, pose_data[2].proj_matrix,
      point_data[0].point_normalized_standard, point_data[2].point_normalized_standard));
    candidates.push_back(TriangulatePoint(
      pose_data[1].proj_matrix, pose_data[2].proj_matrix,
      point_data[1].point_normalized_standard, point_data[2].point_normalized_standard));
  
  } else {
    // This is a non-minimal sample (for refinement in RANSAC)
    std::vector<Eigen::Matrix3x4d> proj_matrices;
    proj_matrices.reserve(point_data.size());
    std::vector<Eigen::Vector2d> points;
    points.reserve(point_data.size());
    for (size_t i = 0; i < point_data.size(); ++i) {
      proj_matrices.push_back(pose_data[i].proj_matrix);
      points.push_back(point_data[i].point_normalized_standard);
    }
    candidates.push_back(TriangulateMultiViewPoint(proj_matrices, points));
  }



  std::vector<Eigen::Vector3d> output;

  for(Eigen::Vector3d &xyz : candidates) {
    // check cheirality for each camera (or half-plane constraint for radial cameras)
    bool cheiral_ok = true;
    for (size_t i = 0; i < pose_data.size(); ++i) {
      if(pose_data[i].camera->ModelId() == Radial1DCameraModel::model_id || pose_data[i].camera->ModelId() == ImplicitDistortionModel::model_id) {
        Eigen::Vector2d n = pose_data[i].proj_matrix.topRows<2>() * xyz.homogeneous();
        cheiral_ok &= n.dot(point_data[i].point_normalized_standard) > 0;
        // cheiral_ok &= n.dot(point_data[i].point_normalized) > 0;
      } 
      else {
      cheiral_ok &= HasPointPositiveDepth(pose_data[i].proj_matrix_standard, xyz);
      // cheiral_ok &= HasPointPositiveDepth(pose_data[i].proj_matrix, xyz);
      }
    }

    if(!cheiral_ok)
      continue;

    bool tri_angle_ok = false;
    // check point-based triangulation angle
    for (size_t i = 0; i < pose_data.size(); ++i) {
      // if(pose_data[i].camera->ModelId() == Radial1DCameraModel::model_id) {
      //   continue;
      // }
      for (size_t j = 0; j < i; ++j) {
        // if(pose_data[j].camera->ModelId() == Radial1DCameraModel::model_id) {
        //   continue;
        // }
        const double tri_angle = CalculateTriangulationAngle(
            pose_data[i].proj_center, pose_data[j].proj_center, xyz);
        tri_angle_ok |= (tri_angle >= min_tri_angle_);
      }
    }
    if(!tri_angle_ok) {
      // we have only pinhole-like cameras and poor triangulatiopn angle
      continue;
    }
    // TODO: triangulation angle equivalent for radial cameras
    output.push_back(xyz);
  }

  return output;
}





std::vector<TriangulationEstimator::M_t> TriangulationEstimator::EstimateInitial(
    const std::vector<X_t>& point_data,
    const std::vector<Y_t>& pose_data) const {
  CHECK_GE(point_data.size(), 2);
  CHECK_EQ(point_data.size(), pose_data.size());

  // Check if any of the cameras in this sample is a 1D radial camera
  bool have_radial = false;
  for(const auto &p : pose_data) { 
    if(p.camera->ModelId() == Radial1DCameraModel::model_id || p.camera->ModelId() == ImplicitDistortionModel::model_id){
      have_radial = true;
      // full triangulation
      // have_radial = false;
      break;
    }
  }

  std::vector<Eigen::Vector3d> candidates;
  if(!have_radial) { // standard point-based triangulation

    if(point_data.size() == 3) {
      // This is a minimal sample and not LO-step in RANSAC
      // We run all possible combinations of three points

      candidates.push_back(TriangulatePoint(
        pose_data[0].proj_matrix, pose_data[1].proj_matrix,
        point_data[0].point_normalized, point_data[1].point_normalized));  
      candidates.push_back(TriangulatePoint(
        pose_data[0].proj_matrix, pose_data[2].proj_matrix,
        point_data[0].point_normalized, point_data[2].point_normalized));
      candidates.push_back(TriangulatePoint(
        pose_data[1].proj_matrix, pose_data[2].proj_matrix,
        point_data[1].point_normalized, point_data[2].point_normalized));
    
    } else {
      // This is a non-minimal sample (for refinement in RANSAC)
      std::vector<Eigen::Matrix3x4d> proj_matrices;
      proj_matrices.reserve(point_data.size());
      std::vector<Eigen::Vector2d> points;
      points.reserve(point_data.size());
      for (size_t i = 0; i < point_data.size(); ++i) {
        proj_matrices.push_back(pose_data[i].proj_matrix);
        points.push_back(point_data[i].point_normalized);
      }
      candidates.push_back(TriangulateMultiViewPoint(proj_matrices, points));
    }

  } else { // triangulation with radial cameras involved

    std::vector<Eigen::Vector3d> lines;
    std::vector<Eigen::Matrix3x4d> proj_matrices;

    for (size_t i = 0; i < point_data.size(); ++i) {
      proj_matrices.push_back(pose_data[i].proj_matrix);
      if(pose_data[i].camera->ModelId() == Radial1DCameraModel::model_id || pose_data[i].camera->ModelId() == ImplicitDistortionModel::model_id){
        // add line corresponding to the radial line
        lines.emplace_back(point_data[i].point_normalized(1), -point_data[i].point_normalized(0), 0.0);
      } else {        
        // add two lines corresponding to x/y coordinates
        lines.emplace_back(1.0, 0.0, -point_data[i].point_normalized(0));
        proj_matrices.push_back(pose_data[i].proj_matrix);        
        lines.emplace_back(0.0, 1.0, -point_data[i].point_normalized(1));
      }        
    }

    candidates.push_back(TriangulateMultiViewPointFromLines(proj_matrices, lines));

  }

  std::vector<Eigen::Vector3d> output;

  for(Eigen::Vector3d &xyz : candidates) {
    // check cheirality for each camera (or half-plane constraint for radial cameras)
    bool cheiral_ok = true;
    for (size_t i = 0; i < pose_data.size(); ++i) {
      if(pose_data[i].camera->ModelId() == Radial1DCameraModel::model_id || pose_data[i].camera->ModelId() == ImplicitDistortionModel::model_id){
        Eigen::Vector2d n = pose_data[i].proj_matrix.topRows<2>() * xyz.homogeneous();
        cheiral_ok &= n.dot(point_data[i].point_normalized) > 0;
      } else {
        cheiral_ok &= HasPointPositiveDepth(pose_data[i].proj_matrix, xyz);
      }
    }

    if(!cheiral_ok)
      continue;

    bool tri_angle_ok = false;
    // check point-based triangulation angle
    for (size_t i = 0; i < pose_data.size(); ++i) {
      if(pose_data[i].camera->ModelId() == Radial1DCameraModel::model_id || pose_data[i].camera->ModelId() == ImplicitDistortionModel::model_id) {
        continue;
      }
      for (size_t j = 0; j < i; ++j) {
        if(pose_data[j].camera->ModelId() == Radial1DCameraModel::model_id || pose_data[j].camera->ModelId() == ImplicitDistortionModel::model_id){
          continue;
        }
        const double tri_angle = CalculateTriangulationAngle(
            pose_data[i].proj_center, pose_data[j].proj_center, xyz);
        tri_angle_ok |= (tri_angle >= min_tri_angle_);
      }
    }

    if(!have_radial && !tri_angle_ok) {
      // we have only pinhole-like cameras and poor triangulatiopn angle
      continue;
    }

    // TODO: triangulation angle equivalent for radial cameras

    output.push_back(xyz);
  }

  return output;
}

void TriangulationEstimator::Residuals(const std::vector<X_t>& point_data,
                                       const std::vector<Y_t>& pose_data,
                                       const M_t& xyz,
                                       std::vector<double>* residuals) const {
  CHECK_EQ(point_data.size(), pose_data.size());

  residuals->resize(point_data.size());

  for (size_t i = 0; i < point_data.size(); ++i) {
    if (residual_type_ == ResidualType::REPROJECTION_ERROR) {
      (*residuals)[i] = CalculateSquaredReprojectionError(
          point_data[i].point, xyz, pose_data[i].proj_matrix,
          *pose_data[i].camera);
    } else if (residual_type_ == ResidualType::ANGULAR_ERROR) {
      const double angular_error = CalculateAngularError(
          point_data[i].point, xyz, pose_data[i].proj_matrix,
          *pose_data[i].camera);
      (*residuals)[i] = angular_error * angular_error;
    } else if(residual_type_== ResidualType::ANGULAR_ERROR_SPLITTING){
      const double angular_error = CalculateAngularErrorSplitting(
          point_data[i].point, xyz, pose_data[i].proj_matrix, point_data[i].focal_length,
          *pose_data[i].camera);
      (*residuals)[i] = angular_error * angular_error;
    } else {
      LOG(FATAL) << "Invalid residual type";

    }
  }
}

bool EstimateTriangulation(
    const EstimateTriangulationOptions& options,
    const std::vector<TriangulationEstimator::PointData>& point_data,
    const std::vector<TriangulationEstimator::PoseData>& pose_data,
    std::vector<char>* inlier_mask, Eigen::Vector3d* xyz, bool initial, bool standard_triangulation) {
  CHECK_NOTNULL(inlier_mask);
  CHECK_NOTNULL(xyz);
  CHECK_GE(point_data.size(), 2);
  CHECK_EQ(point_data.size(), pose_data.size());
  options.Check();
  if (point_data.size() == 2) {
    TriangulationEstimator estimator;
    std::vector<Eigen::Vector3d> xyzs = estimator.Estimate(point_data, pose_data, initial);
    if (xyzs.empty()) {
      return false;
    }
    else {
      *xyz = xyzs[0];
      inlier_mask->resize(point_data.size(), true);
      return true;
    }
  }
  // std::coust << "Initial in EstimateTriangulation: " << initial << std::endl;

  // Robustly estimate track using LORANSAC.
  LORANSAC<TriangulationEstimator, TriangulationEstimator,
           InlierSupportMeasurer, CombinationSampler>
      ransac(options.ransac_options);
  ransac.estimator.SetMinTriAngle(options.min_tri_angle);
  ransac.estimator.SetResidualType(options.residual_type);
  ransac.local_estimator.SetMinTriAngle(options.min_tri_angle);
  ransac.local_estimator.SetResidualType(options.residual_type);
  if(options.residual_type == TriangulationEstimator::ResidualType::ANGULAR_ERROR_SPLITTING) {
    initial = true;
  }
   
  // decltype(auto) report = standard_triangulation ?
  //                       ransac.EstimateStandard(point_data, pose_data, initial) :
  //                       ransac.Estimate(point_data, pose_data, initial);
  decltype(auto) report = ransac.Estimate(point_data, pose_data, initial);
  if (!report.success) {
    return false;
  }

  *inlier_mask = report.inlier_mask;
  *xyz = report.model;

  return report.success;
}


}  // namespace colmap
