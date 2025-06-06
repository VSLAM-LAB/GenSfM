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

#ifndef COLMAP_SRC_BASE_COST_FUNCTIONS_H_
#define COLMAP_SRC_BASE_COST_FUNCTIONS_H_

#include <Eigen/Core>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <ceres/jet.h>
#include <type_traits>
#include "base/spline.h"
#include <fstream>
#include <iterator>
#include <iostream>

namespace colmap {


  template <typename CameraModel>
  class BundleAdjustmentCostFunction {
  public:
    explicit BundleAdjustmentCostFunction(const Eigen::Vector2d& point2D)
      : observed_x_(point2D(0)), observed_y_(point2D(1)) {    
}

    static ceres::CostFunction* Create(const Eigen::Vector2d& point2D) {
      return (new ceres::AutoDiffCostFunction<
        BundleAdjustmentCostFunction<CameraModel>, 2, 4, 3, 3,
        CameraModel::kNumParams>(
          new BundleAdjustmentCostFunction(point2D)));
    }

    template <typename T>
    bool operator()(const T* const qvec, const T* const tvec,
      const T* const point3D, const T* const camera_params,
      T* residuals) const {
      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += tvec[0];
      projection[1] += tvec[1];
      projection[2] += tvec[2];

      // Project to image plane.
      projection[0] /= projection[2];
      projection[1] /= projection[2];

      // Distort and transform to pixel space.
      CameraModel::WorldToImage(camera_params, projection[0], projection[1],
        &residuals[0], &residuals[1]);

      // Re-projection error.
      residuals[0] -= T(observed_x_);
      residuals[1] -= T(observed_y_);

      return true;
    }

  private:
    const double observed_x_;
    const double observed_y_;
  };

  // 1D Radial bundle adjustment cost function for variable
  // camera pose and calibration and point parameters.
  template <>
  class BundleAdjustmentCostFunction<Radial1DCameraModel> {
  public:
    explicit BundleAdjustmentCostFunction(const Eigen::Vector2d& point2D)
      : observed_x_(point2D(0)), observed_y_(point2D(1)) {
    }

    static ceres::CostFunction* Create(const Eigen::Vector2d& point2D) {
      return (new ceres::AutoDiffCostFunction<
        BundleAdjustmentCostFunction<Radial1DCameraModel>, 2, 4, 3, 3,
        // Radial1DCameraModel::kNumParams>(
        ImplicitDistortionModel::kNumParams>(
          new BundleAdjustmentCostFunction(point2D)));
    }

    template <typename T>
    bool operator()(const T* const qvec, const T* const tvec,
      const T* const point3D, const T* const camera_params,
      T* residuals) const {
      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += tvec[0];
      projection[1] += tvec[1];

      T x_c, y_c;

      // Subtract principal point from image point
      Radial1DCameraModel::ImageToWorld(camera_params, T(observed_x_),
        T(observed_y_), &x_c, &y_c);

      // Compute radial reprojection error
      T dot_product = projection[0] * x_c + projection[1] * y_c;
      T alpha = dot_product /
        (projection[0] * projection[0] + projection[1] * projection[1]);

      // Re-projection error.
      residuals[0] = alpha * projection[0] - x_c;
      residuals[1] = alpha * projection[1] - y_c;

      return true;
    }

  private:
    const double observed_x_;
    const double observed_y_;
  };


  template <>
  class BundleAdjustmentCostFunction<ImplicitDistortionModel> {
  public:
    explicit BundleAdjustmentCostFunction(const Eigen::Vector2d& point2D)
      : observed_x_(point2D(0)), observed_y_(point2D(1)) {
    }

    static ceres::CostFunction* Create(const Eigen::Vector2d& point2D) {
      return (new ceres::AutoDiffCostFunction<
        BundleAdjustmentCostFunction<ImplicitDistortionModel>, 2, 4, 3, 3,
        ImplicitDistortionModel::kNumParams>(
          new BundleAdjustmentCostFunction(point2D)));
    }

    template <typename T>
    bool operator()(const T* const qvec, const T* const tvec,
      const T* const point3D, const T* const camera_params,
      T* residuals) const {
      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += tvec[0];
      projection[1] += tvec[1];
      projection[2] += tvec[2];

      projection[0] /= projection[2];
      projection[1] /= projection[2];

      T x_c, y_c;

      int num_control_points = (ImplicitDistortionModel::kNumParams - 2) / 2;
      std::vector<T> sample_x = {};
      for (int i = 2; i < 2 + num_control_points; i++) {
        sample_x.push_back(camera_params[i]);
      }
      std::vector<T> sample_y = {};
      for (int i = 2 + num_control_points; i < ImplicitDistortionModel::kNumParams; i++) {
        sample_y.push_back(camera_params[i]);
      }

      tk::spline<T> spline_focal_lengths;
      spline_focal_lengths.set_points(sample_x, sample_y);

      // Subtract principal point from image point
      ImplicitDistortionModel::ImageToWorld(camera_params, T(observed_x_),
        T(observed_y_), &x_c, &y_c);

      if (sample_x[0] != T(350)) {

        T radius = sqrt(x_c * x_c + y_c * y_c);

        T rho = sqrt(projection[0] * projection[2] * projection[0] * projection[2] + projection[1] * projection[2] * projection[1] * projection[2]);
        T theta = ceres::atan2(rho, projection[2]);
        T r_calculated = spline_focal_lengths(theta);
        T focal_length = r_calculated / tan(theta);
        if (theta >= sample_x[0] && theta <= sample_x[sample_x.size() - 1]) {
          residuals[0] = projection[0] * T(focal_length) + T(camera_params[0]) - T(observed_x_);
          residuals[1] = projection[1] * T(focal_length) + T(camera_params[1]) - T(observed_y_);
        }
        else {
          projection[0] *= projection[2];
          projection[1] *= projection[2];
          T dot_product = projection[0] * x_c + projection[1] * y_c;
          T alpha = dot_product /
            (projection[0] * projection[0] + projection[1] * projection[1]);

          // Re-projection error.

          residuals[0] = alpha * projection[0] - x_c;
          residuals[1] = alpha * projection[1] - y_c;
        }
      }
      else {
        //compute radial reproj error
        projection[0] *= projection[2];
        projection[1] *= projection[2];
        T dot_product = projection[0] * x_c + projection[1] * y_c;
        T alpha = dot_product /
          (projection[0] * projection[0] + projection[1] * projection[1]);

        // Re-projection error.
        residuals[0] = alpha * projection[0] - x_c;
        residuals[1] = alpha * projection[1] - y_c;

      }

      return true;
    }

  private:
    const double observed_x_;
    const double observed_y_;
    // const std::vector<double> camera_parameters_;
  };

  // Bundle adjustment cost function for variable
  // camera calibration and point parameters, and fixed camera pose.

  template <typename CameraModel>
  class BundleAdjustmentConstantPoseCostFunction {
  public:
    BundleAdjustmentConstantPoseCostFunction(const Eigen::Vector4d& qvec,
      const Eigen::Vector3d& tvec,
      const Eigen::Vector2d& point2D)
      : qw_(qvec(0)),
      qx_(qvec(1)),
      qy_(qvec(2)),
      qz_(qvec(3)),
      tx_(tvec(0)),
      ty_(tvec(1)),
      tz_(tvec(2)),
      observed_x_(point2D(0)),
      observed_y_(point2D(1)) {
    }

    static ceres::CostFunction* Create(const Eigen::Vector4d& qvec,
      const Eigen::Vector3d& tvec,
      const Eigen::Vector2d& point2D) {
      return (new ceres::AutoDiffCostFunction<
        BundleAdjustmentConstantPoseCostFunction<CameraModel>, 2, 3,
        CameraModel::kNumParams>(
          new BundleAdjustmentConstantPoseCostFunction(qvec, tvec, point2D)));
    }

    template <typename T>
    bool operator()(const T* const point3D, const T* const camera_params,
      T* residuals) const {
      const T qvec[4] = { T(qw_), T(qx_), T(qy_), T(qz_) };

      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += T(tx_);
      projection[1] += T(ty_);
      projection[2] += T(tz_);

      // Project to image plane.
      projection[0] /= projection[2];
      projection[1] /= projection[2];

      // Distort and transform to pixel space.
      CameraModel::WorldToImage(camera_params, projection[0], projection[1],
        &residuals[0], &residuals[1]);

      // Re-projection error.
      residuals[0] -= T(observed_x_);
      residuals[1] -= T(observed_y_);

      return true;
    }

  private:
    const double qw_;
    const double qx_;
    const double qy_;
    const double qz_;
    const double tx_;
    const double ty_;
    const double tz_;
    const double observed_x_;
    const double observed_y_;
  };

  // Bundle adjustment cost function for variable
  // camera calibration and point parameters, and fixed camera pose.
  // Specialized for the 1D radial camera model
  template <>
  class BundleAdjustmentConstantPoseCostFunction<Radial1DCameraModel> {
  public:
    BundleAdjustmentConstantPoseCostFunction(const Eigen::Vector4d& qvec,
      const Eigen::Vector3d& tvec,
      const Eigen::Vector2d& point2D)
      : qw_(qvec(0)),
      qx_(qvec(1)),
      qy_(qvec(2)),
      qz_(qvec(3)),
      tx_(tvec(0)),
      ty_(tvec(1)),
      tz_(tvec(2)),
      observed_x_(point2D(0)),
      observed_y_(point2D(1)) {    
}

    static ceres::CostFunction* Create(const Eigen::Vector4d& qvec,
      const Eigen::Vector3d& tvec,
      const Eigen::Vector2d& point2D) {
      return (new ceres::AutoDiffCostFunction<
        BundleAdjustmentConstantPoseCostFunction<Radial1DCameraModel>, 2, 3,
        // Radial1DCameraModel::kNumParams>(
        ImplicitDistortionModel::kNumParams>(
          new BundleAdjustmentConstantPoseCostFunction(qvec, tvec, point2D)));
    }

    template <typename T>
    bool operator()(const T* const point3D, const T* const camera_params,
      T* residuals) const {
      const T qvec[4] = { T(qw_), T(qx_), T(qy_), T(qz_) };

      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += T(tx_);
      projection[1] += T(ty_);

      T x_c, y_c;

      // Subtract principal point from image point
      Radial1DCameraModel::ImageToWorld(camera_params, T(observed_x_),
        T(observed_y_), &x_c, &y_c);

      // Compute radial reprojection error
      T dot_product = projection[0] * x_c + projection[1] * y_c;
      T alpha = dot_product /
        (projection[0] * projection[0] + projection[1] * projection[1]);

      // Re-projection error.
      residuals[0] = alpha * projection[0] - x_c;
      residuals[1] = alpha * projection[1] - y_c;

      return true;
    }

  private:
    const double qw_;
    const double qx_;
    const double qy_;
    const double qz_;
    const double tx_;
    const double ty_;
    const double tz_;
    const double observed_x_;
    const double observed_y_;
  };

  template <>
  class BundleAdjustmentConstantPoseCostFunction<ImplicitDistortionModel> {
  public:
    BundleAdjustmentConstantPoseCostFunction(const Eigen::Vector4d& qvec,
      const Eigen::Vector3d& tvec,
      const Eigen::Vector2d& point2D)
      : qw_(qvec(0)),
      qx_(qvec(1)),
      qy_(qvec(2)),
      qz_(qvec(3)),
      tx_(tvec(0)),
      ty_(tvec(1)),
      tz_(tvec(2)),
      observed_x_(point2D(0)),
      observed_y_(point2D(1)) {
    }

    static ceres::CostFunction* Create(const Eigen::Vector4d& qvec,
      const Eigen::Vector3d& tvec,
      const Eigen::Vector2d& point2D) {
      return (new ceres::AutoDiffCostFunction<
        BundleAdjustmentConstantPoseCostFunction<ImplicitDistortionModel>, 2, 3,
        ImplicitDistortionModel::kNumParams>(
          new BundleAdjustmentConstantPoseCostFunction(qvec, tvec, point2D)));
    }

    template <typename T>
    bool operator()(const T* const point3D, const T* const camera_params,
      T* residuals) const {
      const T qvec[4] = { T(qw_), T(qx_), T(qy_), T(qz_) };

      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += T(tx_);
      projection[1] += T(ty_);
      projection[2] += T(tz_);
      // Project to image plane 
      projection[0] /= projection[2];
      projection[1] /= projection[2];

      T x_c, y_c;
      int num_control_points = (ImplicitDistortionModel::kNumParams - 2) / 2;
      std::vector<T> sample_x = {};
      for (int i = 2; i < num_control_points + 2; i++) {
        sample_x.push_back(camera_params[i]);
      }

      std::vector<T> sample_y = {};
      for (int i = num_control_points + 2; i < ImplicitDistortionModel::kNumParams; i++) {
        sample_y.push_back(camera_params[i]);
      }
      tk::spline<T> spline_focal_lengths;
      spline_focal_lengths.set_points(sample_x, sample_y);
      ImplicitDistortionModel::ImageToWorld(camera_params, T(observed_x_),
        T(observed_y_), &x_c, &y_c);
      if (sample_x[0] != T(350)) {
        T radius = sqrt(x_c * x_c + y_c * y_c);
        T rho = sqrt(projection[0] * projection[2] * projection[0] * projection[2] + projection[1] * projection[2] * projection[1] * projection[2]);
        T theta = ceres::atan2(rho, projection[2]);
        T r_calculated = spline_focal_lengths(theta);
        T focal_length = r_calculated / tan(theta);
        if (theta >= sample_x[0] && theta <= sample_x[sample_x.size() - 1]) {
          residuals[0] = projection[0] * T(focal_length) + camera_params[0] - T(observed_x_);
          residuals[1] = projection[1] * T(focal_length) + camera_params[1] - T(observed_y_);
        }
        else {
          projection[0] *= projection[2];
          projection[1] *= projection[2];
          T dot_product = projection[0] * x_c + projection[1] * y_c;
          T alpha = dot_product /
            (projection[0] * projection[0] + projection[1] * projection[1]);

          // Re-projection error.
          residuals[0] = alpha * projection[0] - x_c;
          residuals[1] = alpha * projection[1] - y_c;
        }
      }
      else {
        projection[0] *= projection[2];
        projection[1] *= projection[2];
        //compute radial reproj error
        T dot_product = projection[0] * x_c + projection[1] * y_c;
        T alpha = dot_product /
          (projection[0] * projection[0] + projection[1] * projection[1]);
        // Re-projection error.
        residuals[0] = alpha * projection[0] - x_c;
        residuals[1] = alpha * projection[1] - y_c;
      }

      return true;
    }

  private:
    const double qw_;
    const double qx_;
    const double qy_;
    const double qz_;
    const double tx_;
    const double ty_;
    const double tz_;
    const double observed_x_;
    const double observed_y_;
  };

  // Rig bundle adjustment cost function for variable camera pose and calibration
  // and point parameters. Different from the standard bundle adjustment function,
  // this cost function is suitable for camera rigs with consistent relative poses
  // of the cameras within the rig. The cost function first projects points into
  // the local system of the camera rig and then into the local system of the
  // camera within the rig.
  template <typename CameraModel>
  class RigBundleAdjustmentCostFunction {
  public:
    explicit RigBundleAdjustmentCostFunction(const Eigen::Vector2d& point2D)
      : observed_x_(point2D(0)), observed_y_(point2D(1)) {
    }

    static ceres::CostFunction* Create(const Eigen::Vector2d& point2D) {
      return (new ceres::AutoDiffCostFunction<
        RigBundleAdjustmentCostFunction<CameraModel>, 2, 4, 3, 4, 3, 3,
        CameraModel::kNumParams>(
          new RigBundleAdjustmentCostFunction(point2D)));
    }

    template <typename T>
    bool operator()(const T* const rig_qvec, const T* const rig_tvec,
      const T* const rel_qvec, const T* const rel_tvec,
      const T* const point3D, const T* const camera_params,
      T* residuals) const {
      // Concatenate rotations.
      T qvec[4];
      ceres::QuaternionProduct(rel_qvec, rig_qvec, qvec);

      // Concatenate translations.
      T tvec[3];
      ceres::UnitQuaternionRotatePoint(rel_qvec, rig_tvec, tvec);
      tvec[0] += rel_tvec[0];
      tvec[1] += rel_tvec[1];
      tvec[2] += rel_tvec[2];

      // Rotate and translate.
      T projection[3];
      ceres::UnitQuaternionRotatePoint(qvec, point3D, projection);
      projection[0] += tvec[0];
      projection[1] += tvec[1];
      projection[2] += tvec[2];

      // Project to image plane.
      projection[0] /= projection[2];
      projection[1] /= projection[2];

      // Distort and transform to pixel space.
      CameraModel::WorldToImage(camera_params, projection[0], projection[1],
        &residuals[0], &residuals[1]);

      // Re-projection error.
      residuals[0] -= T(observed_x_);
      residuals[1] -= T(observed_y_);

      return true;
    }

  private:
    const double observed_x_;
    const double observed_y_;
  };

  // Cost function for refining two-view geometry based on the Sampson-Error.
  //
  // First pose is assumed to be located at the origin with 0 rotation. Second
  // pose is assumed to be on the unit sphere around the first pose, i.e. the
  // pose of the second camera is parameterized by a 3D rotation and a
  // 3D translation with unit norm. `tvec` is therefore over-parameterized as is
  // and should be down-projected using `HomogeneousVectorParameterization`.
  class RelativePoseCostFunction {
  public:
    RelativePoseCostFunction(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2)
      : x1_(x1(0)), y1_(x1(1)), x2_(x2(0)), y2_(x2(1)) {    
}

    static ceres::CostFunction* Create(const Eigen::Vector2d& x1,
      const Eigen::Vector2d& x2) {
      return (new ceres::AutoDiffCostFunction<RelativePoseCostFunction, 1, 4, 3>(
        new RelativePoseCostFunction(x1, x2)));
    }

    template <typename T>
    bool operator()(const T* const qvec, const T* const tvec,
      T* residuals) const {
      Eigen::Matrix<T, 3, 3, Eigen::RowMajor> R;
      ceres::QuaternionToRotation(qvec, R.data());

      // Matrix representation of the cross product t x R.
      Eigen::Matrix<T, 3, 3> t_x;
      t_x << T(0), -tvec[2], tvec[1], tvec[2], T(0), -tvec[0], -tvec[1], tvec[0],
        T(0);

      // Essential matrix.
      const Eigen::Matrix<T, 3, 3> E = t_x * R;

      // Homogeneous image coordinates.
      const Eigen::Matrix<T, 3, 1> x1_h(T(x1_), T(y1_), T(1));
      const Eigen::Matrix<T, 3, 1> x2_h(T(x2_), T(y2_), T(1));

      // Squared sampson error.
      const Eigen::Matrix<T, 3, 1> Ex1 = E * x1_h;
      const Eigen::Matrix<T, 3, 1> Etx2 = E.transpose() * x2_h;
      const T x2tEx1 = x2_h.transpose() * Ex1;
      residuals[0] = x2tEx1 * x2tEx1 /
        (Ex1(0) * Ex1(0) + Ex1(1) * Ex1(1) + Etx2(0) * Etx2(0) +
          Etx2(1) * Etx2(1));

      return true;
    }

  private:
    const double x1_;
    const double y1_;
    const double x2_;
    const double y2_;
  };

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_COST_FUNCTIONS_H_