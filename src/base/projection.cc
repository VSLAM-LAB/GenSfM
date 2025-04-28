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

#include "base/projection.h"

#include "base/pose.h"
#include "base/camera_models.h"
#include "util/matrix.h"


namespace colmap {

  Eigen::Matrix3x4d ComposeProjectionMatrix(const Eigen::Vector4d& qvec,
    const Eigen::Vector3d& tvec) {
    Eigen::Matrix3x4d proj_matrix;
    proj_matrix.leftCols<3>() = QuaternionToRotationMatrix(qvec);
    proj_matrix.rightCols<1>() = tvec;
    return proj_matrix;
  }

  Eigen::Matrix3x4d ComposeProjectionMatrix(const Eigen::Matrix3d& R,
    const Eigen::Vector3d& T) {
    Eigen::Matrix3x4d proj_matrix;
    proj_matrix.leftCols<3>() = R;
    proj_matrix.rightCols<1>() = T;
    return proj_matrix;
  }

  Eigen::Matrix3x4d InvertProjectionMatrix(const Eigen::Matrix3x4d& proj_matrix) {
    Eigen::Matrix3x4d inv_proj_matrix;
    inv_proj_matrix.leftCols<3>() = proj_matrix.leftCols<3>().transpose();
    inv_proj_matrix.rightCols<1>() = ProjectionCenterFromMatrix(proj_matrix);
    return inv_proj_matrix;
  }

  Eigen::Matrix3d ComputeClosestRotationMatrix(const Eigen::Matrix3d& matrix) {
    const Eigen::JacobiSVD<Eigen::Matrix3d> svd(
      matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixU() * (svd.matrixV().transpose());
    if (R.determinant() < 0.0) {
      R *= -1.0;
    }
    return R;
  }

  bool DecomposeProjectionMatrix(const Eigen::Matrix3x4d& P, Eigen::Matrix3d* K,
    Eigen::Matrix3d* R, Eigen::Vector3d* T) {
    Eigen::Matrix3d RR;
    Eigen::Matrix3d QQ;
    DecomposeMatrixRQ(P.leftCols<3>().eval(), &RR, &QQ);

    *R = ComputeClosestRotationMatrix(QQ);

    const double det_K = RR.determinant();
    if (det_K == 0) {
      return false;
    }
    else if (det_K > 0) {
      *K = RR;
    }
    else {
      *K = -RR;
    }

    for (int i = 0; i < 3; ++i) {
      if ((*K)(i, i) < 0.0) {
        K->col(i) = -K->col(i);
        R->row(i) = -R->row(i);
      }
    }

    *T = K->triangularView<Eigen::Upper>().solve(P.col(3));
    if (det_K < 0) {
      *T = -(*T);
    }

    return true;
  }

  Eigen::Vector2d ProjectPointToImage(const Eigen::Vector3d& point3D,
    const Eigen::Matrix3x4d& proj_matrix,
    const Camera& camera) {
    CHECK(camera.ModelId() != Radial1DCameraModel::model_id &&
      camera.ModelId() != ImplicitDistortionModel::model_id);
    const Eigen::Vector3d world_point = proj_matrix * point3D.homogeneous();
    return camera.WorldToImage(world_point.hnormalized());
  }

  double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D,
    const Eigen::Vector3d& point3D,
    const Eigen::Vector4d& qvec,
    const Eigen::Vector3d& tvec,
    const Camera& camera) {
    const Eigen::Vector3d proj_point3D =
      QuaternionRotatePoint(qvec, point3D) + tvec;

    Eigen::Vector2d proj_point2D;
    if (!camera.IsFullyCalibrated(point2D)) {
      // if(false) {       
      const Eigen::Vector2d n = proj_point3D.topRows<2>().normalized();
      const Eigen::Vector2d point2D_center = camera.ImageToWorld(point2D);
      const double dot_product = n.dot(point2D_center);

      // check that we project onto the correct half-plane
      // if(fabs(dot_product) < std::numeric_limits<double>::epsilon()) {
      // check if we have a focal length that can apply
      if (dot_product < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }
      proj_point2D = camera.WorldToImage(dot_product * n);

    }
    else if (camera.ModelId() == ImplicitDistortionModel::model_id) {
      double focal = camera.EvalFocalLength(point2D);

      // Check that point is infront of camera.
      if (focal * proj_point3D.z() < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }

      proj_point2D = camera.WorldToImage(focal * proj_point3D.hnormalized());
    }
    else {
      // Check that point is infront of camera.
      if (proj_point3D.z() * camera.FocalLength() < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }
      proj_point2D = camera.WorldToImage(proj_point3D.hnormalized());
    }

    return (proj_point2D - point2D).squaredNorm();
  }

  double CalculateSquaredReprojectionErrorFinal(const Eigen::Vector2d& point2D,
    const Eigen::Vector3d& point3D,
    const Eigen::Vector4d& qvec,
    const Eigen::Vector3d& tvec,
    const std::vector<double>& radii,
    const std::vector<double>& focal_lengths,
    const std::vector<double>& theta,
    const Camera& camera) {

    const Eigen::Vector3d proj_point3D =
      QuaternionRotatePoint(qvec, point3D) + tvec; //the 3D point in camera coordinate system

    Eigen::Vector2d proj_point2D;

    if (camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id) {

      const Eigen::Vector2d n = proj_point3D.topRows<2>().normalized(); //the unit direction pointing from camera center to 3D point
      const Eigen::Vector2d point2D_center = camera.ImageToWorld(point2D);
      const double dot_product = n.dot(point2D_center);

      // check that we project onto the correct half-plane
      if (dot_product < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }
      
      proj_point2D = camera.WorldToImage(dot_product * n);
      // check if radii is not empty
      if (radii.size() >= 4) {
        Eigen::Vector2d pp = Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY());
        double radius = (point2D - pp).norm();
        // //////////////// Calculating the theta of the point3D for mapping to the radius and thus focal length ///////////////////////
        double rho = proj_point3D.topRows<2>().norm();
        double z = proj_point3D(2);
        double theta_this = std::atan2(rho, z);
        double focal_length = 0;
        if (camera.ModelId() == ImplicitDistortionModel::model_id) {
          tk::spline<double> s = camera.GetSpline();
          if (!std::isnan(s(radius))) {
            focal_length = s(radius);
          }
        }
        if (int(focal_length) != 0) {
          proj_point2D = proj_point3D.hnormalized() * focal_length + pp;
        }
        else {
          proj_point2D = camera.WorldToImage(dot_product * n);
        }
    
      }
    }
    else {
      // Check that point is infront of camera.
      if (proj_point3D.z() * camera.FocalLength() < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }
      proj_point2D = camera.WorldToImage(proj_point3D.hnormalized());
    }
    // std::cout << "reprojection error: " << (proj_point2D - point2D).squaredNorm() << std::endl;
    return (proj_point2D - point2D).squaredNorm();
  }

  double CalculateSquaredReprojectionError(const Eigen::Vector2d& point2D,
    const Eigen::Vector3d& point3D,
    const Eigen::Matrix3x4d& proj_matrix,
    const Camera& camera) {
    Eigen::Matrix3d R = proj_matrix.leftCols<3>();
    Eigen::Vector3d T = proj_matrix.rightCols<1>();
    const Eigen::Vector3d proj_point3D = R * point3D + T;

    Eigen::Vector2d proj_point2D;

    if (!camera.IsFullyCalibrated(point2D)) {
      const Eigen::Vector2d n = (proj_matrix.topRows<2>() * point3D.homogeneous()).normalized();
      const Eigen::Vector2d point2D_center = camera.ImageToWorld(point2D);
      const double dot_product = n.dot(point2D_center);

      // check that we project onto the correct half-plane
      if (dot_product < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }
      proj_point2D = camera.WorldToImage(dot_product * n);

    }
    else if (camera.ModelId() == ImplicitDistortionModel::model_id) {
      double focal = camera.EvalFocalLength(proj_point3D);

      const double proj_z = proj_matrix.row(2).dot(point3D.homogeneous());

      // Check that point is infront of camera.
      if (proj_z * focal < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }

      const double proj_x = proj_matrix.row(0).dot(point3D.homogeneous());
      const double proj_y = proj_matrix.row(1).dot(point3D.homogeneous());
      const double inv_proj_z = 1.0 / proj_z;

      proj_point2D = camera.WorldToImage(focal * Eigen::Vector2d(inv_proj_z * proj_x, inv_proj_z * proj_y));
    }
    else {
      const double proj_z = proj_matrix.row(2).dot(point3D.homogeneous());

      // Check that point is infront of camera.
      if (proj_z < std::numeric_limits<double>::epsilon()) {
        return std::numeric_limits<double>::max();
      }

      const double proj_x = proj_matrix.row(0).dot(point3D.homogeneous());
      const double proj_y = proj_matrix.row(1).dot(point3D.homogeneous());
      const double inv_proj_z = 1.0 / proj_z;

      proj_point2D = camera.WorldToImage(Eigen::Vector2d(inv_proj_z * proj_x, inv_proj_z * proj_y));

    }

    return (proj_point2D - point2D).squaredNorm();
  }

  double CalculateAngularError(const Eigen::Vector2d& point2D,
    const Eigen::Vector3d& point3D,
    const Eigen::Vector4d& qvec,
    const Eigen::Vector3d& tvec,
    const Camera& camera) {
    Eigen::Matrix3x4d proj_matrix;
    proj_matrix << QuaternionToRotationMatrix(qvec), tvec;
    return CalculateAngularError(point2D, point3D, proj_matrix, camera);
  }

  double CalculateAngularError(const Eigen::Vector2d& point2D,
    const Eigen::Vector3d& point3D,
    const Eigen::Matrix3x4d& proj_matrix,
    const Camera& camera) {
    const Eigen::Vector2d point2D_norm = camera.ImageToWorld(point2D);
    const Eigen::Vector3d point3D_cam = proj_matrix * point3D.homogeneous();
    if (camera.ModelId() == ImplicitDistortionModel::model_id && camera.IsFullyCalibrated(point2D)) {
      // double focal = camera.EvalFocalLength(point2D);
      double focal = camera.EvalFocalLength(point3D_cam);

      const Eigen::Vector2d point2D_norm_full = point2D_norm / focal;
      double angle_error = std::acos(point2D_norm_full.homogeneous().normalized().dot(point3D_cam.normalized()));
      if (focal >= 0) {
        return angle_error;
      }
      else
        return M_PI - angle_error;
    }

    if (camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id) {
      double angle = std::acos(point2D_norm.normalized().dot(point3D_cam.topRows<2>().normalized()));
      if (angle > M_PI / 2.0) {
        return std::acos(point2D_norm.normalized().dot(-point3D_cam.topRows<2>().normalized()));
      }
      else {
        return angle;
      }

    }
    else {
      return std::acos(point2D_norm.homogeneous().normalized().dot(point3D_cam.normalized()));
    }
  }

  double CalculateAngularErrorSplitting(const Eigen::Vector2d& point2D,
    const Eigen::Vector3d& point3D,
    const Eigen::Matrix3x4d& proj_matrix,
    const double& focal_length,
    const Camera& camera) {
    // std::cout << "Using splitting " << std::endl;
    const Eigen::Vector2d point2D_norm = camera.ImageToWorld(point2D) / focal_length;

    const Eigen::Vector3d point3D_cam = proj_matrix * point3D.homogeneous();

    if (camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id) {
      return std::acos(point2D_norm.homogeneous().normalized().dot(point3D_cam.normalized()));
    }
    else {
      return std::acos(point2D_norm.homogeneous().normalized().dot(point3D_cam.normalized()));
    }
  }

  double CalculateDepth(const Eigen::Matrix3x4d& proj_matrix,
    const Eigen::Vector3d& point3D) {
    const double proj_z = proj_matrix.row(2).dot(point3D.homogeneous());
    return proj_z * proj_matrix.col(2).norm();
  }

  bool HasPointPositiveDepth(const Eigen::Matrix3x4d& proj_matrix,
    const Eigen::Vector3d& point3D) {
    return proj_matrix.row(2).dot(point3D.homogeneous()) >=
      std::numeric_limits<double>::epsilon();
  }

}  // namespace colmap
