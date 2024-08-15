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
  } else if (det_K > 0) {
    *K = RR;
  } else {
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
  // if(camera.ModelId()==ImplicitDistortionModel::model_id){
  //     return CalculateSquaredReprojectionErrorFinal(point2D, point3D, qvec, tvec, camera.GetRawRadii(), camera.GetFocalLengthParams(), camera.GetTheta(), camera);
  //   }
  
  if(camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id){
    // if(false) {       
    const Eigen::Vector2d n = proj_point3D.topRows<2>().normalized();
    const Eigen::Vector2d point2D_center = camera.ImageToWorld(point2D);
    const double dot_product = n.dot(point2D_center);

    // check that we project onto the correct half-plane
    // if(fabs(dot_product) < std::numeric_limits<double>::epsilon()) {
      if(dot_product < std::numeric_limits<double>::epsilon()) {
      // if(camera.ModelId() == Radial1DCameraModel::model_id){
        return std::numeric_limits<double>::max();
      // }
    }
    proj_point2D = camera.WorldToImage(dot_product * n);


    // using same triangulation:
    // if (proj_point3D.z() < std::numeric_limits<double>::epsilon()) {
    //   return std::numeric_limits<double>::max();
    // }
    // proj_point2D = camera.WorldToImage(proj_point3D.hnormalized());
    // if(camera.ModelId()==ImplicitDistortionModel::model_id){
    //   return CalculateSquaredReprojectionErrorFinal(point2D, point3D, qvec, tvec, camera.GetRawRadii(), camera.GetFocalLengthParams(), camera.GetTheta(), camera);
    // }

  } else {
    // Check that point is infront of camera.
    if (proj_point3D.z() < std::numeric_limits<double>::epsilon()) {
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
  
  if(camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id){
    
    const Eigen::Vector2d n = proj_point3D.topRows<2>().normalized(); //the unit direction pointing from camera center to 3D point
    const Eigen::Vector2d point2D_center = camera.ImageToWorld(point2D);
    const double dot_product = n.dot(point2D_center);

    // check that we project onto the correct half-plane
    if(dot_product < std::numeric_limits<double>::epsilon()) {
      if(camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id){
      return std::numeric_limits<double>::max();}
    }
    proj_point2D = camera.WorldToImage(dot_product * n);
    // check if radii is not empty
    if(radii.size()>=4) {
      
      // std::cout << "Using full reprojection error" << std::endl;
      Eigen::Vector2d pp = Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY());
      // double radius = (point2D - pp).norm();
      // //////////////// Calculating the theta of the point3D for mapping to the radius and thus focal length ///////////////////////
      double rho = proj_point3D.topRows<2>().norm();
      double z = proj_point3D(2);
      double theta_this = std::atan2(rho,z);
      double radius;
      // std::cout << "theta size: " << theta.size() << std::endl;
      // std::cout << "radii size: " << radii.size() << std::endl;

      if (theta_this <= theta[0]) {
            // special case, extrpolate from first points
            double dtheta = theta[1] - theta[0];
            double alpha = (theta_this - theta[0]) / dtheta;
            radius = (1.0 - alpha) * radii[0] + alpha * radii[1];
            
        } else if (theta_this >= theta[theta.size()-1]) {
            // special case, extrapolate from last points
            double dtheta = theta[theta.size()-1] - theta[theta.size()-2];
            double alpha = (theta_this - theta[theta.size()-2]) / dtheta;
            radius = (1.0 - alpha) * radii[theta.size()-2] + alpha * radii[theta.size()-1];
        
        } else {
            // we interpolate between the two neighboring calibration points
            // concatenate theta and radii into a vector of pairs: theta_r
            std::vector<std::pair<double,double>> theta_r;
            for (size_t i = 0; i < radii.size(); i++) {
                theta_r.push_back(std::make_pair(theta[i], radii[i]));
            }
            auto it = std::lower_bound(theta_r.begin(), theta_r.end(), std::make_pair(theta_this, std::numeric_limits<double>::lowest()));
            
            const std::pair<double,double> &theta_r1 = *it;
            const std::pair<double,double> &theta_r0 = *(--it);
            
            // Interpolate radius           
            double dtheta = theta_r1.first - theta_r0.first;
            double alpha = (theta_this - theta_r0.first) / dtheta;
            radius = (1.0 - alpha) * theta_r0.second + alpha * theta_r1.second;
        }
      // std::cout << "calculated radius: " << radius << std::endl;
      // std::cout << "plain radius: " << (point2D - pp).norm() << std::endl;

      // get the focal length for this radius
      radius = (point2D - pp).norm();
      double focal_length = 0;
      for (int i = 0; i < radii.size() - 1; i++) {
        // std::cout << "radii " << i << " "<<radii[i] << std::endl;
        if (radius >= radii[i] && radius <= radii[i+1]) {
          // interpolate the focal length

          focal_length = focal_lengths[i] + (focal_lengths[i+1] - focal_lengths[i]) * (radius - radii[i]) / (radii[i+1] - radii[i] + std::numeric_limits<double>::epsilon());
          break;
        }
      }
      // if the radius is smaller than the smallest radius, set the focal length to the smallest focal length
      if (radius < radii[0]) {
        focal_length = focal_lengths[0];
      }
      // if the radius is larger than the largest radius, set the focal length to the largest focal length
      if (radius > radii[radii.size() - 1]) {
        focal_length = focal_lengths[radii.size() - 1];
      }
      // if (camera.ModelId()==ImplicitDistortionModel::model_id){
      //   tk::spline s = camera.GetSpline();
      //   focal_length = s(radius);
      // }
      // std::cout << "focal length: " << focal_length << std::endl;
      // proj_point2D(0) = proj_point3D.hnormalized()(0) * focal_length + pp(0);
      // proj_point2D(1) = proj_point3D.hnormalized()(1) * focal_length + pp(1);
      if (focal_length > 100){
      proj_point2D= proj_point3D.hnormalized() * focal_length + pp;
      // std::cout<<"proj_point2D: "<<proj_point2D<<std::endl;
      }else{
        proj_point2D = camera.WorldToImage(proj_point3D.hnormalized());
      }
      // std::cout << "pp: " << pp << std::endl;
      // std::cout << "proj_point2D: " << proj_point2D << std::endl;
      // std::cout << "point2D: " << point2D << std::endl;
      // std::cout << "calculated reprojection error: " << (proj_point2D - point2D).squaredNorm() << std::endl;
      // std::cout << "used radius: " << radius << std::endl;
      // std::cout << "used focal length: " << focal_length << std::endl;
      // double actual_radius = (proj_point2D - pp).norm();
      // double actual_f ;
      // for (int i = 0; i < radii.size() - 1; i++) {
      //   // std::cout << "radii " << i << " "<<radii[i] << std::endl;
      //   if (actual_radius >= radii[i] && actual_radius <= radii[i+1]) {
      //     // interpolate the focal length

      //     actual_f = focal_lengths[i] + (focal_lengths[i+1] - focal_lengths[i]) * (radius - radii[i]) / (radii[i+1] - radii[i]+ std::numeric_limits<double>::epsilon());
      //     break;
      //   }
      // }

      // std::cout << "actual radius: " << (proj_point2D - pp).norm() << std::endl;
      // std::cout << "actual focal length: " << actual_f << std::endl;
      // check that point is infront of camera.
      // if (proj_point3D.z() < std::numeric_limits<double>::epsilon()) {
      //   return std::numeric_limits<double>::max();
      // }
    }
  } else {
    // Check that point is infront of camera.
    if (proj_point3D.z() < std::numeric_limits<double>::epsilon()) {
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

  Eigen::Vector2d proj_point2D;

  if(camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id){        
    const Eigen::Vector2d n = (proj_matrix.topRows<2>() * point3D.homogeneous()).normalized();
    const Eigen::Vector2d point2D_center = camera.ImageToWorld(point2D);
    const double dot_product = n.dot(point2D_center);

    // check that we project onto the correct half-plane
    if(dot_product < std::numeric_limits<double>::epsilon()) {
    return std::numeric_limits<double>::max();
    }
    proj_point2D = camera.WorldToImage(dot_product * n);

  } else {
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

  if(camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id){    
    return std::acos(point2D_norm.normalized().dot(point3D_cam.topRows<2>().normalized()));
    // return std::acos(point2D_norm.homogeneous().normalized().dot(point3D_cam.normalized()));
  } else {
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

  if(camera.ModelId() == Radial1DCameraModel::model_id || camera.ModelId() == ImplicitDistortionModel::model_id){ 

    // return 0.5*(std::acos(point2D_norm.homogeneous().normalized().dot(point3D_cam.normalized()))) + 0.5*(std::acos((point2D_norm * focal_length).normalized().dot(point3D_cam.topRows<2>().normalized())));
    return std::acos(point2D_norm.homogeneous().normalized().dot(point3D_cam.normalized())); 
  } else {
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
