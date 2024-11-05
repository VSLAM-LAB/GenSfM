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

#ifndef COLMAP_SRC_BASE_IMAGE_H_
#define COLMAP_SRC_BASE_IMAGE_H_

#include <string>
#include <vector>
#include<optional>
#include <algorithm>
#include <random>

#include <Eigen/Core>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp> 
#include <boost/math/interpolators/cubic_hermite.hpp>
#include <unsupported/Eigen/Splines>
#include "base/camera.h"
#include "base/point2d.h"
#include "base/visibility_pyramid.h"
#include "base/spline.h"
#include "util/alignment.h"
#include "util/logging.h"
#include "util/math.h"
#include "util/types.h"

namespace colmap {

// Class that holds information about an image. An image is the product of one
// camera shot at a certain location (parameterized as the pose). An image may
// share a camera with multiple other images, if its intrinsics are the same.
class Image {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Image();
  typedef Eigen::Spline<double, 1> Spline1D;

  // Setup / tear down the image and necessary internal data structures before
  // and after being used in reconstruction.
  void SetUp(const Camera& camera);
  void TearDown();
  mutable bool use_radial_ = false;
  inline bool UseRadial() const;
  inline void SetUseRadial(const bool use_radial_);

  // Access the unique identifier of the image.
  inline image_t ImageId() const;
  inline void SetImageId(const image_t image_id);

  // Access the name of the image.
  inline const std::string& Name() const;
  inline std::string& Name();
  inline void SetName(const std::string& name);

  // Access the unique identifier of the camera. Note that multiple images
  // might share the same camera.
  inline camera_t CameraId() const;
  inline void SetCameraId(const camera_t camera_id);
  // Check whether identifier of camera has been set.
  inline bool HasCamera() const;

  // Check if image is registered.
  inline bool IsRegistered() const;
  inline void SetRegistered(const bool registered);

  // Get the number of image points.
  inline point2D_t NumPoints2D() const;

  // Get the number of triangulations, i.e. the number of points that
  // are part of a 3D point track.
  inline point2D_t NumPoints3D() const;

  // Get the number of observations, i.e. the number of image points that
  // have at least one correspondence to another image.
  inline point2D_t NumObservations() const;
  inline void SetNumObservations(const point2D_t num_observations);

  // Get the number of correspondences for all image points.
  inline point2D_t NumCorrespondences() const;
  inline void SetNumCorrespondences(const point2D_t num_observations);

  // Get the number of observations that see a triangulated point, i.e. the
  // number of image points that have at least one correspondence to a
  // triangulated point in another image.
  inline point2D_t NumVisiblePoints3D() const;

  // Get the score of triangulated observations. In contrast to
  // `NumVisiblePoints3D`, this score also captures the distribution
  // of triangulated observations in the image. This is useful to select
  // the next best image in incremental reconstruction, because a more
  // uniform distribution of observations results in more robust registration.
  inline size_t Point3DVisibilityScore() const;

  // Parameters for retrieving point-wise focal length from implicit_distortion.
  inline std::vector<double> GetFocalLengthParams() const;
  inline void SetFocalLengthParams(const std::vector<double>& focal_length_params) const;
  inline std::vector<double> GetRawRadii() const; 
  inline void SetRawRadii(const std::vector<double>& raw_radii) const;
  inline std::vector<double> GetTheta() const;
  inline void SetTheta(const std::vector<double>& theta) const;
  // inline boost::math::interpolators::cubic_hermite<std::vector<double>> GetSpline() const;
  // inline void SetSpline(const boost::math::interpolators::cubic_hermite<std::vector<double>>& spline) const;

  // Access quaternion vector as (qw, qx, qy, qz) specifying the rotation of the
  // pose which is defined as the transformation from world to image space.
  inline const Eigen::Vector4d& Qvec() const;
  inline Eigen::Vector4d& Qvec();
  inline double Qvec(const size_t idx) const;
  inline double& Qvec(const size_t idx);
  inline void SetQvec(const Eigen::Vector4d& qvec);

  // Quaternion prior, e.g. given by EXIF gyroscope tag.
  inline const Eigen::Vector4d& QvecPrior() const;
  inline Eigen::Vector4d& QvecPrior();
  inline double QvecPrior(const size_t idx) const;
  inline double& QvecPrior(const size_t idx);
  inline bool HasQvecPrior() const;
  inline void SetQvecPrior(const Eigen::Vector4d& qvec);

  // Access translation vector as (tx, ty, tz) specifying the translation of the
  // pose which is defined as the transformation from world to image space.
  inline const Eigen::Vector3d& Tvec() const;
  inline Eigen::Vector3d& Tvec();
  inline double Tvec(const size_t idx) const;
  inline double& Tvec(const size_t idx);
  inline void SetTvec(const Eigen::Vector3d& tvec);

  // Translation prior, e.g. given by EXIF GPS tag.
  inline const Eigen::Vector3d& TvecPrior() const;
  inline Eigen::Vector3d& TvecPrior();
  inline double TvecPrior(const size_t idx) const;
  inline double& TvecPrior(const size_t idx);
  inline bool HasTvecPrior() const;
  inline void SetTvecPrior(const Eigen::Vector3d& tvec);

  // Access the coordinates of image points.
  inline const class Point2D& Point2D(const point2D_t point2D_idx) const;
  inline class Point2D& Point2D(const point2D_t point2D_idx);
  inline const std::vector<class Point2D>& Points2D() const;
  void SetPoints2D(const std::vector<Eigen::Vector2d>& points);
  void SetPoints2D(const std::vector<class Point2D>& points);

  // Set the point as triangulated, i.e. it is part of a 3D point track.
  void SetPoint3DForPoint2D(const point2D_t point2D_idx,
                            const point3D_t point3D_id);

  // Set the point as not triangulated, i.e. it is not part of a 3D point track.
  void ResetPoint3DForPoint2D(const point2D_t point2D_idx);

  // Check whether an image point has a correspondence to an image point in
  // another image that has a 3D point.
  inline bool IsPoint3DVisible(const point2D_t point2D_idx) const;

  // Check whether one of the image points is part of the 3D point track.
  bool HasPoint3D(const point3D_t point3D_id) const;

  // Indicate that another image has a point that is triangulated and has
  // a correspondence to this image point. Note that this must only be called
  // after calling `SetUp`.
  void IncrementCorrespondenceHasPoint3D(const point2D_t point2D_idx);

  // Fitting the spline model for the focal length - radii mapping
  // inline void FitSpline(std::vector<double>& radii,  std::vector<double>& focal_lengths) const;

  // inline double EvalFocalLength(double radius) const;

  // Indicate that another image has a point that is not triangulated any more
  // and has a correspondence to this image point. This assumes that
  // `IncrementCorrespondenceHasPoint3D` was called for the same image point
  // and correspondence before. Note that this must only be called
  // after calling `SetUp`.
  void DecrementCorrespondenceHasPoint3D(const point2D_t point2D_idx);

  // Normalize the quaternion vector.
  void NormalizeQvec();

  // Compose the projection matrix from world to image space.
  Eigen::Matrix3x4d ProjectionMatrix() const;

  // Compose the inverse projection matrix from image to world space
  Eigen::Matrix3x4d InverseProjectionMatrix() const;

  // Compose rotation matrix from quaternion vector.
  Eigen::Matrix3d RotationMatrix() const;

  // Extract the projection center in world space.
  Eigen::Vector3d ProjectionCenter() const;

  // Extract the viewing direction of the image.
  Eigen::Vector3d ViewingDirection() const;

  // The number of levels in the 3D point multi-resolution visibility pyramid.
  static const int kNumPoint3DVisibilityPyramidLevels;

 private:
  // Identifier of the image, if not specified `kInvalidImageId`.
  image_t image_id_;

  // The name of the image, i.e. the relative path.
  std::string name_;

  // The identifier of the associated camera. Note that multiple images might
  // share the same camera. If not specified `kInvalidCameraId`.
  camera_t camera_id_;

  // Whether the image is successfully registered in the reconstruction.
  bool registered_;

  // Return the spline value for a new radius
  

  // The number of 3D points the image observes, i.e. the sum of its `points2D`
  // where `point3D_id != kInvalidPoint3DId`.
  point2D_t num_points3D_;

  // The number of image points that have at least one correspondence to
  // another image.
  point2D_t num_observations_;

  // The sum of correspondences per image point.
  point2D_t num_correspondences_;

  // The number of 2D points, which have at least one corresponding 2D point in
  // another image that is part of a 3D point track, i.e. the sum of `points2D`
  // where `num_tris > 0`.
  point2D_t num_visible_points3D_;

  // The pose of the image, defined as the transformation from world to image.
  Eigen::Vector4d qvec_;
  Eigen::Vector3d tvec_;

  // The pose prior of the image, e.g. extracted from EXIF tags.
  Eigen::Vector4d qvec_prior_;
  Eigen::Vector3d tvec_prior_;

  // Parameters for retrieving point-wise focal length from implicit_distortion.
  mutable std::vector<double> focal_length_params_;
  mutable std::vector<double> raw_radii_;
  mutable std::vector<double> theta_;
  // mutable boost::math::interpolators::cubic_hermite<std::vector<double>> spline_;

  // All image points, including points that are not part of a 3D point track.
  std::vector<class Point2D> points2D_;
  
  // Per image point, the number of correspondences that have a 3D point.
  std::vector<image_t> num_correspondences_have_point3D_;

  // Data structure to compute the distribution of triangulated correspondences
  // in the image. Note that this structure is only usable after `SetUp`.
  VisibilityPyramid point3D_visibility_pyramid_;

  // The 1D spline for the focal length
  // mutable Spline1D spline_;
  // mutable tk::spline spline_;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

image_t Image::ImageId() const { return image_id_; }

void Image::SetImageId(const image_t image_id) { image_id_ = image_id; }

const std::string& Image::Name() const { return name_; }

std::string& Image::Name() { return name_; }

void Image::SetName(const std::string& name) { name_ = name; }

inline camera_t Image::CameraId() const { return camera_id_; }

inline void Image::SetCameraId(const camera_t camera_id) {
  CHECK_NE(camera_id, kInvalidCameraId);
  camera_id_ = camera_id;
}

inline bool Image::HasCamera() const { return camera_id_ != kInvalidCameraId; }

bool Image::IsRegistered() const { return registered_; }

void Image::SetRegistered(const bool registered) { registered_ = registered; }

point2D_t Image::NumPoints2D() const {
  return static_cast<point2D_t>(points2D_.size());
}
inline bool Image::UseRadial() const { return use_radial_; }
inline void Image::SetUseRadial(const bool use_radial) { use_radial_ = use_radial; }
point2D_t Image::NumPoints3D() const { return num_points3D_; }

point2D_t Image::NumObservations() const { return num_observations_; }

void Image::SetNumObservations(const point2D_t num_observations) {
  num_observations_ = num_observations;
}

point2D_t Image::NumCorrespondences() const { return num_correspondences_; }

void Image::SetNumCorrespondences(const point2D_t num_correspondences) {
  num_correspondences_ = num_correspondences;
}

point2D_t Image::NumVisiblePoints3D() const { return num_visible_points3D_; }

size_t Image::Point3DVisibilityScore() const {
  return point3D_visibility_pyramid_.Score();
}

const Eigen::Vector4d& Image::Qvec() const { return qvec_; }

Eigen::Vector4d& Image::Qvec() { return qvec_; }

inline double Image::Qvec(const size_t idx) const { return qvec_(idx); }

inline double& Image::Qvec(const size_t idx) { return qvec_(idx); }

inline std::vector<double> Image::GetFocalLengthParams() const {
  return focal_length_params_;
}

inline void Image::SetFocalLengthParams(const std::vector<double>& focal_length_params) const {
  focal_length_params_ = focal_length_params;
}

// inline void Image::FitSpline(std::vector<double>& radii, std::vector<double>& focal_lengths) const{
//   // Convert std::vector to Eigen vectors
//   assert(radii.size() == focal_lengths.size());

//   std::vector<int> increasing_indices;
//   for (int i = 0; i < radii.size() - 1; i++) {
//     if(radii[i+1] > radii[i]) {
//       increasing_indices.push_back(i);
//     }
//   }
//   std::vector<double> new_radii;
//   std::vector<double> new_focal_lengths;
//   if (!increasing_indices.empty() && increasing_indices.back() != radii.size() - 1) {
//         increasing_indices.push_back(radii.size() - 1);
//     }

//   for (int i = 0; i < increasing_indices.size(); i++) {
//     new_radii.push_back(radii[increasing_indices[i]]);
//     new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
//   }
//   std::cout << "new_radii size: " << new_radii.size() << std::endl;
//   std::cout << "new_focal_lengths size: " << new_focal_lengths.size() << std::endl;

//   // use ransac technique to fit the spline
//   int best_inliers_count = 0;
//   std::vector<double> best_coeffs;
//   std::random_device rd;
//   std::mt19937 gen(rd());
//   std::uniform_int_distribution<> dis(0, new_radii.size() - 1);
//   tk::spline best_spline;
//   const int max_iterations = 20;
//   const double threshold = 5.0;
//   int degree = 9;
//   for (int i = 0; i < max_iterations; ++i){
//     std::vector<std::pair<double, double>> samples;
//     for (int j = 0; j < degree; ++j) {
//       int idx = dis(gen);
//       samples.emplace_back(new_radii[idx], new_focal_lengths[idx]);
//     }
//     std::sort(samples.begin(), samples.end());
//     std::vector<double> sample_x, sample_y;
//         for (const auto& pair : samples) {
//             sample_x.push_back(pair.first);
//             sample_y.push_back(pair.second);
//         }
//     tk::spline s;
//     s.set_points(sample_x, sample_y);
//     int inliers_count = 0;
//     for (size_t j = 0; j < new_radii.size(); ++j) {
//       // std::cout << "new_radii[j]: " << new_radii[j] << std::endl;
//       // check if s is empty
//       // check the range of sample_x
//       auto min_max_x = std::minmax_element(sample_x.begin(), sample_x.end());
//       auto x = s.get_x();
//       auto min_max_x_s = std::minmax_element(x.begin(), x.end());
//       std::cout << "Min x value: " << *min_max_x.first << ", Max x value: " << *min_max_x.second << std::endl;
//       std::cout << "Min x value for s: " << *min_max_x_s.first << ", Max x value for s: " << *min_max_x_s.second << std::endl;
//         // if (s.get_x().empty()) {
//         //     continue;
//         // }
//         double y_est;
//         try{
//          y_est = s(new_radii[j]);
//         }catch(const std::exception& e){
//           std::cout << "Exception caught: " << e.what() << std::endl;
//           // y_est = new_focal_lengths[j];
//           continue;
//         }
        
//         // double y_est = s(new_radii[j]);
//         std::cout << "j: " << j << "new_radii[j]" << new_radii[j]<< std::endl;
//         std::cout << "y_est for s: " << y_est << std::endl;
//         if (fabs(y_est - new_focal_lengths[j]) < threshold) {
//             ++inliers_count;
//         }
//     }
//     if (inliers_count >= best_inliers_count) {
//             best_inliers_count = inliers_count;
//             best_spline = s;
//         }
//   }

//   std::vector<double> inliers_x, inliers_y;
//     for (size_t j = 0; j < new_radii.size(); ++j) {
      
//         double y_est = best_spline(new_radii[j]);
//         std::cout << "y_est for best_spline: " << y_est << std::endl;
//         if (fabs(y_est - new_focal_lengths[j]) < threshold) {
//             inliers_x.push_back(new_radii[j]);
//             inliers_y.push_back(new_focal_lengths[j]);
//         }
//     }
//   // if(inliers_x.size() >= 4){
//   //  spline_.set_points(inliers_x, inliers_y);
//   // }else{
//   //   spline_.set_points(new_radii, new_focal_lengths);
//   // }
//   spline_ = best_spline;
// }

// double Image::EvalFocalLength(double radius) const {
//   return spline_(radius);
  
// }

inline std::vector<double> Image::GetRawRadii() const {
  return raw_radii_;
}

inline void Image::SetRawRadii(const std::vector<double>& raw_radii) const {
  raw_radii_ = raw_radii;
}

inline std::vector<double> Image::GetTheta() const {
  return theta_;
}

inline void Image::SetTheta(const std::vector<double>& theta) const {
  theta_ = theta;
}

// inline boost::math::interpolators::cubic_hermite<std::vector<double>> Image::GetSpline() const {
//   return spline_;
// }

// inline void Image::SetSpline(const boost::math::interpolators::cubic_hermite<std::vector<double>>& spline) const {
//   spline_ = spline;
// }
void Image::SetQvec(const Eigen::Vector4d& qvec) { qvec_ = qvec; }

const Eigen::Vector4d& Image::QvecPrior() const { return qvec_prior_; }

Eigen::Vector4d& Image::QvecPrior() { return qvec_prior_; }

inline double Image::QvecPrior(const size_t idx) const {
  return qvec_prior_(idx);
}

inline double& Image::QvecPrior(const size_t idx) { return qvec_prior_(idx); }

inline bool Image::HasQvecPrior() const { return !IsNaN(qvec_prior_.sum()); }

void Image::SetQvecPrior(const Eigen::Vector4d& qvec) { qvec_prior_ = qvec; }

const Eigen::Vector3d& Image::Tvec() const { return tvec_; }

Eigen::Vector3d& Image::Tvec() { return tvec_; }

inline double Image::Tvec(const size_t idx) const { return tvec_(idx); }

inline double& Image::Tvec(const size_t idx) { return tvec_(idx); }

void Image::SetTvec(const Eigen::Vector3d& tvec) { tvec_ = tvec; }

const Eigen::Vector3d& Image::TvecPrior() const { return tvec_prior_; }

Eigen::Vector3d& Image::TvecPrior() { return tvec_prior_; }

inline double Image::TvecPrior(const size_t idx) const {
  return tvec_prior_(idx);
}

inline double& Image::TvecPrior(const size_t idx) { return tvec_prior_(idx); }

inline bool Image::HasTvecPrior() const { return !IsNaN(tvec_prior_.sum()); }

void Image::SetTvecPrior(const Eigen::Vector3d& tvec) { tvec_prior_ = tvec; }

const class Point2D& Image::Point2D(const point2D_t point2D_idx) const {
  return points2D_.at(point2D_idx);
}

class Point2D& Image::Point2D(const point2D_t point2D_idx) {
  return points2D_.at(point2D_idx);
}

const std::vector<class Point2D>& Image::Points2D() const { return points2D_; }

bool Image::IsPoint3DVisible(const point2D_t point2D_idx) const {
  return num_correspondences_have_point3D_.at(point2D_idx) > 0;
}

}  // namespace colmap

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION_CUSTOM(colmap::Image)

#endif  // COLMAP_SRC_BASE_IMAGE_H_
