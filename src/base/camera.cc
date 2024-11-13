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

#include "base/camera.h"

#include <iomanip>

#include "base/camera_models.h"
#include "util/logging.h"
#include "util/misc.h"

namespace colmap {

Camera::Camera()
    : camera_id_(kInvalidCameraId),
      model_id_(kInvalidCameraModelId),
      width_(0),
      height_(0),
      prior_focal_length_(false) {}

std::string Camera::ModelName() const { return CameraModelIdToName(model_id_); }

void Camera::SetModelId(const int model_id) {
  CHECK(ExistsCameraModelWithId(model_id));
  model_id_ = model_id;
  params_.resize(CameraModelNumParams(model_id_), 0);
}

void Camera::SetModelIdFromName(const std::string& model_name) {
  CHECK(ExistsCameraModelWithName(model_name));
  model_id_ = CameraModelNameToId(model_name);
  params_.resize(CameraModelNumParams(model_id_), 0);
}

const std::vector<size_t>& Camera::FocalLengthIdxs() const {
  return CameraModelFocalLengthIdxs(model_id_);
}

const std::vector<size_t>& Camera::PrincipalPointIdxs() const {
  return CameraModelPrincipalPointIdxs(model_id_);
}

const std::vector<size_t>& Camera::ExtraParamsIdxs() const {
  return CameraModelExtraParamsIdxs(model_id_);
}

const std::vector<size_t>& Camera::XParamsIdx() const {
  return CameraModelXParamsIdxs(model_id_);
}

Eigen::Matrix3d Camera::CalibrationMatrix() const {
  Eigen::Matrix3d K = Eigen::Matrix3d::Identity();

  const std::vector<size_t>& idxs = FocalLengthIdxs();
  if (idxs.size() == 1) {
    K(0, 0) = params_[idxs[0]];
    K(1, 1) = params_[idxs[0]];
  } else if (idxs.size() == 2) {
    K(0, 0) = params_[idxs[0]];
    K(1, 1) = params_[idxs[1]];
  } else {
    LOG(FATAL)
        << "Camera model must either have 1 or 2 focal length parameters.";
  }

  K(0, 2) = PrincipalPointX();
  K(1, 2) = PrincipalPointY();

  return K;
}

std::string Camera::ParamsInfo() const {
  return CameraModelParamsInfo(model_id_);
}

double Camera::MeanFocalLength() const {
  const auto& focal_length_idxs = FocalLengthIdxs();
  double focal_length = 0;
  for (const auto idx : focal_length_idxs) {
    focal_length += params_[idx];
  }
  return focal_length / focal_length_idxs.size();
}

double Camera::FocalLength() const {
  const std::vector<size_t>& idxs = FocalLengthIdxs();
  CHECK_EQ(idxs.size(), 1);
  return params_[idxs[0]];
}

double Camera::FocalLengthX() const {
  const std::vector<size_t>& idxs = FocalLengthIdxs();
  CHECK_EQ(idxs.size(), 2);
  return params_[idxs[0]];
}

double Camera::FocalLengthY() const {
  const std::vector<size_t>& idxs = FocalLengthIdxs();
  CHECK_EQ(idxs.size(), 2);
  return params_[idxs[1]];
}

void Camera::SetFocalLength(const double focal_length) {
  const std::vector<size_t>& idxs = FocalLengthIdxs();
  for (const auto idx : idxs) {
    params_[idx] = focal_length;
  }
}

void Camera::SetFocalLengthX(const double focal_length_x) {
  const std::vector<size_t>& idxs = FocalLengthIdxs();
  CHECK_EQ(idxs.size(), 2);
  params_[idxs[0]] = focal_length_x;
}

void Camera::SetFocalLengthY(const double focal_length_y) {
  const std::vector<size_t>& idxs = FocalLengthIdxs();
  CHECK_EQ(idxs.size(), 2);
  params_[idxs[1]] = focal_length_y;
}

double Camera::PrincipalPointX() const {
  const std::vector<size_t>& idxs = PrincipalPointIdxs();
  CHECK_EQ(idxs.size(), 2);
  return params_[idxs[0]];
}

double Camera::PrincipalPointY() const {
  const std::vector<size_t>& idxs = PrincipalPointIdxs();
  CHECK_EQ(idxs.size(), 2);
  return params_[idxs[1]];
}

void Camera::SetPrincipalPointX(const double ppx) {
  const std::vector<size_t>& idxs = PrincipalPointIdxs();
  CHECK_EQ(idxs.size(), 2);
  params_[idxs[0]] = ppx;
}

void Camera::SetPrincipalPointY(const double ppy) {
  const std::vector<size_t>& idxs = PrincipalPointIdxs();
  CHECK_EQ(idxs.size(), 2);
  params_[idxs[1]] = ppy;
}

std::string Camera::ParamsToString() const { return VectorToCSV(params_); }

bool Camera::SetParamsFromString(const std::string& string) {
  const std::vector<double> new_camera_params = CSVToVector<double>(string);
  if (!CameraModelVerifyParams(model_id_, new_camera_params)) {
    return false;
  }

  params_ = new_camera_params;
  return true;
}

bool Camera::VerifyParams() const {
  return CameraModelVerifyParams(model_id_, params_);
}

bool Camera::HasBogusParams(const double min_focal_length_ratio,
                            const double max_focal_length_ratio,
                            const double max_extra_param) const {
  return CameraModelHasBogusParams(model_id_, params_, width_, height_,
                                   min_focal_length_ratio,
                                   max_focal_length_ratio, max_extra_param);
}

void Camera::InitializeWithId(const int model_id, const double focal_length,
                              const size_t width, const size_t height) {
  CHECK(ExistsCameraModelWithId(model_id));
  model_id_ = model_id;
  width_ = width;
  height_ = height;
  params_ = CameraModelInitializeParams(model_id, focal_length, width, height);
}

void Camera::InitializeWithName(const std::string& model_name,
                                const double focal_length, const size_t width,
                                const size_t height) {
  InitializeWithId(CameraModelNameToId(model_name), focal_length, width,
                   height);
}

Eigen::Vector2d Camera::ImageToWorld(const Eigen::Vector2d& image_point) const {
  Eigen::Vector2d world_point;
  CameraModelImageToWorld(model_id_, params_, image_point(0), image_point(1),
                          &world_point(0), &world_point(1));
  return world_point;
}

double Camera::ImageToWorldThreshold(const double threshold) const {
  return CameraModelImageToWorldThreshold(model_id_, params_, threshold);
}

Eigen::Vector2d Camera::WorldToImage(const Eigen::Vector2d& world_point) const {
  Eigen::Vector2d image_point;
  CameraModelWorldToImage(model_id_, params_, world_point(0), world_point(1),
                          &image_point(0), &image_point(1));
  return image_point;
}

void Camera::Rescale(const double scale) {
  CHECK_GT(scale, 0.0);
  const double scale_x =
      std::round(scale * width_) / static_cast<double>(width_);
  const double scale_y =
      std::round(scale * height_) / static_cast<double>(height_);
  width_ = static_cast<size_t>(std::round(scale * width_));
  height_ = static_cast<size_t>(std::round(scale * height_));
  SetPrincipalPointX(scale_x * PrincipalPointX());
  SetPrincipalPointY(scale_y * PrincipalPointY());
  if (FocalLengthIdxs().size() == 1) {
    SetFocalLength((scale_x + scale_y) / 2.0 * FocalLength());
  } else if (FocalLengthIdxs().size() == 2) {
    SetFocalLengthX(scale_x * FocalLengthX());
    SetFocalLengthY(scale_y * FocalLengthY());
  } else {
    LOG(FATAL)
        << "Camera model must either have 1 or 2 focal length parameters.";
  }
}

void Camera::Rescale(const size_t width, const size_t height) {
  const double scale_x =
      static_cast<double>(width) / static_cast<double>(width_);
  const double scale_y =
      static_cast<double>(height) / static_cast<double>(height_);
  width_ = width;
  height_ = height;
  SetPrincipalPointX(scale_x * PrincipalPointX());
  SetPrincipalPointY(scale_y * PrincipalPointY());
  if (FocalLengthIdxs().size() == 1) {
    SetFocalLength((scale_x + scale_y) / 2.0 * FocalLength());
  } else if (FocalLengthIdxs().size() == 2) {
    SetFocalLengthX(scale_x * FocalLengthX());
    SetFocalLengthY(scale_y * FocalLengthY());
  } else {
    LOG(FATAL)
        << "Camera model must either have 1 or 2 focal length parameters.";
  }
}

bool Camera::FitPIeceWiseSpline_binary(std::vector<double>& radii, std::vector<double>& focal_lengths, std::vector<double>& principal_point) {
  assert(radii.size() == focal_lengths.size());
  // Ensure that the radii are sorted in ascending order
  std::vector<int> increasing_indices;
  std::vector<double> new_radii;
  std::vector<double> new_focal_lengths;
  std::vector<std::vector<double>> uncalibrated_areas = {};
  for (int i = 0; i < radii.size() - 1; i++) {
    if(radii[i+1] > radii[i]) {
      increasing_indices.push_back(i);
    }
  }
  if (!increasing_indices.empty() && increasing_indices.back() != radii.size() - 1) {
        increasing_indices.push_back(radii.size() - 1);
    }
  for(int i = 0; i < increasing_indices.size(); i++){
    new_radii.push_back(radii[increasing_indices[i]]);
    new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
  } 
  std::vector<double> intervals = {};
  for (int i = 0; i < new_radii.size()-1; i++) {
    intervals.push_back(new_radii[i+1] - new_radii[i]);
  }
  double mean_interval = std::accumulate(intervals.begin(), intervals.end(), 0.0)/intervals.size();
  double std_interval = 0.0;
  for (int i = 0; i < intervals.size(); i++) {
    std_interval += pow(intervals[i] - mean_interval, 2);
  }
  std_interval = sqrt(std_interval/intervals.size());
  std::vector<std::vector<double>> radii_segments = {};
  std::vector<std::vector<double>> focal_lengths_segments = {};
  double threshold = mean_interval + std_interval;
  double std_threshold = 0.5*std_interval;
  std::cout << "!!! Original threshold: " << threshold;
  threshold = std::max(threshold, DegToRad(0.1));
  std::cout << ", New threshold: " << threshold << std::endl;
  recursiveSplit(new_radii, new_focal_lengths, radii_segments, focal_lengths_segments, threshold, std_threshold);

  if (radii_segments.size() == 0) {
    return false;
  }
  // print the beginning and end of the longest segment
  // find the longest segment
  int longest_segment = 0;
  int longest_segment_size = 0;
  for(int i = 0; i < radii_segments.size(); i++){
    if(radii_segments[i].size() > longest_segment_size){
      longest_segment = i;
      longest_segment_size = radii_segments[i].size();
    }
  }
  
  std::vector<double> calibrated_range = {focal_lengths_segments[longest_segment].front(), focal_lengths_segments[longest_segment].back()};
  std::vector<double> radii_calibrated = radii_segments[longest_segment];
  std::vector<double> focal_lengths_calibrated = focal_lengths_segments[longest_segment];
  int max_it = 80;
  int degree = (params_.size() - 2) / 2;
  double threshold_ransac = 5.0;
  tk::spline<double> best_spline = ransac_spline(max_it, degree, threshold_ransac, radii_calibrated, focal_lengths_calibrated);
  std::vector<double> used_x;
  std::vector<double> used_y;

  // std::cout << "inliers_count: " << best_spline.get_inliers_count() << std::endl;
  int inliers_count = 0;
  for (size_t j = 0; j < radii_calibrated.size(); ++j) {
    double y_est = best_spline(radii_calibrated[j]);
    if (fabs(y_est - focal_lengths_calibrated[j]) < threshold_ransac) {
        ++inliers_count;
    }
  }

  std::ofstream out("radii_vs_focal.txt");
  for (size_t j = 0; j < radii_calibrated.size(); ++j) {
    out << radii_calibrated[j] << " " << focal_lengths_calibrated[j] << std::endl;
  }
  out.close();

  std::cout << "inliers_count / total count " << inliers_count << " / " << radii_calibrated.size() << std::endl;

  // If only a few points are inliers, then the calibration is not good
  if (inliers_count < 0.9 * radii_calibrated.size()) {
    return false;
  }

  // If only a few points are used to fit the spline, then the calibration is not good
  if (radii_calibrated.size() < 500) {
    return false;
  }

  if (ModelId()==ImplicitDistortionModel::model_id) {
    used_x = best_spline.get_x();
    used_y = best_spline.get_y();
    std::vector<double> updated_params;
    updated_params.push_back(principal_point[0]);
    updated_params.push_back(principal_point[1]);
    for (int i = 0; i < degree; i++) {
      updated_params.push_back(used_x[i]);
    }
    for(int i = 0; i < degree; i++) {
      updated_params.push_back(used_y[i]);
    }
    updated_params_ = updated_params;
    params_ = updated_params_;
    SetParams(updated_params);
  }
  // Directly use theta-r mapping
  spline_ = best_spline;
  return true;
}
}  // namespace colmap
