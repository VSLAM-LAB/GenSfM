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

#ifndef COLMAP_SRC_BASE_CAMERA_H_
#define COLMAP_SRC_BASE_CAMERA_H_

#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <Eigen/Core>
#include "base/spline.h"
#include "base/camera_models.h"

#include "util/math.h"
#include "util/types.h"
#include "estimators/spline_fitting.h"

namespace colmap {

// Camera class that holds the intrinsic parameters. Cameras may be shared
// between multiple images, e.g., if the same "physical" camera took multiple
// pictures with the exact same lens and intrinsics (focal length, etc.).
// This class has a specific distortion model defined by a camera model class.

constexpr int MIN_NUM_IMAGES_FOR_UPGRADE = 1;

class Camera {
 public:
  Camera();

  // Access the unique identifier of the camera.
  inline camera_t CameraId() const;
  inline void SetCameraId(const camera_t camera_id);

  // Access the camera model.
  inline int ModelId() const;
  std::string ModelName() const;
  void SetModelId(const int model_id);
  void SetModelIdFromName(const std::string& model_name);

  // Access dimensions of the camera sensor.
  inline size_t Width() const;
  inline size_t Height() const;
  inline void SetWidth(const size_t width);
  inline void SetHeight(const size_t height);

  // Access focal length parameters.
  double MeanFocalLength() const;
  double FocalLength() const;
  double FocalLengthX() const;
  double FocalLengthY() const;
  void SetFocalLength(const double focal_length);
  void SetFocalLengthX(const double focal_length_x);
  void SetFocalLengthY(const double focal_length_y);

  // Check if camera has prior focal length.
  inline bool HasPriorFocalLength() const;
  inline void SetPriorFocalLength(const bool prior);

  // Access principal point parameters. Only works if there are two
  // principal point parameters.
  double PrincipalPointX() const;
  double PrincipalPointY() const;
  void SetPrincipalPointX(const double ppx);
  void SetPrincipalPointY(const double ppy);

  // Get the indices of the parameter groups in the parameter vector.
  const std::vector<size_t>& FocalLengthIdxs() const;
  const std::vector<size_t>& PrincipalPointIdxs() const;
  const std::vector<size_t>& ExtraParamsIdxs() const;
  const std::vector<size_t>& XParamsIdx() const;

  // Get intrinsic calibration matrix composed from focal length and principal
  // point parameters, excluding distortion parameters.
  Eigen::Matrix3d CalibrationMatrix() const;

  // Get human-readable information about the parameter vector ordering.
  std::string ParamsInfo() const;

  // Access the raw parameter vector.
  inline size_t NumParams() const;
  inline const std::vector<double>& Params() const;
  inline std::vector<double>& Params();
  inline double Params(const size_t idx) const;
  inline double& Params(const size_t idx);
  inline const double* ParamsData() const;
  inline double* ParamsData();
  inline void SetParams(const std::vector<double>& params);
  inline void UpdateParams();

  // Concatenate parameters as comma-separated list.
  std::string ParamsToString() const;

  // Set camera parameters from comma-separated list.
  bool SetParamsFromString(const std::string& string);

  // Check whether parameters are valid, i.e. the parameter vector has
  // the correct dimensions that match the specified camera model.
  bool VerifyParams() const;

  // Check whether camera has bogus parameters.
  bool HasBogusParams(const double min_focal_length_ratio,
                      const double max_focal_length_ratio,
                      const double max_extra_param) const;

  // Initialize parameters for given camera model and focal length, and set
  // the principal point to be the image center.
  void InitializeWithId(const int model_id, const double focal_length,
                        const size_t width, const size_t height);
  void InitializeWithName(const std::string& model_name,
                          const double focal_length, const size_t width,
                          const size_t height);

  // Project point in image plane to world / infinity.
  Eigen::Vector2d ImageToWorld(const Eigen::Vector2d& image_point) const;

  // Convert pixel threshold in image plane to world space.
  double ImageToWorldThreshold(const double threshold) const;

  // Project point from world / infinity to image plane.
  Eigen::Vector2d WorldToImage(const Eigen::Vector2d& world_point) const;

  // Rescale camera dimensions and accordingly the focal length and
  // and the principal point.
  void Rescale(const double scale);
  void Rescale(const size_t width, const size_t height);



  ////// --------------------% TODO: Add the following functions into new camera model %------------------ //////
  inline std::vector<double> GetFocalLengthParams() const;
  inline void SetFocalLengthParams(const std::vector<double>& focal_length_params) const;
  inline std::vector<double> GetRawRadii() const; 
  inline tk::spline<double> GetSpline() const;
  inline std::vector<std::vector<double>> GetIntervals() const;
  inline void SetRawRadii(const std::vector<double>& raw_radii) const;
  inline std::vector<double> GetTheta() const;
  inline void SetTheta(const std::vector<double>& theta) const;
  inline void recursiveSplit(const std::vector<double>& radii, const std::vector<double>& focal_lengths,
                        std::vector<std::vector<double>>& radii_segments,
                        std::vector<std::vector<double>>& focal_lengths_segments, 
                        double threshold = 0.05, double stddev_threshold = 0.01) const;

  bool FitPIeceWiseSpline_binary(std::vector<double>& radii, std::vector<double>& focal_lengths, std::vector<double>& principal_point) ;
  inline void FitSpline(std::vector<double>& radii,  std::vector<double>& focal_lengths) const ;
  inline void FitSpline_theta_r(std::vector<double>& radii,  std::vector<double>& focal_lengths, std::vector<double> & principle_point) ;
  inline void FitPieceWiseSpline(std::vector<double>& radii,  std::vector<double>& focal_lengths) const ;
  inline void FitGridSpline(std::vector<double>& radii, std::vector<double>& focal_lengths) const;
  inline double EvalFocalLength(double radius) const;
  inline double EvalPieceFocalLength(double radius) const;
  inline double EvalGridFocalLength(double radius) const;

  inline double EvalFocalLength(const Eigen::Vector2d& image_point) const;
  inline double EvalFocalLength(const Eigen::Vector3d& point3d) const;
  inline bool IsFullyCalibrated(const Eigen::Vector2d& image_point) const;
  inline void SetCalibrated(bool calibrated);
  inline bool IsCalibrated() const;
  inline bool SetSplineFromParams();

 private:
  // The unique identifier of the camera. If the identifier is not specified
  // it is set to `kInvalidCameraId`.
  camera_t camera_id_;

  // The identifier of the camera model. If the camera model is not specified
  // the identifier is `kInvalidCameraModelId`.
  int model_id_;

  // The dimensions of the image, 0 if not initialized.
  size_t width_;
  size_t height_;

  // The focal length, principal point, and extra parameters. If the camera
  // model is not specified, this vector is empty.
  std::vector<double> params_;

  // Whether there is a safe prior for the focal length,
  // e.g. manually provided or extracted from EXIF
  bool prior_focal_length_;

  bool is_fully_calibrated_ = false;

  ////// --------------------% TODO: Add the following functions into new camera model %------------------ //////
  mutable std::vector<double> focal_length_params_;
  mutable std::vector<double> raw_radii_;
  mutable std::vector<double> theta_;
  mutable tk::spline<double> spline_;
  mutable std::vector<tk::spline<double>> piece_splines_; 
  mutable std::vector<tk::spline<double>> grid_splines_;
  mutable std::vector<double> updated_params_;
  mutable std::vector<std::vector<double>> intervals_;
  mutable std::vector<double> grids_; 
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

camera_t Camera::CameraId() const { return camera_id_; }

void Camera::SetCameraId(const camera_t camera_id) { camera_id_ = camera_id; }

int Camera::ModelId() const { return model_id_; }

size_t Camera::Width() const { return width_; }

size_t Camera::Height() const { return height_; }

void Camera::SetWidth(const size_t width) { width_ = width; }

void Camera::SetHeight(const size_t height) { height_ = height; }

bool Camera::HasPriorFocalLength() const { return prior_focal_length_; }

void Camera::SetPriorFocalLength(const bool prior) {
  prior_focal_length_ = prior;
}

size_t Camera::NumParams() const { return params_.size(); }

const std::vector<double>& Camera::Params() const { return params_; }

std::vector<double>& Camera::Params() { return params_; }

double Camera::Params(const size_t idx) const { return params_[idx]; }

double& Camera::Params(const size_t idx) { return params_[idx]; }

const double* Camera::ParamsData() const { return params_.data(); }

double* Camera::ParamsData() { return params_.data(); }

void Camera::SetParams(const std::vector<double>& params) { params_ = params; }
void Camera::UpdateParams() { 
  
  if(ModelId()==ImplicitDistortionModel::model_id &&!updated_params_.empty()){
  params_ = updated_params_; }
  }

////// --------------------% TODO: Add the following functions into new camera model %------------------ //////

inline std::vector<double> Camera::GetFocalLengthParams() const {return focal_length_params_;}
inline void Camera::SetFocalLengthParams(const std::vector<double>& focal_length_params) const {focal_length_params_ = focal_length_params;}
inline void Camera::FitSpline(std::vector<double>& radii, std::vector<double>& focal_lengths) const{
  // Convert std::vector to Eigen vectors
  assert(radii.size() == focal_lengths.size());

  std::vector<int> increasing_indices;
  for (int i = 0; i < radii.size() - 1; i++) {
    if(radii[i+1] > radii[i]) {
      increasing_indices.push_back(i);
    }
  }
  std::vector<double> new_radii;
  std::vector<double> new_focal_lengths;
  if (!increasing_indices.empty() && increasing_indices.back() != radii.size() - 1) {
        increasing_indices.push_back(radii.size() - 1);
    }

  for (int i = 0; i < increasing_indices.size(); i++) {
    new_radii.push_back(radii[increasing_indices[i]]);
    new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
  }
  // std::cout << "new_radii size: " << new_radii.size() << std::endl;
  // std::cout << "new_focal_lengths size: " << new_focal_lengths.size() << std::endl;
  
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, new_radii.size() - 1);
  tk::spline<double> best_spline;
  const int max_iterations = 80;
  const double threshold = 10.0;
  int degree = 10;
  
  for (int i = 0; i < max_iterations; ++i){
    std::vector<std::pair<double, double>> samples;
    for (int j = 0; j < degree; ++j) {
      int idx = dis(gen);
      samples.emplace_back(new_radii[idx], new_focal_lengths[idx]);
    }
    std::sort(samples.begin(), samples.end());
    std::vector<double> sample_x, sample_y;
        for (const auto& pair : samples) {
            sample_x.push_back(pair.first);
            sample_y.push_back(pair.second);
        }
    tk::spline<double> s;
    s.set_points(sample_x, sample_y);
    int inliers_count = 0;
    for (size_t j = 0; j < new_radii.size(); ++j) {
      // std::cout << "new_radii[j]: " << new_radii[j] << std::endl;
      // check if s is empty
      // check the range of sample_x
      auto min_max_x = std::minmax_element(sample_x.begin(), sample_x.end());
      auto x = s.get_x();
      auto min_max_x_s = std::minmax_element(x.begin(), x.end());
      double y_est;
      try{
        y_est = s(new_radii[j]);
      }catch(const std::exception& e){
        std::cout << "Exception caught: " << e.what() << std::endl;
        // y_est = new_focal_lengths[j];
        continue;
      }
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            ++inliers_count;
        }
    }
    if (inliers_count >= best_inliers_count) {
            best_inliers_count = inliers_count;
            best_spline = s;
        }
  }

  std::vector<double> inliers_x, inliers_y;
    for (size_t j = 0; j < new_radii.size(); ++j) {
      
        double y_est = best_spline(new_radii[j]);
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            inliers_x.push_back(new_radii[j]);
            inliers_y.push_back(new_focal_lengths[j]);
        }
    }
  spline_ = best_spline;

  std::vector<double> used_x;
  std::vector<double> used_y;
  if (ModelId() == ImplicitDistortionModel::model_id) {
    used_x = best_spline.get_x();
    used_y = best_spline.get_y();
    std::vector<double> updated_params;
    updated_params.push_back(PrincipalPointX());
    updated_params.push_back(PrincipalPointY());
    for (int i = 0; i < degree; i++) {
      updated_params.push_back(used_x[i]);
    }
    for(int i = 0; i < degree; i++) {
      updated_params.push_back(used_y[i]);
    }
    // updated_params_ = updated_params;
    // SetParams(updated_params);
  }
}

inline void Camera::FitSpline_theta_r(std::vector<double>& radii, std::vector<double>& focal_lengths, std::vector<double> & principle_point) {
  // Convert std::vector to Eigen vectors
  assert(radii.size() == focal_lengths.size());

  std::vector<int> increasing_indices;
  for (int i = 0; i < radii.size() - 1; i++) {
    if(radii[i+1] > radii[i]) {
      increasing_indices.push_back(i);
    }
  }
  std::vector<double> new_radii;
  std::vector<double> new_focal_lengths;
  if (!increasing_indices.empty() && increasing_indices.back() != radii.size() - 1) {
        increasing_indices.push_back(radii.size() - 1);
    }

  for (int i = 0; i < increasing_indices.size(); i++) {
    new_radii.push_back(radii[increasing_indices[i]]);
    new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
  }
  // std::cout << "new_radii size: " << new_radii.size() << std::endl;
  // std::cout << "new_focal_lengths size: " << new_focal_lengths.size() << std::endl;
  
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, new_radii.size() - 1);
  tk::spline<double> best_spline;
  const int max_iterations = 80;
  const double threshold =5.0;
  int degree = 10;
  std::vector<int> indices ={};
  for (int i = 0; i < max_iterations; ++i){
    std::vector<std::pair<double, double>> samples;
    indices={}; 
    for (int j = 0; j < degree; ++j) {
      int idx = dis(gen);
      // check if the idx has already been explored
      while(std::find(indices.begin(), indices.end(), idx) != indices.end()){
        idx = dis(gen);
      }
      indices.push_back(idx);
      samples.emplace_back(new_radii[idx], new_focal_lengths[idx]);
    }
    std::sort(samples.begin(), samples.end());
    std::vector<double> sample_x, sample_y;
        for (const auto& pair : samples) {
            sample_x.push_back(pair.first);
            sample_y.push_back(pair.second);
        }
    tk::spline<double> s;
    s.set_points(sample_x, sample_y);
    int inliers_count = 0;
    for (size_t j = 0; j < new_radii.size(); ++j) {
      // std::cout << "new_radii[j]: " << new_radii[j] << std::endl;
      // check if s is empty
      // check the range of sample_x
      auto min_max_x = std::minmax_element(sample_x.begin(), sample_x.end());
      auto x = s.get_x();
      auto min_max_x_s = std::minmax_element(x.begin(), x.end());
      double y_est;
      try{
        y_est = s(new_radii[j]);
      }catch(const std::exception& e){
        std::cout << "Exception caught: " << e.what() << std::endl;
        // y_est = new_focal_lengths[j];
        continue;
      }
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            ++inliers_count;
        }
    }
    if (inliers_count >= best_inliers_count) {
            best_inliers_count = inliers_count;
            best_spline = s;
        }
  }

  std::vector<double> inliers_x, inliers_y;
    for (size_t j = 0; j < new_radii.size(); ++j) {
      
        double y_est = best_spline(new_radii[j]);
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            inliers_x.push_back(new_radii[j]);
            inliers_y.push_back(new_focal_lengths[j]);
        }
    }
  // spline_ = best_spline;

  std::vector<double> used_x;
  std::vector<double> used_y;
  if (ModelId() == ImplicitDistortionModel::model_id) {
    updated_params_={};
    used_x = best_spline.get_x();
    used_y = best_spline.get_y();
    std::vector<double> updated_params;
    updated_params.push_back(principle_point[0]);
    updated_params.push_back(principle_point[1]);
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
}
inline void Camera::recursiveSplit(const std::vector<double>& radii, const std::vector<double>& focal_lengths,
                        std::vector<std::vector<double>>& radii_segments,
                        std::vector<std::vector<double>>& focal_lengths_segments, 
                        double threshold, double stddev_threshold)const {
        // Find the largest interval
        double max_gap = 0;
        size_t index_of_max_gap = 0;
        for (size_t i = 0; i < radii.size() - 1; i++) {
            double gap = radii[i + 1] - radii[i];
            if (gap > max_gap) {
                max_gap = gap;
                index_of_max_gap = i;
            }
        }

        // Calculate the standard deviation of the intervals
        double mean = max_gap / (radii.size() - 1);
        double sum_sq_diff = 0;
        for (size_t i = 0; i < radii.size() - 1; i++) {
            double diff = radii[i + 1] - radii[i] - mean;
            sum_sq_diff += diff * diff;
        }
        double stddev = std::sqrt(sum_sq_diff / (radii.size() - 1));

        // Base case for recursion
        // if (max_gap < threshold || stddev < stddev_threshold) {
        if (radii.size() < NUM_CONTROL_POINTS)
          return;
        if (max_gap < threshold) {
            radii_segments.push_back(radii);
            focal_lengths_segments.push_back(focal_lengths);
            return;
        }

        // Recursive case: split at the largest gap
        std::vector<double> left_radii(radii.begin(), radii.begin() + index_of_max_gap + 1);
        std::vector<double> left_focal_lengths(focal_lengths.begin(), focal_lengths.begin() + index_of_max_gap + 1);
        std::vector<double> right_radii(radii.begin() + index_of_max_gap + 1, radii.end());
        std::vector<double> right_focal_lengths(focal_lengths.begin() + index_of_max_gap + 1, focal_lengths.end());

        recursiveSplit(left_radii, left_focal_lengths, radii_segments, focal_lengths_segments, threshold, stddev_threshold);
        recursiveSplit(right_radii, right_focal_lengths, radii_segments, focal_lengths_segments, threshold, stddev_threshold);
    }
inline void Camera::FitPieceWiseSpline(std::vector<double>& radii, std::vector<double>& focal_lengths) const{
  // Convert std::vector to Eigen vectors
  assert(radii.size() == focal_lengths.size());

  std::vector<int> increasing_indices;
  for (int i = 0; i < radii.size() - 1; i++) {
    if(radii[i+1] > radii[i]) {
      increasing_indices.push_back(i);
    }
  }
  std::vector<double> new_radii;
  std::vector<double> new_focal_lengths;
  if (!increasing_indices.empty() && increasing_indices.back() != radii.size() - 1) {
        increasing_indices.push_back(radii.size() - 1);
    }

  for (int i = 0; i < increasing_indices.size(); i++) {
    new_radii.push_back(radii[increasing_indices[i]]);
    new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
  }
//   std::cout << "new_radii size: " << new_radii.size() << std::endl;
// //   for (int i = 0; i < new_radii.size(); i++) {
// //     std::cout << "new_radii: " << new_radii[i] << std::endl;
// //   }
//   std::cout << "new_focal_lengths size: " << new_focal_lengths.size() << std::endl;
  
  
  // use ransac technique to fit the spline
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, new_radii.size() - 1);
  tk::spline<double> best_spline;
  const int max_iterations = 40;
  const double threshold = 20.0;
  int degree = 10;
  std::vector<int> indices;
  
  for (int i = 0; i < max_iterations; ++i){
    std::vector<std::pair<double, double>> samples;
    indices.clear();
    for (int j = 0; j < degree; ++j) {
      int idx = dis(gen);
    //   ensure that the samples are unique
        while(std::find(indices.begin(), indices.end(), idx) != indices.end()){
            idx = dis(gen);
        }
        indices.push_back(idx);
      samples.emplace_back(new_radii[idx], new_focal_lengths[idx]);
    }
    std::sort(samples.begin(), samples.end());
    std::vector<double> sample_x, sample_y;
        for (const auto& pair : samples) {
            sample_x.push_back(pair.first);
            // std::cout << "sample_x: " << pair.first << std::endl;
            sample_y.push_back(pair.second);
        }
    tk::spline<double> s;
    s.set_points(sample_x, sample_y);
    // std::cout << "s fitted" << std::endl;
    int inliers_count = 0;
    for (size_t j = 0; j < new_radii.size(); ++j) {
    
      auto min_max_x = std::minmax_element(sample_x.begin(), sample_x.end());
      auto x = s.get_x();
      auto min_max_x_s = std::minmax_element(x.begin(), x.end());
     
        double y_est;
        try{
         y_est = s(new_radii[j]);
        }catch(const std::exception& e){
          std::cout << "Exception caught: " << e.what() << std::endl;
          // y_est = new_focal_lengths[j];
          continue;
        }
        
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            ++inliers_count;
        }
    }
    if (inliers_count >= best_inliers_count) {
            best_inliers_count = inliers_count;
            best_spline = s;
        }
  }

  std::vector<double> inliers_x, inliers_y;
    for (size_t j = 0; j < new_radii.size(); ++j) {
      
        double y_est = best_spline(new_radii[j]);
 
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            inliers_x.push_back(new_radii[j]);
            inliers_y.push_back(new_focal_lengths[j]);
        }
    }

  spline_= best_spline;
  // calculate an confidence band for the spline based on the distance from the data points
  std::vector<double> errors;
  double mean_error=0.0;
  double std_error=0.0;
  for (size_t j = 0; j < new_radii.size(); ++j) {
    double y_est = best_spline(new_radii[j]);
    errors.push_back(fabs(y_est - new_focal_lengths[j]));
    mean_error += fabs(y_est - new_focal_lengths[j]);
    std_error += pow(fabs(y_est - new_focal_lengths[j]), 2);
  }
  mean_error /= new_radii.size();
  // std::cout << "mean_error: " << mean_error << std::endl;
  // std::cout << "std_error before division: " << std_error <<"new_radii size:"<<new_radii.size() << std::endl;
  std_error = sqrt(std_error / new_radii.size() );
  // std::cout << "std_error: " << std_error << std::endl;   

  double outlier_threshold = mean_error + 0.5*std_error;
  // Determinig the intervals where a new piece-wise spline should be fitted
  std::vector<double> outlier_xs;
  for (size_t j = 0; j < new_radii.size(); ++j) {
    if(errors[j] > outlier_threshold){
      outlier_xs.push_back(new_radii[j]);
    }
  }
  // std::cout << "outlier_xs size: " << outlier_xs.size() << std::endl;
  std::vector<std::vector<double>> new_radii_segments;
  new_radii_segments.clear();
  std::vector<std::vector<double>> new_focal_lengths_segments;
  new_focal_lengths_segments.clear();

bool start_new_spline = false;
bool update_segment = false;    
std::vector<double> segment_radii ;
segment_radii.clear();
std::vector<double> segment_focal_lengths;
segment_focal_lengths.clear();
// std::vector<double> new_radii_segments;
// std::vector<double> new_focal_lengths_segments;
const double outlier_fraction_threshold = 0.3;
int outlier_count = 0;
int segment_count = 0;
int inlier_count = 0;   
std::vector<double> outliers_indices;
outliers_indices.clear();
for (int i = 0; i < new_radii.size(); i++) {
  if(abs(best_spline(new_radii[i]) - new_focal_lengths[i]) > outlier_threshold){
    outliers_indices.push_back(i);
  }
}
for(int i = 0; i < outliers_indices.size()-1; i++){
    if(i == 0 && outliers_indices[i] != 0){
        for (int j = 0; j < outliers_indices[i]; j++) {
            segment_radii.push_back(new_radii[j]);
            segment_focal_lengths.push_back(new_focal_lengths[j]);
        }
        new_radii_segments.push_back(segment_radii);
        new_focal_lengths_segments.push_back(segment_focal_lengths);
        segment_radii.clear();
        segment_focal_lengths.clear();
    }
    if(outliers_indices[i+1] - outliers_indices[i] <5){
        segment_radii.push_back(new_radii[outliers_indices[i]]);
        segment_focal_lengths.push_back(new_focal_lengths[outliers_indices[i]]);
        start_new_spline = true;
        update_segment = true;
    }else{
        new_radii_segments.push_back(segment_radii);
        new_focal_lengths_segments.push_back(segment_focal_lengths);
        segment_radii.clear();
        segment_focal_lengths.clear();
        for(int j = outliers_indices[i]; j < outliers_indices[i+1]; j++){
            segment_radii.push_back(new_radii[j]);
            segment_focal_lengths.push_back(new_focal_lengths[j]);
        }
        new_radii_segments.push_back(segment_radii);
        new_focal_lengths_segments.push_back(segment_focal_lengths);
        segment_radii.clear();
        segment_focal_lengths.clear();
    }
}
// std::cout << "new_radii_segments size: " << new_radii_segments.size() << std::endl;

//   populate the intervals_
    intervals_.clear();
    for (size_t i = 0; i < new_radii_segments.size(); ++i) {
      if (new_radii_segments[i].empty()) {
        std::cerr << "Error: new_radii_segments[" << i << "] is empty." << std::endl;
        continue; // Skip empty segments
    }
        std::vector<double> interval;
        interval.push_back(new_radii_segments[i].front());
        interval.push_back(new_radii_segments[i].back());
        intervals_.push_back(interval);
    }
    // std::cout << "intervals_ set: "  << std::endl; 
    if(intervals_.empty()){
        piece_splines_ = {best_spline};
    }

  // Fit the piece-wise splines
  // std::vector<<tk::spline<double>>> splines;
  std::vector<tk::spline<double>> splines;
  int piecewise_max_it = 10;
  for (size_t i = 0; i < new_radii_segments.size(); ++i) {
      // tk::spline s;
      std::vector<int>piece_indices;

      // s.set_points(new_radii_segments[i], new_focal_lengths_segments[i]);
      // splines.push_back(s);
      if (new_radii_segments[i].size() > degree) {
          // using ransac technique to fit the spline
          std::uniform_int_distribution<> dis_piece(0, new_radii_segments[i].size() - 1);
          int best_piece_inliers_count = 0;
          tk::spline<double> best_piece_spline;
          for (int j = 0; j < piecewise_max_it; ++j) {
              std::vector<std::pair<double, double>> samples_piece;
              piece_indices.clear();
              for (int k = 0; k < degree; ++k) {
                  int idx = dis_piece(gen);
                    while(std::find(piece_indices.begin(), piece_indices.end(), idx) != piece_indices.end()){
                        idx = dis_piece(gen);
                    }
                    piece_indices.push_back(idx);
                    // std::cout << "idx: " << idx << std::endl;
                    // std::cout << "new_radii_segments[i].size: " << new_radii_segments[i].size() << std::endl;
                  samples_piece.emplace_back(new_radii_segments[i][idx], new_focal_lengths_segments[i][idx]);
              }
              std::sort(samples_piece.begin(), samples_piece.end());
              for (int k = 0; k < piece_indices.size(); k++) {
                  // std::cout << "new_radii_segments[i].size: " << new_radii_segments[i].size() << std::endl;
                  // std::cout << "piece_indices: " << piece_indices[k] << std::endl;
                  // std::cout << "samples_piece: " << samples_piece[k].first << std::endl;

              }
              std::vector<double> sample_x_piece, sample_y_piece;
              for (const auto& pair : samples_piece) {
                  sample_x_piece.push_back(pair.first);
                  sample_y_piece.push_back(pair.second);
              }
              tk::spline<double> s_piece;

              s_piece.set_points(sample_x_piece, sample_y_piece);
              // std::cout << "s_piece fitted" << std::endl;
              int piece_inliers_count = 0;
              for (size_t j = 0; j < new_radii_segments[i].size(); ++j) {
                  double y_est_piece = s_piece(new_radii_segments[i][j]);
                  if (fabs(y_est_piece - new_focal_lengths_segments[i][j]) < outlier_threshold) {
                      ++piece_inliers_count;
                  }
              }
              if (piece_inliers_count >= best_piece_inliers_count) {
                  best_piece_inliers_count = piece_inliers_count;
                  best_piece_spline = s_piece;
              }
          }
          splines.push_back(best_piece_spline);
      }else if (new_radii_segments[i].size() <= degree && new_radii_segments[i].size() > 2){
        tk::spline<double> s;
        s.set_points(new_radii_segments[i], new_focal_lengths_segments[i]);
        splines.push_back(s);
      }else{
        std::cout << "Not enough data points to fit a spline" << std::endl;
        continue;
      }
  }

  if(!splines.empty()){
    piece_splines_ = splines;
}
  std::cout << "piece_splines_ size: " << piece_splines_.size() << std::endl;
 
}

// TODO: THIS FUNCTION IS INCORRECT
inline double Camera::EvalPieceFocalLength(double radius) const {
  if (intervals_.empty() ) {
      return spline_(radius);
    }
    if(radius <= intervals_.front()[0]){
        return piece_splines_.front()(radius);
    }else if (radius >= intervals_.back()[1])
    {
        return piece_splines_.back()(radius);
    }else{
    
    for (size_t i = 0; i < intervals_.size(); ++i) {
        if (radius >= intervals_[i][0] && radius <= intervals_[i][1]) {
            if (piece_splines_.size() > (i+1)) {
                return piece_splines_[i](radius);
              }else{
                return piece_splines_.back()(radius);
              }
              }
    }
    }
    return spline_(radius);

}

inline void Camera::FitGridSpline(std::vector<double>& radii, std::vector<double>& focal_lengths) const{
// Convert std::vector to Eigen vectors
  assert(radii.size() == focal_lengths.size());

  std::vector<int> increasing_indices;
  for (int i = 0; i < radii.size() - 1; i++) {
    if(radii[i+1] > radii[i]) {
      increasing_indices.push_back(i);
    }
  }
  std::vector<double> new_radii;
  std::vector<double> new_focal_lengths;
  if (!increasing_indices.empty() && increasing_indices.back() != radii.size() - 1) {
        increasing_indices.push_back(radii.size() - 1);
    }

  for (int i = 0; i < increasing_indices.size(); i++) {
    new_radii.push_back(radii[increasing_indices[i]]);
    new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
  }
  // std::cout << "new_radii size: " << new_radii.size() << std::endl;
//   for (int i = 0; i < new_radii.size(); i++) {
//     std::cout << "new_radii: " << new_radii[i] << std::endl;
//   }
  // std::cout << "new_focal_lengths size: " << new_focal_lengths.size() << std::endl;
  
  
  // use ransac technique to fit the spline
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, new_radii.size() - 1);
  tk::spline<double> best_spline;
  const int max_iterations = 40;
  const double threshold = 20.0;
  int degree = 10;
  std::vector<int> indices;
  
  for (int i = 0; i < max_iterations; ++i){
    std::vector<std::pair<double, double>> samples;
    indices.clear();
    for (int j = 0; j < degree; ++j) {
      int idx = dis(gen);
    //   ensure that the samples are unique
        while(std::find(indices.begin(), indices.end(), idx) != indices.end()){
            idx = dis(gen);
        }
        indices.push_back(idx);
      samples.emplace_back(new_radii[idx], new_focal_lengths[idx]);
    }
    std::sort(samples.begin(), samples.end());
    std::vector<double> sample_x, sample_y;
        for (const auto& pair : samples) {
            sample_x.push_back(pair.first);
            // std::cout << "sample_x: " << pair.first << std::endl;
            sample_y.push_back(pair.second);
        }
    tk::spline<double> s;
    s.set_points(sample_x, sample_y);
    // std::cout << "s fitted" << std::endl;
    int inliers_count = 0;
    for (size_t j = 0; j < new_radii.size(); ++j) {
    
      auto min_max_x = std::minmax_element(sample_x.begin(), sample_x.end());
      auto x = s.get_x();
      auto min_max_x_s = std::minmax_element(x.begin(), x.end());
     
        double y_est;
        try{
         y_est = s(new_radii[j]);
        }catch(const std::exception& e){
          std::cout << "Exception caught: " << e.what() << std::endl;
          // y_est = new_focal_lengths[j];
          continue;
        }
        
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            ++inliers_count;
        }
    }
    if (inliers_count >= best_inliers_count) {
            best_inliers_count = inliers_count;
            best_spline = s;
        }
  }

  std::vector<double> inliers_x, inliers_y;
    for (size_t j = 0; j < new_radii.size(); ++j) {
      
        double y_est = best_spline(new_radii[j]);
 
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            inliers_x.push_back(new_radii[j]);
            inliers_y.push_back(new_focal_lengths[j]);
        }
    }

  std::vector<double> grid_points;
  grid_points.clear();
  // grid_points.push_back(new_radii[0]);
  double interval = (new_radii.back() - new_radii[0])/4;
  double mean_count = new_radii.size()/4; 
  double radial_threshold = 0.3*mean_count;
  int grid_max_it = 10;
  int grid_degree = 10;
  int grid_threshold = 2;
  std::vector<double> gird_count = {0, 0, 0, 0};
  for (size_t i = 0; i < new_radii.size(); ++i) {
    if(new_radii[i] < new_radii[0] + interval){
      gird_count[0] += 1;
    }else if(new_radii[i] < new_radii[0] + 2*interval){
      gird_count[1] += 1;
    }else if(new_radii[i] < new_radii[0] + 3*interval){
      gird_count[2] += 1;
    }else{
      gird_count[3] += 1;
    }
  }
  // std::cout << "grid_count: " << gird_count[0] << " " << gird_count[1] << " " << gird_count[2] << " " << gird_count[3] << std::endl;
  // std::cout << "grid points:"<< new_radii[0] << " " << new_radii[0] + interval << " " << new_radii[0] + 2*interval << " " << new_radii[0] + 3*interval << " " << new_radii[0] + 4*interval << std::endl;
  std::vector<bool> calibrated = {true, true, true, true};
  for (int i = 0; i < gird_count.size(); ++i) {
    if(gird_count[i] < radial_threshold || gird_count[i] < 3){
      calibrated[i] = false;
    }
  }
  for (int i = 0; i < calibrated.size(); ++i) {
    if(i == 0 && calibrated[i]){
      grid_points.push_back(new_radii[0]);
    }
    if(calibrated[i]){
      grid_points.push_back(new_radii[0] + (i+1)*interval);
    }
  }
  grids_ = grid_points;
  // fit a separate spline for each grid point
  if (grid_points.size() > 1) {
    std::vector<tk::spline<double>> grid_splines;
    for (size_t i = 0; i < grid_points.size() - 1; ++i) {
      std::vector<double> grid_radii;
      std::vector<double> grid_focal_lengths;
      for (size_t j = 0; j < new_radii.size(); ++j) {
        if (new_radii[j] >= grid_points[i] && new_radii[j] < grid_points[i+1]) {
          grid_radii.push_back(new_radii[j]);
          grid_focal_lengths.push_back(new_focal_lengths[j]);
        }
      }
      if(grid_radii.size() >= grid_degree){
      // using ransac technique to fit the spline
      std::uniform_int_distribution<> dis_grid(0, grid_radii.size() - 1);
      int best_grid_inliers_count = 0;
      tk::spline<double> best_grid_spline;
      std::vector<int> grid_indices;
      for (int j = 0; j < grid_max_it; ++j) {
        std::vector<std::pair<double, double>> samples_grid;
        grid_indices.clear();
        std::cout<<"sample begains" << std::endl;
        for (int k = 0; k < grid_degree; ++k) {
          int idx = dis_grid(gen);
          while(std::find(grid_indices.begin(), grid_indices.end(), idx) != grid_indices.end()){
            idx = dis_grid(gen);
          }
          // std::cout << "k: " << k<<" idx: "<<idx << std::endl;
          grid_indices.push_back(idx);
          samples_grid.emplace_back(grid_radii[idx], grid_focal_lengths[idx]);
        }
        std::sort(samples_grid.begin(), samples_grid.end());
        std::vector<double> sample_x_grid, sample_y_grid;
        for (const auto& pair : samples_grid) {
          sample_x_grid.push_back(pair.first);
          sample_y_grid.push_back(pair.second);
        }
        tk::spline<double> s_grid;
        s_grid.set_points(sample_x_grid, sample_y_grid);
        int grid_inliers_count = 0;
        for (size_t k = 0; k < grid_radii.size(); ++k) {
          double y_est_grid = s_grid(grid_radii[k]);
          if (fabs(y_est_grid - grid_focal_lengths[k]) < grid_threshold) {
            ++grid_inliers_count;
          }
        }
        if (grid_inliers_count >= best_grid_inliers_count) {
          best_grid_inliers_count = grid_inliers_count;
          best_grid_spline = s_grid;
        }
      }grid_splines.push_back(best_grid_spline);
      }else{
        tk::spline<double> best_grid_spline;
        best_grid_spline.set_points(grid_radii, grid_focal_lengths);
        grid_splines.push_back(best_grid_spline);  
      }
      
    }
    grid_splines_ = grid_splines;
  }else{
    grid_splines_.push_back(best_spline); 
  }
  // std::cout << "grid_splines_ size: " << grid_splines_.size() << std::endl;  
}

// TODO: THIS FUNCTION IS INCORRECT
inline double Camera::EvalGridFocalLength(double radius) const{
  if(grids_.empty() ){
    return spline_(radius);
  }else{
    if(radius <= grids_.front()){
        return grid_splines_.front()(radius);
    }else if (radius >= grids_.back())
    {
        return grid_splines_.back()(radius);
    }else{
    
    for (size_t i = 0; i < grids_.size(); ++i) {
        if (radius >= grids_[i] && radius <= grids_[i+1]) {
            if (grid_splines_.size() > (i+1)) {
                return grid_splines_[i](radius);
              }else{
                return grid_splines_.back()(radius);
              }
              }
    }
    }
    return spline_(radius);
  
  }
}


// double Camera::EvalFocalLength(double radius) const {return spline_(radius);}
double Camera::EvalFocalLength(double radius) const {
  if (model_id_ == Radial1DCameraModel::model_id)
    return 1.;
  else if (model_id_ != ImplicitDistortionModel::model_id) {
    return FocalLength();
  } else {
    if (!is_fully_calibrated_)
      return 1.;

    int num_control_points = (ImplicitDistortionModel::kNumParams - 2) / 2;
    // First check whether it is within the calibrated region
    if (radius < params_[2 + num_control_points] || radius > params_[2 * num_control_points + 1]) {
      return 1.;
    }

    auto it = std::upper_bound(params_.begin() + 2 + num_control_points, params_.end(), radius);
    size_t idx = std::max<int>((it - params_.begin()) - 1 - 2 + num_control_points, 0); // (params[idx + 2 + num_control_points] <= radius)

    double theta = params_[2 + idx];
    double residual = spline_(theta) - radius;
    int num_iter = 0;
    while (abs(residual) > 1e-6 && num_iter < 10) {
      // Use newton's method to find the theta
      theta = theta - residual / spline_.deriv(1, theta);
      residual = spline_(theta) - radius;
      num_iter++;
    }

    // Convert the theta to focal length
    return radius / std::tan(theta);
  }
  return 1.;
}
inline std::vector<double> Camera::GetRawRadii() const {return raw_radii_;}
inline void Camera::SetRawRadii(const std::vector<double>& raw_radii) const {raw_radii_ = raw_radii;}
inline std::vector<double> Camera::GetTheta() const {return theta_;}
inline void Camera::SetTheta(const std::vector<double>& theta) const {theta_ = theta;}
inline tk::spline<double> Camera::GetSpline() const {return spline_;}
inline std::vector<std::vector<double>> Camera::GetIntervals() const {return intervals_;}

inline double Camera::EvalFocalLength(const Eigen::Vector2d& image_point) const {
  return EvalFocalLength(ImageToWorld(image_point).norm());
};

inline double Camera::EvalFocalLength(const Eigen::Vector3d& point3d) const {
  double theta = std::atan2(point3d.topRows<2>().norm(), point3d[2]);
  double radius = spline_(theta);
  return radius / std::tan(theta);
}

// 
inline bool Camera::IsFullyCalibrated(const Eigen::Vector2d& image_point) const {
  // If it is with Radial1DCameraModel, it is not fully calibrated
  if (model_id_ == Radial1DCameraModel::model_id)
    return false;
  // If it is parametric camera model, it is fully calibrated
  if (model_id_ != ImplicitDistortionModel::model_id)
    return true;

  // If it is the case with ImplicitDistortionModel, need to check whether it is within calibrated region
  Eigen::Vector2d image_point_norm = ImageToWorld(image_point);
  double radius = image_point_norm.norm();

  // if(raw_radii_.size() < 20 || (raw_radii_.size() > 20 && raw_radii_[0] < 0))
  //   return false;
  int num_control_points = (ImplicitDistortionModel::kNumParams - 2) / 2;
  if (!is_fully_calibrated_ || (radius < params_[2 + num_control_points] || radius > params_[2 * num_control_points + 1]))
    return false;
  else
    return true;
}

void Camera::SetCalibrated(bool calibrated) {
  is_fully_calibrated_ = calibrated;
}

bool Camera::IsCalibrated() const {
  return (model_id_ != ImplicitDistortionModel::model_id && model_id_ != Radial1DCameraModel::model_id) 
   || is_fully_calibrated_ ;
}

bool Camera::SetSplineFromParams() {
  if (model_id_ != ImplicitDistortionModel::model_id || !is_fully_calibrated_)
    return true;
    
  std::vector<double> sample_x, sample_y;
  int num_control_points = (ImplicitDistortionModel::kNumParams - 2) / 2;
  for (int i = 2; i < 2 + num_control_points; i++) {
    sample_x.push_back(params_[i]);
  }
  
  // Enforce the monotonicity of y
  bool is_incresing = sample_x[0] < sample_x[1];
  double extreme_y = params_[2 + num_control_points];
  double diagonal = sqrt(pow(width_, 2) + pow(height_, 2)) / 2;
  for(int i = 2 + num_control_points; i < ImplicitDistortionModel::kNumParams; i++) {
    if (is_incresing) {
      extreme_y = std::max(params_[i], extreme_y + 1e-3);
    } else {
      extreme_y = std::min(params_[i], extreme_y - 1e-3);
    }
    extreme_y = std::min(extreme_y, diagonal + i * 1e-3);
    extreme_y = std::max(extreme_y, -i * 1e-3);
    params_[i] = extreme_y;
    sample_y.push_back(extreme_y);
  }
  spline_.set_points(sample_x, sample_y);
  return true;
}

}  // namespace colmap

#endif  // COLMAP_SRC_BASE_CAMERA_H_
