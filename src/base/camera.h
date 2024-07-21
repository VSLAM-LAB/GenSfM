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

#include "util/types.h"

namespace colmap {

// Camera class that holds the intrinsic parameters. Cameras may be shared
// between multiple images, e.g., if the same "physical" camera took multiple
// pictures with the exact same lens and intrinsics (focal length, etc.).
// This class has a specific distortion model defined by a camera model class.
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
  inline tk::spline GetSpline() const;
  inline void SetRawRadii(const std::vector<double>& raw_radii) const;
  inline std::vector<double> GetTheta() const;
  inline void SetTheta(const std::vector<double>& theta) const;
  inline void FitSpline(std::vector<double>& radii,  std::vector<double>& focal_lengths) const ;
  inline void FitPieceWiseSpline(std::vector<double>& radii,  std::vector<double>& focal_lengths) const ;
  inline double EvalFocalLength(double radius) const;


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

  ////// --------------------% TODO: Add the following functions into new camera model %------------------ //////
  mutable std::vector<double> focal_length_params_;
  mutable std::vector<double> raw_radii_;
  mutable std::vector<double> theta_;
  mutable tk::spline spline_;
  mutable std::vector<double> updated_params_;
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
  std::cout << "new_radii size: " << new_radii.size() << std::endl;
  std::cout << "new_focal_lengths size: " << new_focal_lengths.size() << std::endl;
  
  
  // //// uniformly sample the radii ////
  // std::vector<std::pair<double,double>> paired_data(new_radii.size());
  // for(int i = 0; i < new_radii.size(); i++){
  //   paired_data[i] = std::make_pair(new_radii[i], new_focal_lengths[i]);
  // }

  // int low_quantile_index = (0.05 * new_radii.size());
  // int high_quantile_index = (0.95 * new_radii.size());
  // std::vector<double> sample_x, sample_y;
  // int degree = 10;
  // double step = (new_radii[high_quantile_index] - new_radii[low_quantile_index]) / (degree-1);
  // for (int i = 0; i < degree; ++i) {
  //       int idx = (low_quantile_index + i * step);
  //       sample_x.push_back(new_radii[idx]);
  //       sample_y.push_back(new_focal_lengths[idx]);
  //   }

  // spline_.set_points(sample_x, sample_y);
  // tk::spline best_spline = spline_;

  // use ransac technique to fit the spline
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, new_radii.size() - 1);
  tk::spline best_spline;
  const int max_iterations = 40;
  const double threshold = 20.0;
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
    tk::spline s;
    s.set_points(sample_x, sample_y);
    int inliers_count = 0;
    for (size_t j = 0; j < new_radii.size(); ++j) {
      // std::cout << "new_radii[j]: " << new_radii[j] << std::endl;
      // check if s is empty
      // check the range of sample_x
      auto min_max_x = std::minmax_element(sample_x.begin(), sample_x.end());
      auto x = s.get_x();
      auto min_max_x_s = std::minmax_element(x.begin(), x.end());
      // std::cout << "Min x value: " << *min_max_x.first << ", Max x value: " << *min_max_x.second << std::endl;
      // std::cout << "Min x value for s: " << *min_max_x_s.first << ", Max x value for s: " << *min_max_x_s.second << std::endl;
        // if (s.get_x().empty()) {
        //     continue;
        // }
        double y_est;
        try{
         y_est = s(new_radii[j]);
        }catch(const std::exception& e){
          std::cout << "Exception caught: " << e.what() << std::endl;
          // y_est = new_focal_lengths[j];
          continue;
        }
        
        // double y_est = s(new_radii[j]);
        // std::cout << "j: " << j << "new_radii[j]" << new_radii[j]<< "new_focal_lengths[j]"<<new_focal_lengths[j]<< std::endl;
        // std::cout << "y_est for s: " << y_est << std::endl;
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
        // std::cout << "y_est for best_spline: " << y_est << std::endl;
        // std::cout << "new_radii[j]: " << new_radii[j] << std::endl;
        // std::cout << "new_focal_lengths[j]: " << new_focal_lengths[j] << std::endl;
        if (fabs(y_est - new_focal_lengths[j]) < threshold) {
            inliers_x.push_back(new_radii[j]);
            inliers_y.push_back(new_focal_lengths[j]);
        }
    }
  // if(inliers_x.size() >= 4){
  //  spline_.set_points(inliers_x, inliers_y);
  // }else{
  //   spline_.set_points(new_radii, new_focal_lengths);
  // }
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
    updated_params_ = updated_params;
    // SetParams(updated_params);
  }
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
  std::cout << "new_radii size: " << new_radii.size() << std::endl;
  std::cout << "new_focal_lengths size: " << new_focal_lengths.size() << std::endl;
  
  
  // use ransac technique to fit the spline
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, new_radii.size() - 1);
  tk::spline best_spline;
  const int max_iterations = 40;
  const double threshold = 20.0;
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
    tk::spline s;
    s.set_points(sample_x, sample_y);
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
    updated_params_ = updated_params;
    // SetParams(updated_params);
  }
}

double Camera::EvalFocalLength(double radius) const {return spline_(radius);}
inline std::vector<double> Camera::GetRawRadii() const {return raw_radii_;}
inline void Camera::SetRawRadii(const std::vector<double>& raw_radii) const {raw_radii_ = raw_radii;}
inline std::vector<double> Camera::GetTheta() const {return theta_;}
inline void Camera::SetTheta(const std::vector<double>& theta) const {theta_ = theta;}
inline tk::spline Camera::GetSpline() const {return spline_;}


}  // namespace colmap

#endif  // COLMAP_SRC_BASE_CAMERA_H_
