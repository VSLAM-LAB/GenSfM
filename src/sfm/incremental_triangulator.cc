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

#include "sfm/incremental_triangulator.h"
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <map>
#include <tuple>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/math/interpolators/cubic_hermite.hpp>




#include "base/projection.h"
#include "estimators/triangulation.h"
#include "estimators/implicit_camera_pose.h"
#include "estimators/implicit_cost_matrix.h"
#include "estimators/implicit_pose_refinement.h"
#include "util/math.h"
#include "util/misc.h"
// include boost


namespace colmap {

bool IncrementalTriangulator::Options::Check() const {
  CHECK_OPTION_GE(max_transitivity, 0);
  CHECK_OPTION_GT(create_max_angle_error, 0);
  CHECK_OPTION_GT(continue_max_angle_error, 0);
  CHECK_OPTION_GT(merge_max_reproj_error, 0);
  CHECK_OPTION_GT(complete_max_reproj_error, 0);
  CHECK_OPTION_GE(complete_max_transitivity, 0);
  CHECK_OPTION_GT(re_max_angle_error, 0);
  CHECK_OPTION_GE(re_min_ratio, 0);
  CHECK_OPTION_LE(re_min_ratio, 1);
  CHECK_OPTION_GE(re_max_trials, 0);
  CHECK_OPTION_GT(min_angle, 0);
  return true;
}

IncrementalTriangulator::IncrementalTriangulator(
    const CorrespondenceGraph* correspondence_graph,
    Reconstruction* reconstruction)
    : correspondence_graph_(correspondence_graph),
      reconstruction_(reconstruction) {}

size_t IncrementalTriangulator::TriangulateImage(const Options& options,
                                                 const image_t image_id, bool initial, bool standard_triangulation) {
  CHECK(options.Check());

  size_t num_tris = 0;

  ClearCaches();
  // std::cout<<"Entering triangulator::tiangulateImage: "<<std::endl;
  // std::cout << "Initial in IncrementalTriangulator::Triangulateimage: " << initial << std::endl;

  const Image& image = reconstruction_->Image(image_id);
  if (!image.IsRegistered()) {
    return num_tris;
  }
  // std::cout<<"Image is registered"<<std::endl;

  const Camera& camera = reconstruction_->Camera(image.CameraId());
  CorrData ref_corr_data;
  ref_corr_data.image_id = image_id;
  ref_corr_data.image = &image;
  ref_corr_data.camera = &camera;

  // Container for correspondences from reference observation to other images.
  std::vector<CorrData> corrs_data;
  bool update_calibration = true;

  // Try to triangulate all image observations.
  for (point2D_t point2D_idx = 0; point2D_idx < image.NumPoints2D();
       ++point2D_idx) {
   
    if(camera.GetRawRadii().size() <20){
      update_calibration = true;
    }
    const size_t num_triangulated =
        Find(options, image_id, point2D_idx,
             static_cast<size_t>(options.max_transitivity), &corrs_data);
    if (corrs_data.empty()) {
      continue;
    }
    
    const Point2D& point2D = image.Point2D(point2D_idx);
    ref_corr_data.point2D_idx = point2D_idx;
    ref_corr_data.point2D = &point2D; 

    if (num_triangulated == 0) {
      corrs_data.push_back(ref_corr_data);
     
      num_tris += Create(options, corrs_data, initial, standard_triangulation, update_calibration);
    } else {
      // Continue correspondences to existing 3D points.
      num_tris += Continue(options, ref_corr_data, corrs_data);
      // Create points from correspondences that are not continued.
      corrs_data.push_back(ref_corr_data);
      num_tris += Create(options, corrs_data, initial, standard_triangulation, update_calibration);
    }
    update_calibration = false;
  }

  return num_tris;
}

size_t IncrementalTriangulator::CompleteImage(const Options& options,
                                              const image_t image_id, bool standard_triangulation) {
  CHECK(options.Check());

  size_t num_tris = 0;

  ClearCaches();

  const Image& image = reconstruction_->Image(image_id);
  if (!image.IsRegistered()) {
    return num_tris;
  }

  const Camera& camera = reconstruction_->Camera(image.CameraId());
  // if (HasCameraBogusParams(options, camera)) {
  //   return num_tris;
  // }
  if(camera.GetRawRadii().size() == 0){
    return num_tris;
  }

  // Setup estimation options.
  EstimateTriangulationOptions tri_options;
  tri_options.min_tri_angle = DegToRad(options.min_angle);
  tri_options.residual_type =
      TriangulationEstimator::ResidualType::REPROJECTION_ERROR;
  tri_options.ransac_options.max_error = options.complete_max_reproj_error;
  tri_options.ransac_options.confidence = 0.9999;
  tri_options.ransac_options.min_inlier_ratio = 0.02;
  tri_options.ransac_options.max_num_trials = 10000;

  // Correspondence data for reference observation in given image. We iterate
  // over all observations of the image and each observation once becomes
  // the reference correspondence.
  CorrData ref_corr_data;
  ref_corr_data.image_id = image_id;
  ref_corr_data.image = &image;
  ref_corr_data.camera = &camera;

  // Container for correspondences from reference observation to other images.
  std::vector<CorrData> corrs_data;

  for (point2D_t point2D_idx = 0; point2D_idx < image.NumPoints2D();
       ++point2D_idx) {
    const Point2D& point2D = image.Point2D(point2D_idx);
    if (point2D.HasPoint3D()) {
      // Complete existing track.
      num_tris += Complete(options, point2D.Point3DId());
      continue;
    }

    if (options.ignore_two_view_tracks &&
        correspondence_graph_->IsTwoViewObservation(image_id, point2D_idx)) {
      continue;
    }

    const size_t num_triangulated =
        Find(options, image_id, point2D_idx,
             static_cast<size_t>(options.max_transitivity), &corrs_data);
    if (num_triangulated || corrs_data.empty()) {
      continue;
    }

    ref_corr_data.point2D = &point2D;
    ref_corr_data.point2D_idx = point2D_idx;
    corrs_data.push_back(ref_corr_data);

    // Setup data for triangulation estimation.
    std::vector<TriangulationEstimator::PointData> point_data;
    point_data.resize(corrs_data.size());
    std::vector<TriangulationEstimator::PoseData> pose_data;
    pose_data.resize(corrs_data.size());

    std::vector<Camera> cameras_temp(corrs_data.size());
    bool is_calibrated = false;
    for (size_t i = 0; i < corrs_data.size(); ++i) {
      const CorrData& corr_data = corrs_data[i];
      point_data[i].point = corr_data.point2D->XY();
      point_data[i].point_normalized =
          corr_data.camera->ImageToWorld(point_data[i].point);
      pose_data[i].proj_matrix = corr_data.image->ProjectionMatrix();
      pose_data[i].proj_center = corr_data.image->ProjectionCenter();
      pose_data[i].camera = corr_data.camera;

      double radius = point_data[i].point_normalized.norm();
      if (!pose_data[i].camera->IsFullyCalibrated(point_data[i].point)) {
        point_data[i].is_radial = true;
        continue;
      }
      if (!standard_triangulation){
        point_data[i].is_radial = true;
        continue;
      }
      
      is_calibrated = true;

      double focal_length = 0;
      focal_length = corr_data.camera->EvalFocalLength(radius);
      point_data[i].focal_length = focal_length;
      point_data[i].point_normalized_standard = point_data[i].point_normalized / point_data[i].focal_length;
      Eigen::Matrix3d K;
      K << focal_length, 0, pose_data[i].camera->PrincipalPointX(),
          0, focal_length, pose_data[i].camera->PrincipalPointY(),
          0, 0, 1;

      point_data[i].point_normalized = point_data[i].point_normalized_standard;
      cameras_temp[i].SetModelId(SimplePinholeCameraModel::model_id);
      cameras_temp[i].SetParams({focal_length, pose_data[i].camera->PrincipalPointX(), pose_data[i].camera->PrincipalPointY()});
    }

    // Enforce exhaustive sampling for small track lengths.
    const size_t kExhaustiveSamplingThreshold = 15;
    if (point_data.size() <= kExhaustiveSamplingThreshold) {
      tri_options.ransac_options.min_num_trials =
          NChooseK(point_data.size(), 2);
    }

    // Estimate triangulation.
    Eigen::Vector3d xyz;
    std::vector<char> inlier_mask;
    // standard_triangulation=false;
    if (!EstimateTriangulation(tri_options, point_data, pose_data, &inlier_mask,
                               &xyz, false, standard_triangulation)) {
      continue;
    }

    // Add inliers to estimated track.
    int num_constraints = 0;
    Track track;
    track.Reserve(corrs_data.size());
    for (size_t i = 0; i < inlier_mask.size(); ++i) {
      if (inlier_mask[i]) {
        const CorrData& corr_data = corrs_data[i];
        track.AddElement(corr_data.image_id, corr_data.point2D_idx);
        num_constraints += (pose_data[i].camera->IsFullyCalibrated(corr_data.point2D->XY())) ? 2 : 1;
      }
    }

    if(num_constraints < 4) {
      // this is a underconstrained point, we do not add it to the reconstruction since it will get filtered later anyways
      continue;
    }

    num_tris += track.Length();

    const point3D_t point3D_id = reconstruction_->AddPoint3D(xyz, track);
    modified_point3D_ids_.insert(point3D_id);
  }

  return num_tris;
}

size_t IncrementalTriangulator::CompleteTracks(
    const Options& options, const std::unordered_set<point3D_t>& point3D_ids) {
  CHECK(options.Check());

  size_t num_completed = 0;

  ClearCaches();

  for (const point3D_t point3D_id : point3D_ids) {
    num_completed += Complete(options, point3D_id);
  }

  return num_completed;
}

size_t IncrementalTriangulator::CompleteAllTracks(const Options& options) {
  CHECK(options.Check());

  size_t num_completed = 0;

  ClearCaches();

  for (const point3D_t point3D_id : reconstruction_->Point3DIds()) {
    num_completed += Complete(options, point3D_id);
  }

  return num_completed;
}

size_t IncrementalTriangulator::MergeTracks(
    const Options& options, const std::unordered_set<point3D_t>& point3D_ids) {
  CHECK(options.Check());

  size_t num_merged = 0;

  ClearCaches();

  for (const point3D_t point3D_id : point3D_ids) {
    num_merged += Merge(options, point3D_id);
  }

  return num_merged;
}

size_t IncrementalTriangulator::MergeAllTracks(const Options& options) {
  CHECK(options.Check());

  size_t num_merged = 0;

  ClearCaches();

  for (const point3D_t point3D_id : reconstruction_->Point3DIds()) {
    num_merged += Merge(options, point3D_id);
  }

  return num_merged;
}

size_t IncrementalTriangulator::Retriangulate(const Options& options) {
  CHECK(options.Check());

  size_t num_tris = 0;
  size_t num_registrations = reconstruction_->NumRegImages();
  bool standard_triangulation = true;
  // bool standard_triangulation = reconstruction_->NumPoints3D() == 0;
  // check the number of registered images per camera
  std::map<camera_t, size_t> num_reg_images_per_camera;
  for (const auto& image_id : reconstruction_->RegImageIds()) {
    const Image& image = reconstruction_->Image(image_id);
    num_reg_images_per_camera[image.CameraId()] += 1;
  }
  for (const auto& pair : num_reg_images_per_camera) {
    if (pair.second < options.min_num_reg_images) {
      standard_triangulation = false;
      break;
    }
  }

  ClearCaches();

  std::cout << "ClearCaches done " << std::endl;

  Options re_options = options;
  re_options.continue_max_angle_error = options.re_max_angle_error;

  for (const auto& image_pair : reconstruction_->ImagePairs()) {
    // Only perform retriangulation for under-reconstructed image pairs.
    const double tri_ratio =
        static_cast<double>(image_pair.second.num_tri_corrs) /
        static_cast<double>(image_pair.second.num_total_corrs);
    if (tri_ratio >= options.re_min_ratio) {
      continue;
    }
    

    // Check if images are registered yet.

    image_t image_id1;
    image_t image_id2;
    Database::PairIdToImagePair(image_pair.first, &image_id1, &image_id2);

    const Image& image1 = reconstruction_->Image(image_id1);
    if (!image1.IsRegistered()) {
      continue;
    }

    const Image& image2 = reconstruction_->Image(image_id2);
    if (!image2.IsRegistered()) {
      continue;
    }

    // Only perform retriangulation for a maximum number of trials.

    int& num_re_trials = re_num_trials_[image_pair.first];
    if (num_re_trials >= options.re_max_trials) {
      continue;
    }
    num_re_trials += 1;

    const Camera& camera1 = reconstruction_->Camera(image1.CameraId());
    const Camera& camera2 = reconstruction_->Camera(image2.CameraId());

    // Find correspondences and perform retriangulation.
    // check the number of registered images for the camera of the first image
    size_t num_reg_camera1 = 0;
    size_t num_reg_camera2 = 0;
    for (const auto& image_id : reconstruction_->RegImageIds()) {
      if (reconstruction_->Image(image_id).CameraId() == image1.CameraId()) {
        num_reg_camera1 += 1;
      }
      if (reconstruction_->Image(image_id).CameraId() == image2.CameraId()) {
        num_reg_camera2 += 1;
      }
    }

    const FeatureMatches& corrs =
        correspondence_graph_->FindCorrespondencesBetweenImages(image_id1,
                                                                image_id2);

    // std::cout << "Correspondences found " << corrs.size() << std::endl;
    for (const auto& corr : corrs) {
      const Point2D& point2D1 = image1.Point2D(corr.point2D_idx1);
      const Point2D& point2D2 = image2.Point2D(corr.point2D_idx2);

      // Two cases are possible here: both points belong to the same 3D point
      // or to different 3D points. In the former case, there is nothing
      // to do. In the latter case, we do not attempt retriangulation,
      // as retriangulated correspondences are very likely bogus and
      // would therefore destroy both 3D points if merged.
      if (point2D1.HasPoint3D() && point2D2.HasPoint3D()) {
        continue;
      }

      CorrData corr_data1;
      corr_data1.image_id = image_id1;
      corr_data1.point2D_idx = corr.point2D_idx1;
      corr_data1.image = &image1;
      corr_data1.camera = &camera1;
      corr_data1.point2D = &point2D1;

      CorrData corr_data2;
      corr_data2.image_id = image_id2;
      corr_data2.point2D_idx = corr.point2D_idx2;
      corr_data2.image = &image2;
      corr_data2.camera = &camera2;
      corr_data2.point2D = &point2D2;

      if (point2D1.HasPoint3D() && !point2D2.HasPoint3D()) {
        const std::vector<CorrData> corrs_data1 = {corr_data1};
        num_tris += Continue(re_options, corr_data2, corrs_data1);
      } else if (!point2D1.HasPoint3D() && point2D2.HasPoint3D()) {
        const std::vector<CorrData> corrs_data2 = {corr_data2};
        num_tris += Continue(re_options, corr_data1, corrs_data2);
      } else if (!point2D1.HasPoint3D() && !point2D2.HasPoint3D()) {
        const std::vector<CorrData> corrs_data = {corr_data1, corr_data2};
        // Do not use larger triangulation threshold as this causes
        // significant drift when creating points (options vs. re_options).

        bool update_calibration = false;
        bool force_update = false;
        if(standard_triangulation){
          update_calibration = true;
          force_update = true;
        }
        num_tris += Create(re_options, corrs_data, false, false, false, false);
      }
      // Else both points have a 3D point, but we do not want to
      // merge points in retriangulation.
    }
  }

  return num_tris;
}

void IncrementalTriangulator::AddModifiedPoint3D(const point3D_t point3D_id) {
  modified_point3D_ids_.insert(point3D_id);
}

const std::unordered_set<point3D_t>&
IncrementalTriangulator::GetModifiedPoints3D() {
  // First remove any missing 3D points from the set.
  for (auto it = modified_point3D_ids_.begin();
       it != modified_point3D_ids_.end();) {
    if (reconstruction_->ExistsPoint3D(*it)) {
      ++it;
    } else {
      modified_point3D_ids_.erase(it++);
    }
  }
  return modified_point3D_ids_;
}

void IncrementalTriangulator::ClearModifiedPoints3D() {
  modified_point3D_ids_.clear();
}

void IncrementalTriangulator::ClearCaches() {
  camera_has_bogus_params_.clear();
  merge_trials_.clear();
}

size_t IncrementalTriangulator::Find(const Options& options,
                                     const image_t image_id,
                                     const point2D_t point2D_idx,
                                     const size_t transitivity,
                                     std::vector<CorrData>* corrs_data) {
  const std::vector<CorrespondenceGraph::Correspondence>& corrs =
      correspondence_graph_->FindTransitiveCorrespondences(
          image_id, point2D_idx, transitivity);

  corrs_data->clear();
  corrs_data->reserve(corrs.size());

  size_t num_triangulated = 0;

  for (const CorrespondenceGraph::Correspondence corr : corrs) {
    const Image& corr_image = reconstruction_->Image(corr.image_id);
    if (!corr_image.IsRegistered()) {
      continue;
    }

    const Camera& corr_camera = reconstruction_->Camera(corr_image.CameraId());
    // if (HasCameraBogusParams(options, corr_camera)) {
    //   continue;
    // }

    CorrData corr_data;
    corr_data.image_id = corr.image_id;
    corr_data.point2D_idx = corr.point2D_idx;
    corr_data.image = &corr_image;
    corr_data.camera = &corr_camera;
    corr_data.point2D = &corr_image.Point2D(corr.point2D_idx);

    corrs_data->push_back(corr_data);

    if (corr_data.point2D->HasPoint3D()) {
      num_triangulated += 1;
    }
  }

  return num_triangulated;
}

void sortRadiiAndFocalLengths(std::vector<double>& raw_radii, std::vector<double>& focal_lengths) {
    if (raw_radii.size() != focal_lengths.size()) {
        std::cerr << "Mismatched sizes: radii and focal lengths vectors must be the same size." << std::endl;
        return;
    }

    std::vector<std::pair<double, double>> radius_focal_pairs;
    for (size_t i = 0; i < raw_radii.size(); ++i) {
        radius_focal_pairs.push_back(std::make_pair(raw_radii[i], focal_lengths[i]));
    }

    std::sort(radius_focal_pairs.begin(), radius_focal_pairs.end(),
              [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                  return a.first < b.first;
              });

    for (size_t i = 0; i < radius_focal_pairs.size(); ++i) {
        raw_radii[i] = radius_focal_pairs[i].first;
        focal_lengths[i] = radius_focal_pairs[i].second;
    }
}

std::vector<double> IdentifyCalibratedArea(Camera& camera, std::vector<double>& radii, std::vector<double>& focal_lengths){
  // ensure sorted in ascending order
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
  threshold = std::max(threshold, DegToRad(0.1));
  camera.recursiveSplit(new_radii, new_focal_lengths, radii_segments, focal_lengths_segments, threshold, std_threshold);
  int longest_segment = 0;
  int longest_segment_size = 0;
  for(int i = 0; i < radii_segments.size(); i++){
    if(radii_segments[i].size() > longest_segment_size){
      longest_segment = i;
      longest_segment_size = radii_segments[i].size();
    }
  }
  if (radii_segments.size() == 0) {
    return {0, 0};
  }
  std::vector<double> calibrated_area = {radii_segments[longest_segment].front(), radii_segments[longest_segment].back()};
  return calibrated_area;
}

std::vector<double> movingAverage(const std::vector<double>& data, int window_size) {
    std::vector<double> result;
    int n = data.size();
    if (n == 0 || window_size <= 0) return result;

    double window_sum = 0;
    int window_start = 0, window_end = 0;

    // Initialize the sum with the first 'window_size' elements
    for (window_end = 0; window_end < window_size && window_end < n; ++window_end) {
        window_sum += data[window_end];
    }

    // Store the first average
    result.push_back(window_sum / window_size);

    // Slide the window; add a new element and remove the old element
    for (; window_end < n; ++window_end) {
        window_sum += data[window_end] - data[window_start++];
        result.push_back(window_sum / window_size);
    }

    return result;
}

int IncrementalTriangulator::CalibrateCamera(const Options& options) {
  // Only calibrate cameras if there are enough registered images.

  // check the number of registrated images per camera 
  std::unordered_map<camera_t, std::vector<image_t>> camera_images;

  int num_updated_cameras = 0;
  for (const auto& image_id_this : reconstruction_->RegImageIds()) {
    camera_t camera_id = reconstruction_->Image(image_id_this).CameraId();
    if (camera_images.find(camera_id) == camera_images.end()) {
      camera_images[camera_id] = std::vector<image_t>();
    }
    camera_images[camera_id].push_back(image_id_this);
  }

  PoseRefinementOptions pose_refinement_options;
  CostMatrixOptions cm_options;
  for (auto &[camera_id, camera_const]: reconstruction_->Cameras()) {
    if (camera_images[camera_id].size() < options.min_num_reg_images) {
      continue;
    }
    Camera& camera = reconstruction_->Camera(camera_id);

    // extracting points3D
    std::vector<std::vector<Eigen::Vector3d>> points3D;
    std::vector<std::vector<int>> points3D_lens;
    
    //populating camera poses
    std::vector<CameraPose> poses;
    for (const image_t image_id : camera_images[camera_id]) {
        Image& image = reconstruction_->Image(image_id);
        // image.NormalizeQvec();
        Eigen::Vector4d q(image.Qvec());
        Eigen::Vector3d t(image.Tvec());
        poses.push_back(CameraPose(q, t));
    }

    std::vector<std::vector<Eigen::Vector2d>> points2D;
    int total_track_length = 0;
    for (const image_t image_id_curr : camera_images[camera_id]) {
      const Image& image = reconstruction_->Image(image_id_curr);
      std::vector<Eigen::Vector2d> imagePoints2D;
      std::vector<Eigen::Vector3d> imagePoints3D;
      std::vector<int> imagePoints3D_lens;

      for (const Point2D& point2D : image.Points2D()) {
          if (!point2D.HasPoint3D()) {
              continue;
          }
          imagePoints2D.push_back(point2D.XY());
          imagePoints3D.push_back(reconstruction_->Point3D(point2D.Point3DId()).XYZ());
          imagePoints3D_lens.push_back(reconstruction_->Point3D(point2D.Point3DId()).Track().Length());
          total_track_length += imagePoints3D_lens.back();
      }

      points2D.push_back(imagePoints2D);
      points3D.push_back(imagePoints3D);
      points3D_lens.push_back(imagePoints3D_lens);
    }
    if(points2D.size() == 0){
      continue;
    }
    // randomly subsample the points, capping the maximum total number of imagepoints to 10000
    std::vector<std::vector<Eigen::Vector2d>> points2D_subsampled;
    std::vector<std::vector<Eigen::Vector3d>> points3D_subsampled;
    int total_num_points = 0;
    for (int i = 0; i < points2D.size(); i++) {
      total_num_points += points2D[i].size();
    }
    if (total_num_points == 0) {
      continue;
    }
    if(total_num_points > 10000){
      double ratio = 10000.0 / total_num_points;
      for (int i = 0; i < points2D.size(); i++) {
        std::vector<Eigen::Vector2d> imagePoints2D;
        std::vector<Eigen::Vector3d> imagePoints3D;
        for (int j = 0; j < points2D[i].size(); j++) {
          if (rand() % 100 < ratio * 100) {
            imagePoints2D.push_back(points2D[i][j]);
            imagePoints3D.push_back(points3D[i][j]);
          }
        }
        points2D_subsampled.push_back(imagePoints2D);
        points3D_subsampled.push_back(imagePoints3D);
      }
    } else {
      points2D_subsampled = points2D;
      points3D_subsampled = points3D;
    }
    
    // principal point from a random camera
    Eigen::Vector2d principal_point(reconstruction_->Camera(camera_id).PrincipalPointX(),
                                    reconstruction_->Camera(camera_id).PrincipalPointY());

    CostMatrix costMatrix = build_cost_matrix_multi(points2D_subsampled, cm_options, principal_point);
    IntrinsicCalib intrinsic_calib = calibrate_multi(points2D_subsampled, points3D_subsampled, costMatrix, principal_point, poses);
    std::cout << "!!! points2D.size(): " << points2D.size() << std::endl;
    std::vector<double> radii;
    std::vector<double> theta;
    std::vector<double> focal_lengths;
    for (const auto& pair : intrinsic_calib.r_f) {
        radii.push_back(pair.first); // Extract the first element of each pair (r)
        focal_lengths.push_back(pair.second); // Extract the second element of each pair (f)
    }
    for (const auto& pair : intrinsic_calib.theta_r) {
        theta.push_back(pair.first); // Extract the first element of each pair (r)
    }
    // identify the calibrated area
    std::vector<double> calibrated_area = IdentifyCalibratedArea(camera, theta, radii);
    std::vector<double> principal_point_new = {principal_point[0], principal_point[1]};
    // Check whether the calibrtion is wrong 
    double original_calibrated_area = camera.Params()[11] - camera.Params()[2];
    bool use_new_calibration = !camera.IsCalibrated();
    use_new_calibration = use_new_calibration
       || (original_calibrated_area < calibrated_area[1] - calibrated_area[0]);
    
    double diagonal = sqrt(camera.Width() * camera.Width() + camera.Height() * camera.Height()) / 2;

    int num_control_points = (camera.Params().size() - 2) / 2;
    for (size_t j = 2 + num_control_points; j < camera.Params().size(); j++) {
      // If the calibrated is manually forced, then use the new calibration
      if (camera.Params()[j] > diagonal || camera.Params()[j] < 0)
        use_new_calibration = true;
      if (j != camera.Params().size() - 1 && std::abs(camera.Params()[j + 1] - camera.Params()[j]) <= 1e-3) {
        use_new_calibration = true;
      }
    }
    // if(calibrated_area[1] - calibrated_area[0] < 0.1){
    //   use_new_calibration = false;
    // }
    if (use_new_calibration) {
      // If it is possible to calibrate
      if (!camera.FitPIeceWiseSpline_binary(theta, radii, principal_point_new))
        continue;
      camera.SetRawRadii(radii);
      camera.SetCalibrated(true);

      num_updated_cameras += 1;
    }
    std::cout << "Calibrated region:" << camera.Params()[2 + num_control_points] << " " << camera.Params()[1 + num_control_points * 2] << ", " << camera.Width() << " " << camera.Height() << std::endl;
  }

  return num_updated_cameras;
}

std::string currentDateTime() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}


size_t IncrementalTriangulator::Create(
    const Options& options, const std::vector<CorrData>& corrs_data, bool initial, bool standard_triangulation, bool update_calibration, bool false_update) {
  //  std::vector<image_t> img_ids;
  bool full_error = true;
  int debug_counter = 0;

  // Extract correspondences without an existing triangulated observation.
  std::vector<CorrData> create_corrs_data;
  create_corrs_data.reserve(corrs_data.size());
  for (const CorrData& corr_data : corrs_data) {
    if (!corr_data.point2D->HasPoint3D()) {
      create_corrs_data.push_back(corr_data);
    }
    
  }
  if (create_corrs_data.size() < 2) {
    // Need at least two observations for triangulation.
    return 0;
  } else if (options.ignore_two_view_tracks && create_corrs_data.size() == 2) {
    const CorrData& corr_data1 = create_corrs_data[0];
    if (correspondence_graph_->IsTwoViewObservation(corr_data1.image_id,
                                                    corr_data1.point2D_idx)) {
      return 0;
    }
  }
  // Setup data for triangulation estimation.
  std::vector<TriangulationEstimator::PointData> point_data;
  point_data.resize(create_corrs_data.size());
  std::vector<TriangulationEstimator::PoseData> pose_data;
  pose_data.resize(create_corrs_data.size());
  int count = 0;
  int radial_count = 0;
  std::vector<Camera> cameras_temp(create_corrs_data.size());
  bool is_calibrated = false;
  for (size_t i = 0; i < create_corrs_data.size(); ++i) {
    
    const CorrData& corr_data = create_corrs_data[i];
    point_data[i].point = corr_data.point2D->XY();
    point_data[i].point_normalized =
        corr_data.camera->ImageToWorld(point_data[i].point);
    pose_data[i].proj_matrix = corr_data.image->ProjectionMatrix();
    pose_data[i].proj_center = corr_data.image->ProjectionCenter();
    
    pose_data[i].camera = corr_data.camera;
    
    pose_data[i].image = corr_data.image;
    double radius = point_data[i].point_normalized.norm();
    bool is_calibrated = corr_data.camera->IsFullyCalibrated(point_data[i].point);
    if (!is_calibrated)
      continue;
    
    double focal_length_splined = corr_data.camera->EvalFocalLength(radius);
    point_data[i].focal_length = focal_length_splined;
  
    
    point_data[i].point_normalized_standard = point_data[i].point_normalized / point_data[i].focal_length;
    Eigen::Matrix3d K;
    K << focal_length_splined, 0, pose_data[i].camera->PrincipalPointX(),
         0, focal_length_splined, pose_data[i].camera->PrincipalPointY(),
         0, 0, 1;
    pose_data[i].proj_matrix_standard = K * pose_data[i].proj_matrix;

    // Construct temporary camera for triangulation.
    point_data[i].point_normalized = point_data[i].point_normalized_standard;
    cameras_temp[i].SetModelId(SimplePinholeCameraModel::model_id);
    cameras_temp[i].SetParams({focal_length_splined, pose_data[i].camera->PrincipalPointX(), pose_data[i].camera->PrincipalPointY()});

    pose_data[i].camera = &cameras_temp[i];

    

  }
  if (count > 0) {
    full_error = false;
  }
  if(radial_count >= 0.5*create_corrs_data.size()){
    standard_triangulation = false;
  }

  EstimateTriangulationOptions tri_options;
  tri_options.min_tri_angle = DegToRad(options.min_angle);
  tri_options.residual_type =
      TriangulationEstimator::ResidualType::ANGULAR_ERROR;
  tri_options.ransac_options.max_error =
      DegToRad(options.create_max_angle_error);
  tri_options.ransac_options.confidence = 0.9999;
  tri_options.ransac_options.min_inlier_ratio = 0.02;
  tri_options.ransac_options.max_num_trials = 10000;

  // Enforce exhaustive sampling for small track lengths.
  const size_t kExhaustiveSamplingThreshold = 15;
  if (point_data.size() <= kExhaustiveSamplingThreshold) {
    tri_options.ransac_options.min_num_trials = NChooseK(point_data.size(), 3);
  }

  // Estimate triangulation.
  Eigen::Vector3d xyz;
  std::vector<char> inlier_mask;
  if (!EstimateTriangulation(tri_options, point_data, pose_data, &inlier_mask,
                             &xyz, initial, false)) {
    return 0;
  }
  // Add inliers to estimated track.
  int num_constraints = 0;
  Track track;
  track.Reserve(create_corrs_data.size());
  for (size_t i = 0; i < inlier_mask.size(); ++i) {
    if (inlier_mask[i]) {
      const CorrData& corr_data = create_corrs_data[i];
      track.AddElement(corr_data.image_id, corr_data.point2D_idx);
      num_constraints += (pose_data[i].camera->IsFullyCalibrated(corr_data.point2D->XY())) ? 2 : 1;
    }
  }
  if(num_constraints < 4) {
    // this is a underconstrained point, we do not add it to the reconstruction since it will get filtered later anyways
    return 0;
  }
  point3D_t point3D_id;
  point3D_id = reconstruction_->AddPoint3D(xyz, track);
  modified_point3D_ids_.insert(point3D_id);

  const size_t kMinRecursiveTrackLength = 3;
  // if(!used_full){
  if (create_corrs_data.size() - track.Length() >= kMinRecursiveTrackLength) {
    update_calibration = false;
    return track.Length() + Create(options, create_corrs_data, initial, standard_triangulation, update_calibration);
  }

  return track.Length();
}

size_t IncrementalTriangulator::Continue(
    const Options& options, const CorrData& ref_corr_data,
    const std::vector<CorrData>& corrs_data) {
  // No need to continue, if the reference observation is triangulated.
  if (ref_corr_data.point2D->HasPoint3D()) {
    return 0;
  }

  double best_angle_error = std::numeric_limits<double>::max();
  size_t best_idx = std::numeric_limits<size_t>::max();

  for (size_t idx = 0; idx < corrs_data.size(); ++idx) {
    const CorrData& corr_data = corrs_data[idx];
    if (!corr_data.point2D->HasPoint3D()) {
      continue;
    }

    const Point3D& point3D =
        reconstruction_->Point3D(corr_data.point2D->Point3DId());

    const double angle_error = CalculateAngularError(
        ref_corr_data.point2D->XY(), point3D.XYZ(), ref_corr_data.image->Qvec(),
        ref_corr_data.image->Tvec(), *ref_corr_data.camera);
    if (angle_error < best_angle_error) {
      best_angle_error = angle_error;
      best_idx = idx;
    }
  }

  const double max_angle_error = DegToRad(options.continue_max_angle_error);
  if (best_angle_error <= max_angle_error &&
      best_idx != std::numeric_limits<size_t>::max()) {
    const CorrData& corr_data = corrs_data[best_idx];
    const TrackElement track_el(ref_corr_data.image_id,
                                ref_corr_data.point2D_idx);
    reconstruction_->AddObservation(corr_data.point2D->Point3DId(), track_el);
    modified_point3D_ids_.insert(corr_data.point2D->Point3DId());
    return 1;
  }

  return 0;
}

size_t IncrementalTriangulator::Merge(const Options& options,
                                      const point3D_t point3D_id) {
  if (!reconstruction_->ExistsPoint3D(point3D_id)) {
    return 0;
  }

  const double max_squared_reproj_error =
      options.merge_max_reproj_error * options.merge_max_reproj_error;

  const auto& point3D = reconstruction_->Point3D(point3D_id);

  for (const auto& track_el : point3D.Track().Elements()) {
    const std::vector<CorrespondenceGraph::Correspondence>& corrs =
        correspondence_graph_->FindCorrespondences(track_el.image_id,
                                                   track_el.point2D_idx);

    for (const auto corr : corrs) {
      const auto& image = reconstruction_->Image(corr.image_id);
      if (!image.IsRegistered()) {
        continue;
      }

      const Point2D& corr_point2D = image.Point2D(corr.point2D_idx);
      if (!corr_point2D.HasPoint3D() ||
          corr_point2D.Point3DId() == point3D_id ||
          merge_trials_[point3D_id].count(corr_point2D.Point3DId()) > 0) {
        continue;
      }

      // Try to merge the two 3D points.

      const Point3D& corr_point3D =
          reconstruction_->Point3D(corr_point2D.Point3DId());

      merge_trials_[point3D_id].insert(corr_point2D.Point3DId());
      merge_trials_[corr_point2D.Point3DId()].insert(point3D_id);

      // Weighted average of point locations, depending on track length.
      const Eigen::Vector3d merged_xyz =
          (point3D.Track().Length() * point3D.XYZ() +
           corr_point3D.Track().Length() * corr_point3D.XYZ()) /
          (point3D.Track().Length() + corr_point3D.Track().Length());

      // Count number of inlier track elements of the merged track.
      bool merge_success = true;
      for (const Track* track : {&point3D.Track(), &corr_point3D.Track()}) {
        for (const auto test_track_el : track->Elements()) {
          const Image& test_image =
              reconstruction_->Image(test_track_el.image_id);
          const Camera& test_camera =
              reconstruction_->Camera(test_image.CameraId());
          const Point2D& test_point2D =
              test_image.Point2D(test_track_el.point2D_idx);
          if (CalculateSquaredReprojectionError(
                  test_point2D.XY(), merged_xyz, test_image.Qvec(),
                  test_image.Tvec(), test_camera) > max_squared_reproj_error) {
            merge_success = false;
            break;
          }
        }
        if (!merge_success) {
          break;
        }
      }

      // Only accept merge if all track elements are inliers.
      if (merge_success) {
        const size_t num_merged =
            point3D.Track().Length() + corr_point3D.Track().Length();

        const point3D_t merged_point3D_id = reconstruction_->MergePoints3D(
            point3D_id, corr_point2D.Point3DId());

        modified_point3D_ids_.erase(point3D_id);
        modified_point3D_ids_.erase(corr_point2D.Point3DId());
        modified_point3D_ids_.insert(merged_point3D_id);

        // Merge merged 3D point and return, as the original points are deleted.
        const size_t num_merged_recursive = Merge(options, merged_point3D_id);
        if (num_merged_recursive > 0) {
          return num_merged_recursive;
        } else {
          return num_merged;
        }
      }
    }
  }

  return 0;
}

size_t IncrementalTriangulator::Complete(const Options& options,
                                         const point3D_t point3D_id) {
  size_t num_completed = 0;

  if (!reconstruction_->ExistsPoint3D(point3D_id)) {
    return num_completed;
  }

  const double max_squared_reproj_error =
      options.complete_max_reproj_error * options.complete_max_reproj_error;

  const Point3D& point3D = reconstruction_->Point3D(point3D_id);

  std::vector<TrackElement> queue = point3D.Track().Elements();

  const int max_transitivity = options.complete_max_transitivity;
  for (int transitivity = 0; transitivity < max_transitivity; ++transitivity) {
    if (queue.empty()) {
      break;
    }

    const auto prev_queue = queue;
    queue.clear();

    for (const TrackElement queue_elem : prev_queue) {
      const std::vector<CorrespondenceGraph::Correspondence>& corrs =
          correspondence_graph_->FindCorrespondences(queue_elem.image_id,
                                                     queue_elem.point2D_idx);

      for (const auto corr : corrs) {
        const Image& image = reconstruction_->Image(corr.image_id);
        if (!image.IsRegistered()) {
          continue;
        }

        const Point2D& point2D = image.Point2D(corr.point2D_idx);
        if (point2D.HasPoint3D()) {
          continue;
        }

        const Camera& camera = reconstruction_->Camera(image.CameraId());
        // if (HasCameraBogusParams(options, camera)) {
        //   continue;
        // }

        if (CalculateSquaredReprojectionError(
                point2D.XY(), point3D.XYZ(), image.Qvec(), image.Tvec(),
                camera) > max_squared_reproj_error) {
          continue;
        }
        const TrackElement track_el(corr.image_id, corr.point2D_idx);
        reconstruction_->AddObservation(point3D_id, track_el);
        modified_point3D_ids_.insert(point3D_id);

        // Recursively complete track for this new correspondence.
        if (transitivity < max_transitivity - 1) {
          queue.emplace_back(corr.image_id, corr.point2D_idx);
        }

        num_completed += 1;
      }
    }
  }

  return num_completed;
}

bool IncrementalTriangulator::HasCameraBogusParams(const Options& options,
                                                   const Camera& camera) {
  const auto it = camera_has_bogus_params_.find(camera.CameraId());
  if (it == camera_has_bogus_params_.end()) {
    const bool has_bogus_params = camera.HasBogusParams(
        options.min_focal_length_ratio, options.max_focal_length_ratio,
        options.max_extra_param);
    camera_has_bogus_params_.emplace(camera.CameraId(), has_bogus_params);
    return has_bogus_params;
  } else {
    return it->second;
  }
}

}  // namespace colmap
