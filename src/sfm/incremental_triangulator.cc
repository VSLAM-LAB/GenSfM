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
  // if (HasCameraBogusParams(options, camera)) {
  //   return num_tris;
  // }
  // std::cout<<"Camera has no bogus params"<<std::endl;

  // Correspondence data for reference observation in given image. We iterate
  // over all observations of the image and each observation once becomes
  // the reference correspondence.
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
    // std::cout<<"NumPoints2D: "<<image.NumPoints2D()<<std::endl;
    // if(point2D_idx%3000 == 0){
    //   update_calibration = true;
    // }
    if(camera.GetRawRadii().size() <20){
      update_calibration = true;
    }
    const size_t num_triangulated =
        Find(options, image_id, point2D_idx,
             static_cast<size_t>(options.max_transitivity), &corrs_data);
    if (corrs_data.empty()) {
      continue;
    }
    // std::cout<<"Corrs data is not empty"<<std::endl;

    const Point2D& point2D = image.Point2D(point2D_idx);
    ref_corr_data.point2D_idx = point2D_idx;
    ref_corr_data.point2D = &point2D;
    // std::cout<<"num_triangulated: "<<num_triangulated<<std::endl; 

    if (num_triangulated == 0) {
      corrs_data.push_back(ref_corr_data);
      // if(standard_triangulation){
      //   num_tris += Create(options, corrs_data, initial, false);
      // } 
      num_tris += Create(options, corrs_data, initial, standard_triangulation, update_calibration);
    } else {
      // Continue correspondences to existing 3D points.
      num_tris += Continue(options, ref_corr_data, corrs_data);
      // Create points from correspondences that are not continued.
      corrs_data.push_back(ref_corr_data);
      // if(standard_triangulation){
      //   num_tris += Create(options, corrs_data, initial, false);
      // }
      num_tris += Create(options, corrs_data, initial, standard_triangulation, update_calibration);
    }
    update_calibration = false;
  }

  return num_tris;
}


// isolating the initial triangulation
// size_t IncrementalTriangulator::TriangulateImageInitial(const Options& options,
//                                                  const image_t image_id) {
//   CHECK(options.Check());

//   size_t num_tris = 0;

//   ClearCaches();

//   const Image& image = reconstruction_->Image(image_id);
//   if (!image.IsRegistered()) {
//     return num_tris;
//   }

//   const Camera& camera = reconstruction_->Camera(image.CameraId());
//   if (HasCameraBogusParams(options, camera)) {
//     return num_tris;
//   }

//   // Correspondence data for reference observation in given image. We iterate
//   // over all observations of the image and each observation once becomes
//   // the reference correspondence.
//   CorrData ref_corr_data;
//   ref_corr_data.image_id = image_id;
//   ref_corr_data.image = &image;
//   ref_corr_data.camera = &camera;

//   // Container for correspondences from reference observation to other images.
//   std::vector<CorrData> corrs_data;

//   // Try to triangulate all image observations.
//   for (point2D_t point2D_idx = 0; point2D_idx < image.NumPoints2D();
//        ++point2D_idx) {
//     const size_t num_triangulated =
//         Find(options, image_id, point2D_idx,
//              static_cast<size_t>(options.max_transitivity), &corrs_data);
//     if (corrs_data.empty()) {
//       continue;
//     }

//     const Point2D& point2D = image.Point2D(point2D_idx);
//     ref_corr_data.point2D_idx = point2D_idx;
//     ref_corr_data.point2D = &point2D;

//     if (num_triangulated == 0) {
//       corrs_data.push_back(ref_corr_data);
//       num_tris += CreateInitial(options, corrs_data);
//     } else {
//       // Continue correspondences to existing 3D points.
//       num_tris += Continue(options, ref_corr_data, corrs_data);
//       // Create points from correspondences that are not continued.
//       corrs_data.push_back(ref_corr_data);
//       num_tris += CreateInitial(options, corrs_data);
//     }
//   }

//   return num_tris;
// }



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
    for (size_t i = 0; i < corrs_data.size(); ++i) {
      const CorrData& corr_data = corrs_data[i];
      point_data[i].point = corr_data.point2D->XY();
      point_data[i].point_normalized =
          corr_data.camera->ImageToWorld(point_data[i].point);
      pose_data[i].proj_matrix = corr_data.image->ProjectionMatrix();
      pose_data[i].proj_center = corr_data.image->ProjectionCenter();
      pose_data[i].camera = corr_data.camera;
      if(corr_data.camera->GetRawRadii().size() == 0){
        standard_triangulation = false;
      }
      if(standard_triangulation){
        std::vector<double> raw_radii = corr_data.camera->GetRawRadii();
        std::vector<double> focal_lengths = corr_data.camera->GetFocalLengthParams();
        // std::cout<<"raw_radii size: "<<raw_radii.size()<<std::endl;
        double radius = point_data[i].point_normalized.norm();
        double focal_length = 0;
        focal_length = corr_data.camera->EvalFocalLength(radius);
        //     if(raw_radii.size() == 0){
        //       continue;}
        // for (int i = 0; i < raw_radii.size() - 1; i++) {
        //   if (radius >= raw_radii[i] && radius <= raw_radii[i+1]) {
        //     // interpolate the focal length
        //     focal_length = focal_lengths[i] + (focal_lengths[i+1] - focal_lengths[i]) * (radius - raw_radii[i]) / (raw_radii[i+1] - raw_radii[i]+1e-6);
        //     break;
        //   }
        // }
        // // if the radius is smaller than the smallest radius, set the focal length to the smallest focal length
        // if (radius < raw_radii[0]) {
        //   focal_length = focal_lengths[0];
        // }
        // // if the radius is larger than the largest radius, set the focal length to the largest focal length
        // if (radius > raw_radii[raw_radii.size() - 1]) {
        //   focal_length = focal_lengths[raw_radii.size() - 1];
        // }
        point_data[i].focal_length = focal_length;
        point_data[i].point_normalized_standard = point_data[i].point_normalized / point_data[i].focal_length;
        Eigen::Matrix3d K;
        K << focal_length, 0, pose_data[i].camera->PrincipalPointX(),
            0, focal_length, pose_data[i].camera->PrincipalPointY(),
            0, 0, 1;
        pose_data[i].proj_matrix_standard = K * pose_data[i].proj_matrix;

      }
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
        if (standard_triangulation){
          num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 2 : 2;
        }else{
          num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 1 : 2;
        }
        // num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id) ? 1 : 2;
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
  bool standard_triangulation = false;
  // bool standard_triangulation = reconstruction_->NumPoints3D() == 0;

  ClearCaches();

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
    // if (HasCameraBogusParams(options, camera1) ||
    //     HasCameraBogusParams(options, camera2)) {
    //   continue;
    // }

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
    // related to min_num_reg_images
    if (num_reg_camera1 >=20 && num_reg_camera2 >= 20) {
      standard_triangulation = true;
    } 

    const FeatureMatches& corrs =
        correspondence_graph_->FindCorrespondencesBetweenImages(image_id1,
                                                                image_id2);

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
        // std::cout << "=============================Retriangulating===========================" << std::endl;
        num_tris += Create(options, corrs_data, false, standard_triangulation, false);
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

    // Step 1: Create a vector of pairs
    std::vector<std::pair<double, double>> radius_focal_pairs;
    for (size_t i = 0; i < raw_radii.size(); ++i) {
        radius_focal_pairs.push_back(std::make_pair(raw_radii[i], focal_lengths[i]));
    }

    // Step 2: Sort the vector of pairs based on the first element of each pair (the radius)
    std::sort(radius_focal_pairs.begin(), radius_focal_pairs.end(),
              [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
                  return a.first < b.first;
              });

    // Step 3: Unpack the pairs back into the original vectors
    for (size_t i = 0; i < radius_focal_pairs.size(); ++i) {
        raw_radii[i] = radius_focal_pairs[i].first;
        focal_lengths[i] = radius_focal_pairs[i].second;
    }
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

std::string currentDateTime() {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}


size_t IncrementalTriangulator::Create(
    const Options& options, const std::vector<CorrData>& corrs_data, bool initial, bool standard_triangulation, bool update_calibration) {
  //  std::vector<image_t> img_ids;
  bool full_error = true;
  int debug_counter = 0;
  // std::cout<<"Entering Create"<<std::endl;
  if(standard_triangulation){
  // std::cout<<"Standard Triangulation"<<std::endl;
  
  // Construct raw_radii and focal lengths for the images in corrs_data
  // std::cout << "================== Starting Standard Triangulation =======================" << std::endl;
  // create a vector for unique image_ids in corrs_data
  std::map<camera_t,std::vector<std::tuple<CameraPose,std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector3d>, Eigen::Vector2d>>> camera_correspondences;
  // std::cout<<"Entering if condition"<<std::endl;
  
  for (const CorrData& corr_data : corrs_data) {
  
    // std::cout<<"Entering CorrData "<<std::endl;
    const Image& image = *corr_data.image;
    const Camera& camera = *corr_data.camera;
    camera_t camera_id = camera.CameraId(); 
    // check if the image has already been set a radii and focal_lengths
    // if (std::find(img_ids.begin(), img_ids.end(), image.ImageId()) != img_ids.end() ) {
    //   continue;
    // }
    // if(image.GetRawRadii().size()>0){
    
    // Retrieve the points2D in this image
    std::vector<class Point2D> points2D_exist = image.Points2D();
    std::vector<class Point2D> points2D;
    std::vector<Eigen::Vector2d> points2D_xy;
    std::vector<class Point3D> points3D;
    std::vector<Eigen::Vector3d> points3D_xyz;
    for (const class Point2D& point2D : points2D_exist) {
      if (point2D.HasPoint3D()) {
        points2D.push_back(point2D);
        points2D_xy.push_back(point2D.XY());
        points3D.push_back(reconstruction_->Point3D(point2D.Point3DId()));
        points3D_xyz.push_back(Eigen::Vector3d(points3D.back().X(), points3D.back().Y(), points3D.back().Z()));
      }
    }
    // std::cout<<"Points2D size: "<<points2D.size()<<std::endl;
    // std::cout<<"Points3D size: "<<points3D.size()<<std::endl;

    CameraPose pose(image.Qvec(), image.Tvec());
    // update the camera_correspondences
    

    // filter the outliers
    if(update_calibration){
      // std::cout<<"Updating calibration"<<std::endl;
      // if(true) {
      // std::cout<<"Updating calibration"<<std::endl;
      //
    // if(true) {
      // std::cout<<"Updating calibration"<<std::endl;
      //
      PoseRefinementOptions pose_refinement_options;
      CostMatrixOptions cm_options;
      // std::cout<<"Building Cost matrix"<<std::endl;
      CostMatrix cost_matrix = build_cost_matrix(points2D_xy, cm_options, Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY()));
      // std::cout<<"Cost matrix built"<<std::endl;
      Eigen::Vector2d pp = Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY());
      // Eigen::Vector2d pp = Eigen::Vector2d(3041.29, 2014.07);
      // refining pp
      CameraPose pose_refining = pose;
      Eigen::Vector2d pp_refining = pp;
      // pose_refining.t[2] = 0.;
      PoseRefinement1DRadialOptions opt;
      // min_num_reg_images related


      // check the number of registrated images per camera 
      std::pair<camera_t, std::vector<image_t>> camera_images;
      for (const auto& image_id_this : reconstruction_->RegImageIds()) {
        if (reconstruction_->Image(image_id_this).CameraId() == camera_id) {
          camera_images.first = camera_id;
          camera_images.second.push_back(image_id_this);
        }
      }
      // std::cout<<"Camera images size: "<<camera_images.second.size()<<std::endl;

      



      // if( reconstruction_->RegImageIds().size() <= 0){
      if( camera_images.second.size() <= 20){
        // std::cout<<"Pose refinement started"<<std::endl;
        pose_refining.t[2] = 0.;
    
      pose_refinement_1D_radial(points2D_xy, points3D_xyz, &pose_refining, &pp_refining,opt);
      }
    // Commented out for point_triangulator 
    // try to refine the pose
      CameraPose pose_refined;
      // try
      // {
      // pose_refined = pose_refinement(points2D_xy, points3D_xyz, cost_matrix, pp,pose, pose_refinement_options);
      // // std::cout<<"Pose refined"<<std::endl;
      // }
      // catch(const std::exception& e)
      // {
      //   std::cerr << e.what() << '\n';
      //   pose_refined = pose;
      // }
      if(points3D_xyz.size() < 20){
        pose_refined = pose;}
      // check the number of registrations
      
      else{
        pose_refined = pose_refinement(points2D_xy, points3D_xyz, cost_matrix, pp_refining,pose, pose_refinement_options);
        // image.SetPose(pose_refined);
        }
      // set the camera pose as the pose_refined
      // image.SetPose(pose_refined);
      // min_num_reg_images related
      // if( reconstruction_->RegImageIds().size() >0){
      if( camera_images.second.size() > 20){
        pose_refined = pose;
      }
      image_t image_id = image.ImageId();
      Image& image_update = reconstruction_->Image(image_id);
      image_update.SetQvec(pose_refined.q_vec);
      image_update.SetTvec(pose_refined.t);
      Camera& camera_update = reconstruction_->Camera(camera_id);
      camera_update.SetPrincipalPointX(pp_refining.x());
      camera_update.SetPrincipalPointY(pp_refining.y()); 
      
      
      // CameraPose pose_refined = pose_refinement(points2D_xy, points3D_xyz, cost_matrix, pp,pose, pose_refinement_options);
      // pose_refined = pose;
      filter_result_pose_refinement(points2D_xy, points3D_xyz, pose_refined, pp_refining, pose_refinement_options);
      // std::cout<<"Pose refined and filtered"<<std::endl;
      std::vector<Eigen::Vector2d> points2D_xy_filtered = points2D_xy;
      std::vector<Eigen::Vector3d> points3D_xyz_filtered = points3D_xyz;
      if(camera_correspondences.find(camera_id) == camera_correspondences.end()){
        camera_correspondences[camera_id] = std::vector<std::tuple<CameraPose,std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector3d>, Eigen::Vector2d>>();
        camera_correspondences[camera_id].push_back(std::make_tuple(pose_refined, points2D_xy, points3D_xyz, pp_refining));
      }else{
        camera_correspondences[camera_id].push_back(std::make_tuple(pose_refined, points2D_xy, points3D_xyz, pp_refining));
      }
    }
    // recover the intrinsic calibration parameters
    
    // CostMatrix cost_matrix_2;
    // cost_matrix_2=build_cost_matrix(points2D_xy_filtered, cm_options, pp);
    // std::cout<<"calibration started"<<std::endl;
    // IntrinsicCalib intrinsic_calib = calibrate(points2D_xy_filtered, points3D_xyz_filtered, cost_matrix_2, pp, pose_refined);
    // std::cout<<"calibration finished"<<std::endl;
    // Retrieve the distance from the principal point to the point2D
    // std::vector<double> radii;
    // std::vector<double> theta;
    // std::vector<double> gt_focal_lengths;
    // std::vector<double> radii_for_theta;
    // for (const auto& pair : intrinsic_calib.r_f) {
    //     radii.push_back(pair.first); // Extract the first element of each pair (r)
    // }
    // Eigen::Matrix3d R = image.RotationMatrix();
    // Eigen::Vector3d t = image.Tvec();
    // for (int i = 0; i < points2D_xy.size(); i++) {
    //   radii.push_back((points2D_xy[i] - Eigen::Vector2d(camera.PrincipalPointX(),camera.PrincipalPointY())).norm());
    // }
    // if(radii.size()>=40){
    // // image.SetRawRadii(radii);
    // camera.SetRawRadii(radii);
    // std::cout<<"radii set"<<std::endl;
    // Retrieve rotation matrix 
   
    // std::vector<Eigen::Vector3d> X_cam;

    // for (int i = 0; i < points3D_xyz.size(); i++) {
    //   Eigen::Vector3d Z = points3D_xyz[i];
    //   Eigen::Vector3d Xcam = R * Z + t;
    //   X_cam.push_back(Xcam);
    // }
    // Retrieve the focal length
    // std::vector<double> focal_lengths;
    // for (const auto& pair : intrinsic_calib.r_f) {
    //     focal_lengths.push_back(pair.second); // Extract the first element of each pair (r)
    // }
    // for (const auto& pair : intrinsic_calib.theta_r) {
    //     theta.push_back(pair.first); // Extract the first element of each pair (r)
    //     // radii_for_theta.push_back(pair.second);
    // }
    // image.SetTheta(theta);
    // camera.SetTheta(theta);
    // for (int i = 0; i < radii.size(); i++) {
    // Eigen::Vector2d offset = points2D_xy[i] - Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY());
    // double dot_product = offset.x() * X_cam[i].x() + offset.y() * X_cam[i].y(); // Calculate dot product of the projection

    // // Make sure to handle division by zero or very small dot_product values to avoid instability
    // if (std::abs(dot_product) > std::numeric_limits<double>::epsilon()) {
    //     double focal_length = radii[i] * radii[i] * X_cam[i].z() / dot_product;
    //     // std::cout << "Calculated Focal length: " << focal_length <<"for radii: " << radii[i] << std::endl;
    //     focal_lengths.push_back(focal_length);
    // } else {
    //     focal_lengths.push_back(1); // or some form of error value or handling
    //   }
    // }
    // double fx = 3411.42;
    // double fy = 3410.02;
    // double cx = 3041.29;
    // double cy = 2014.07;
    // double k1 = 0.21047;
    // double k2 = 0.21102;
    // double p1 = -5.36231e-06;
    // double p2 =  0.00051541;
    // double k3 = -0.158023;
    // double k4 =  0.406856;
    // double s1 =  -8.46499e-05;
    // double s2 =  0.000861313;
    // std::vector<double> r_3d;
    // // from 0 to 1, increment by 0.01
    // for (double i = 0; i < 1; i += 0.002) {
    //   r_3d.push_back(i);
    // }
    // std::vector<double> radii_list;

    // for (int i = 0; i < r_3d.size(); i++){
    //   double theta = atan(r_3d[i]);
    //   double x_d = r_3d[i] * cos(theta);
    //   double y_d = r_3d[i] * sin(theta);
    //   double u_d = theta/r_3d[i] * x_d;
    //   double v_d = theta/r_3d[i] * y_d;
    //   double theta_squared = theta * theta;
    //   double t_r = 1 + k1 * theta_squared + k2 * theta_squared * theta_squared + k3 * theta_squared * theta_squared * theta_squared + k4 * theta_squared * theta_squared * theta_squared * theta_squared;
    //   double u_n = u_d * t_r + 2 * p1 * u_d * v_d + p2 * (theta_squared + 2 * u_d * u_d) + s1 * theta_squared;
    //   double v_n = v_d * t_r + 2 * p2 * u_d * v_d + p1 * (theta_squared + 2 * v_d * v_d) + s2 * theta_squared;
    //   double u = fx * u_n + cx;
    //   double v = fy * v_n + cy;
    //   double radius = sqrt((u - cx) * (u - cx) + (v - cy) * (v - cy));
    //   radii_list.push_back(radius);
    //   double gt_f = (u-cx)/x_d;
    //   gt_focal_lengths.push_back(gt_f);
    // }
    // camera.SetFocalLengthParams(gt_focal_lengths);
    // camera.SetRawRadii(radii_list);


    // image.SetFocalLengthParams(focal_lengths);
    // camera.SetFocalLengthParams(focal_lengths);
    
    
    // std::cout<<"focal lengths set"<<std::endl;
    
    // img_ids.push_back(image.ImageId()); 
    // std::cout<<"img_ids size:"<<img_ids.size()<<std::endl;
    // std::cout<<"theta set"<<std::endl;
    // std::cout << "Radii: " << radii.size() << " Focal Lengths: " << focal_lengths.size() << std::endl;
    
    // -------------% Fitting a parametric spline model for estimated r-f map %------------ //
    // image.FitSpline(radii, focal_lengths);
    // camera.FitSpline(radii, focal_lengths);
    // camera.FitSpline(radii_list, gt_focal_lengths);
  
    
  }

    // if (radii.size() <=40) {
    //     standard_triangulation = false;
    // }
  // }
  
  
  // iterate through the camera_correspondences and perform intrinsic_calib
  // std::cout<<"camera_correspondences size: "<<camera_correspondences.size()<<std::endl;
  // std::cout<<"=====> Entering calibration loop: "<<std::endl;
  debug_counter++;
  // std::cout << "debug_counter: " << debug_counter << std::endl;
  if(update_calibration){
  // if(true){ 
  // std::cout<<"Update_calibration, Calibrating cameras"<<std::endl;
  for (const auto& pair : camera_correspondences) {
    camera_t camera_id = pair.first;
    Camera &camera_this = reconstruction_->Camera(camera_id);
    // if(camera_this.GetRawRadii().size()>30){
    //   continue;
    // }
    std::vector<std::tuple<CameraPose,std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector3d>, Eigen::Vector2d> > correspondences = pair.second;
    std::vector<std::vector<Eigen::Vector2d>> points2D_this;
    std::vector<std::vector<Eigen::Vector3d>> points3D_this;
    std::vector<CameraPose> poses;
    Eigen::Vector2d pp_this;
    // std::cout << "correspondences size: " << correspondences.size() <<"for camera: "<< camera_id << std::endl;
    for (int i = 0; i < correspondences.size(); i++) {
      CameraPose pose = std::get<0>(correspondences[i]);
      std::vector<Eigen::Vector2d> points2D_i = std::get<1>(correspondences[i]);
      std::vector<Eigen::Vector3d> points3D_i = std::get<2>(correspondences[i]);
      points2D_this.push_back(points2D_i);
      points3D_this.push_back(points3D_i);
      poses.push_back(pose);
      pp_this = std::get<3>(correspondences[i]);
      std::cout << "Principal Point: " << pp_this[0] << " " << pp_this[1] << std::endl; 
    }
    // if(points3D_this.size() < 20){
    //   std::cout << "points3D_this size: " << points3D_this.size() << std::endl; 
    //   continue;
    // }
    // std::cout << "Calibrating camera: " << camera_id << std::endl;
    CostMatrixOptions cm_option;
    CostMatrix cost_matrix_this = build_cost_matrix_multi(points2D_this, cm_option, pp_this);
    IntrinsicCalib intrinsic_calib = calibrate_multi(points2D_this, points3D_this, cost_matrix_this, pp_this, poses);
    debug_counter++;
    // std::cout << "debug_counter: " << debug_counter << std::endl;
    if(debug_counter>50){
      std::cout << "Breaking out of potentially infinite loop after 100 iterations" << std::endl;
        break;
    }
    std::vector<double> radii;
    std::vector<double> theta;
    std::vector<double> focal_lengths;
    std::vector<double> principal_point;
    for (const auto& pair : intrinsic_calib.r_f) {
        radii.push_back(pair.first); // Extract the first element of each pair (r)
        // std::cout << "r: " << pair.first << " f: " << pair.second << std::endl;
        focal_lengths.push_back(pair.second); // Extract the second element of each pair (f)
    }
    for (const auto& pair : intrinsic_calib.theta_r) {
        theta.push_back(pair.first); // Extract the first element of each pair (r)
    }
    principal_point = {intrinsic_calib.pp.x(), intrinsic_calib.pp.y()};
    // std::cout << "Principal Point: " << principal_point[0] << " " << principal_point[1] << std::endl;
    // save as a file
    // std::string filename = "calibration_" + std::to_string(camera_id) + "_" + currentDateTime() + ".txt";
    // std::ofstream
    // file(filename);
    // if (file.is_open()) {
    //   file << "Principal Point: " << principal_point[0] << " " << principal_point[1] << std::endl;
    //   file.close();
    // }
    

    // // using gt focal lengths
    // double fx = 3411.42;
    // double fy = 3410.02;
    // double cx = 3041.29;
    // double cy = 2014.07;
    // double k1 = 0.21047;
    // double k2 = 0.21102;
    // double p1 = -5.36231e-06;
    // double p2 =  0.00051541;
    // double k3 = -0.158023;
    // double k4 =  0.406856;
    // double s1 =  -8.46499e-05;
    // double s2 =  0.000861313;
    // std::vector<double> r_3d;
    // // from 0 to 1, increment by 0.01
    // for (double i = 0.002; i < 1; i += 0.002) {
    //   r_3d.push_back(i);
    // }
    // std::vector<double> radii_list;
    // std::vector<double> gt_focal_lengths;

    // for (int i = 0; i < r_3d.size(); i++){
    //   double theta = atan(r_3d[i]);
    //   double x_d = r_3d[i] * cos(theta);
    //   double y_d = r_3d[i] * sin(theta);
    //   double u_d = theta/r_3d[i] * x_d;
    //   double v_d = theta/r_3d[i] * y_d;
    //   double theta_squared = theta * theta;
    //   double t_r = 1 + k1 * theta_squared + k2 * theta_squared * theta_squared + k3 * theta_squared * theta_squared * theta_squared + k4 * theta_squared * theta_squared * theta_squared * theta_squared;
    //   double u_n = u_d * t_r + 2 * p1 * u_d * v_d + p2 * (theta_squared + 2 * u_d * u_d) + s1 * theta_squared;
    //   double v_n = v_d * t_r + 2 * p2 * u_d * v_d + p1 * (theta_squared + 2 * v_d * v_d) + s2 * theta_squared;
    //   double u = fx * u_n + cx;
    //   double v = fy * v_n + cy;
    //   double radius = sqrt((u - cx) * (u - cx) + (v - cy) * (v - cy));
    //   radii_list.push_back(radius);
    //   double gt_f = (u-cx)/x_d;
    //   gt_focal_lengths.push_back(gt_f);
    // }
   
    if(radii.size() > 0.7*camera_this.GetRawRadii().size() && radii[0] > 0){
    std::cout<<"updated calibration"<<std::endl;
    camera_this.SetRawRadii(radii);
    // camera_this.SetRawRadii(radii_list);
    camera_this.SetTheta(theta);
    camera_this.SetFocalLengthParams(focal_lengths);
    // camera_this.SetFocalLengthParams(gt_focal_lengths);
    // std::cout << "Calibrated camera: " << camera_id << std::endl;
    camera_this.FitSpline(radii, focal_lengths);
    camera_this.FitSpline_theta_r(theta, radii, principal_point);
    // camera_this.FitSpline(theta,radii);
    // camera_this.FitSpline(radii_list, gt_focal_lengths);
    // camera_this.FitPIeceWiseSpline_binary(radii, focal_lengths);
    }
    // else{
    //   standard_triangulation = false;
    // }
  }
  }
  }
  // Extract correspondences without an existing triangulated observation.
  std::vector<CorrData> create_corrs_data;
  // std::cout<<"standard triangulation: "<<standard_triangulation<<std::endl;
  create_corrs_data.reserve(corrs_data.size());
  for (const CorrData& corr_data : corrs_data) {
    if (!corr_data.point2D->HasPoint3D()) {
      create_corrs_data.push_back(corr_data);
    }
    
  }
  // std::cout << "Initial in IncrementalTriangulator::Create: " << initial << std::endl;
  // std::cout << "create_corrs_data size: " << create_corrs_data.size() << std::endl;

  if (create_corrs_data.size() < 2) {
    // Need at least two observations for triangulation.
    // std::cout << "Returning 0" << std::endl;
    return 0;
  } else if (options.ignore_two_view_tracks && create_corrs_data.size() == 2) {
    const CorrData& corr_data1 = create_corrs_data[0];
    if (correspondence_graph_->IsTwoViewObservation(corr_data1.image_id,
                                                    corr_data1.point2D_idx)) {
      return 0;
    }
  }
  // std::cout << "Setting up data " << std::endl;
  // Setup data for triangulation estimation.
  std::vector<TriangulationEstimator::PointData> point_data;
  point_data.resize(create_corrs_data.size());
  std::vector<TriangulationEstimator::PoseData> pose_data;
  pose_data.resize(create_corrs_data.size());
  int count = 0;
  int radial_count = 0;
  for (size_t i = 0; i < create_corrs_data.size(); ++i) {
    
    const CorrData& corr_data = create_corrs_data[i];
    point_data[i].point = corr_data.point2D->XY();
    point_data[i].point_normalized =
        corr_data.camera->ImageToWorld(point_data[i].point);
    pose_data[i].proj_matrix = corr_data.image->ProjectionMatrix();
    pose_data[i].proj_center = corr_data.image->ProjectionCenter();
    
    pose_data[i].camera = corr_data.camera;
    
    pose_data[i].image = corr_data.image;
    std::vector<double> raw_radii = corr_data.camera->GetRawRadii();
    std::vector<double> thetas = corr_data.camera->GetTheta();
    // std::cout<<"raw_radii size: "<<raw_radii.size()<<std::endl;
    if(raw_radii.size() <20 ||(raw_radii.size() > 20 && raw_radii[0]<0)){
      // std::cout<<"Standard Triangulation not possible"<<std::endl;
      standard_triangulation = false;
    }

    if(standard_triangulation){
    // std::cout<<"Standard Triangulation possible by raw radii"<<std::endl;
    double radius = point_data[i].point_normalized.norm();
    // generate uniformly distributed 100 points in raw_radii range
    // std::vector<double> raw_radii = corr_data.image->GetRawRadii();
    // std::vector<double> focal_lengths = corr_data.image->GetFocalLengthParams();
    
    std::vector<double> focal_lengths = corr_data.camera->GetFocalLengthParams();
    std::vector<double> std_points;
    std::string timestamp = currentDateTime();
    std::string filename ="radii_focal_lengths_" + std::to_string(corr_data.image_id) + "_" + timestamp + ".txt";
    std::vector<double> focal_lengths_splined;
          
    // std::cout<<"linspacing"<<std::endl;
    Eigen::VectorXd points = Eigen::VectorXd::LinSpaced(100, raw_radii[0], raw_radii[raw_radii.size() - 1]);
    // std::cout<<"linspaced"<<std::endl;  
    std_points = std::vector<double>(points.data(), points.data() + points.size());
    
    for (int i = 0; i < std_points.size(); i++) {
      // focal_lengths_splined.push_back(corr_data.image->EvalFocalLength(std_points[i]));
      focal_lengths_splined.push_back(corr_data.camera->EvalFocalLength(std_points[i]));
    }
    // std::cout << "Radius: " << radius << std::endl;

    // calculate the focal length by interpolating the focal_lengths
    // sort the radii and focal lengths
    // sortRadiiAndFocalLengths(raw_radii, focal_lengths);
    // std::vector<double> smoothed_focal_lengths = movingAverage(focal_lengths, 10);
    
    
    // save the sorted radii and focal lengths as txt file
    

    // if(corr_data.image_id == 93 || corr_data.image_id == 107){
    if(corr_data.image_id == 0 ){
    std::ofstream file(filename);
    for (int i = 0; i < raw_radii.size(); i++) {
      file << raw_radii[i] << " " << focal_lengths[i] << " "<<thetas[i] << std::endl;
    }
    file.close();
    }
    
    // find the two radii that the current radius is between
    double focal_length = 20;
    for (int i = 0; i < raw_radii.size() - 1; i++) {
      if (radius >= raw_radii[i] && radius <= raw_radii[i+1]) {
        // interpolate the focal length
        focal_length = focal_lengths[i] + (focal_lengths[i+1] - focal_lengths[i]) * (radius - raw_radii[i]) / (raw_radii[i+1] - raw_radii[i] + 1e-6);
        break;
      }
    }
    // if the radius is smaller than the smallest radius, set the focal length to the smallest focal length
    if (radius < raw_radii[0]) {
      focal_length = focal_lengths[0];
      count++;
    }
    // if the radius is larger than the largest radius, set the focal length to the largest focal length
    if (radius > raw_radii[raw_radii.size() - 1]) {
      focal_length = focal_lengths[raw_radii.size() - 1];
      count++;
    }

    // double focal_length_splined = corr_data.image->EvalFocalLength(radius);
    double focal_length_splined = corr_data.camera->EvalFocalLength(radius);
    // std::cout << "Estimated Focal Length: " << focal_length << std::endl;
    // save the interpolated focal length to the previous txt file
    
    // if(corr_data.image_id == 93 || corr_data.image_id == 107){
    if(corr_data.image_id == 0 ){
    std::ofstream file2(filename, std::ios_base::app);
    file2 << "Interpolated Focal Length: " << focal_length << std::endl;
    file2 << "Focal Length Spline: " << focal_length_splined << std::endl;
    file2 <<"radii: " << radius << std::endl;
    file2.close();
    }
    // if(corr_data.image_id == 93 || corr_data.image_id == 107){
    if(corr_data.image_id == 0 ){
    std::ofstream file2(filename, std::ios_base::app);
    for (int i = 0; i < std_points.size(); i++) {
      file2 <<"spline: "<<std_points[i] << " " << focal_lengths_splined[i] << std::endl;
    }
    // for (int i = 0; i < std_points.size(); i++) {
    //   file2 <<"piece_spline: "<<std_points[i] << " " << focal_lengths_splined_piecewise[i] << std::endl;
    // }
    // for (int i = 0; i < std_points.size(); i++) {
    //   file2 <<"grid_spline: "<<std_points[i] << " " << focal_lengths_splined_grid[i] << std::endl;
    // }
    file2.close();
    }
    
    
    // set the focal length
    // std::cout << "Interpolated Point wise Focal length: " << focal_length << std::endl;
    if (focal_length > 0) {
    point_data[i].focal_length = focal_length_splined;
    }else{
      point_data[i].focal_length = focal_length_splined;
      // count++;
      // standard_triangulation = false;
    }
    
    point_data[i].point_normalized_standard = point_data[i].point_normalized / point_data[i].focal_length;
    Eigen::Matrix3d K;
    K << focal_length_splined, 0, pose_data[i].camera->PrincipalPointX(),
         0, focal_length_splined, pose_data[i].camera->PrincipalPointY(),
         0, 0, 1;
    pose_data[i].proj_matrix_standard = K * pose_data[i].proj_matrix;
  }
  }
  if (count > 0) {
    full_error = false;
  }
  if(radial_count >= 0.5*create_corrs_data.size()){
    standard_triangulation = false;
  }
  // Setup estimation options for radial estimation
  EstimateTriangulationOptions tri_options;
  tri_options.min_tri_angle = DegToRad(options.min_angle);
  tri_options.residual_type =
      TriangulationEstimator::ResidualType::ANGULAR_ERROR;
  // if(full_error && standard_triangulation){
  // if(standard_triangulation){
  //   tri_options.residual_type = TriangulationEstimator::ResidualType::ANGULAR_ERROR_SPLITTING;
  // // //   // std::cout << "================== Using full error=======================" << std::endl;}
  // }
  tri_options.ransac_options.max_error =
      DegToRad(options.create_max_angle_error);
  tri_options.ransac_options.confidence = 0.9999;
  tri_options.ransac_options.min_inlier_ratio = 0.02;
  tri_options.ransac_options.max_num_trials = 10000;

  // Enforce exhaustive sampling for small track lengths.
  const size_t kExhaustiveSamplingThreshold = 15;
  if (point_data.size() <= kExhaustiveSamplingThreshold) {
    tri_options.ransac_options.min_num_trials = NChooseK(point_data.size(), 2);
  }

  // Setup estimation options for full estimation
  EstimateTriangulationOptions tri_options_full;
  tri_options_full.min_tri_angle = DegToRad(options.min_angle);
  tri_options_full.residual_type =
      TriangulationEstimator::ResidualType::ANGULAR_ERROR_SPLITTING;
  // tri_options_full.residual_type =
      // TriangulationEstimator::ResidualType::ANGULAR_ERROR;
  tri_options_full.ransac_options.max_error =
      DegToRad(options.create_max_angle_error);
  tri_options_full.ransac_options.confidence = 0.9999;
  tri_options_full.ransac_options.min_inlier_ratio = 0.02;
  tri_options_full.ransac_options.max_num_trials = 10000;

  // Enforce exhaustive sampling for small track lengths.
  
  if (point_data.size() <= kExhaustiveSamplingThreshold) {
    tri_options_full.ransac_options.min_num_trials = NChooseK(point_data.size(), 2);
  }
  // std::cout<<"Standard_triangulation: "<<standard_triangulation<<std::endl; 

  // Estimate triangulation.
  Eigen::Vector3d xyz;
  std::vector<char> inlier_mask;
  // standard_triangulation = false;
  if (!EstimateTriangulation(tri_options, point_data, pose_data, &inlier_mask,
                             &xyz, initial, false)) {
    return 0;
  }
  Eigen::Vector3d xyz_full;
  std::vector<char> inlier_mask_full;

  if(standard_triangulation){
    EstimateTriangulation(tri_options_full, point_data,pose_data, &inlier_mask_full,
                          &xyz_full, initial, standard_triangulation);
  }
  // standard_triangulation = true;

  // Add inliers to estimated track.
  int num_constraints = 0;
  Track track;
  track.Reserve(create_corrs_data.size());
  for (size_t i = 0; i < inlier_mask.size(); ++i) {
    if (inlier_mask[i]) {
      const CorrData& corr_data = create_corrs_data[i];
      track.AddElement(corr_data.image_id, corr_data.point2D_idx);
      if(standard_triangulation){
        num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 2 : 2;
        // num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 1 : 2;
      }
      else{
        num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 1 : 2;
      }
    }
  }
  int num_constraints_full = 0;
  Track track_full;

  if(standard_triangulation){
    track_full.Reserve(create_corrs_data.size());
    for (size_t i = 0; i < inlier_mask_full.size(); ++i) {
      if (inlier_mask_full[i]) {
        const CorrData& corr_data = create_corrs_data[i];
        track_full.AddElement(corr_data.image_id, corr_data.point2D_idx);  
        num_constraints_full += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 2 : 2;
        // num_constraints_full += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 1 : 2;
        }
      }
    }
  

  // constraint trial

  if((num_constraints < 4) && (num_constraints_full < 4)) {
    // this is a underconstrained point, we do not add it to the reconstruction since it will get filtered later anyways
    return 0;
  }
  point3D_t point3D_id;
  bool used_full = false;
  // Add estimated point to reconstruction.
  // if (track.Length()>track_full.Length()){
    point3D_id = reconstruction_->AddPoint3D(xyz, track);
  // }else{
    // point3D_id = reconstruction_->AddPoint3D(xyz_full,track_full);
    // used_full = true;
    // std::cout << " ------- used full error track ------- " << std::endl;
  // }
  modified_point3D_ids_.insert(point3D_id);

  const size_t kMinRecursiveTrackLength = 3;
  if(!used_full){
  if (create_corrs_data.size() - track.Length() >= kMinRecursiveTrackLength) {
    update_calibration = false;
    return track.Length() + Create(options, create_corrs_data, initial, standard_triangulation, update_calibration);
  }

  return track.Length();} else{
    if (create_corrs_data.size() - track_full.Length() >= kMinRecursiveTrackLength) {
      update_calibration = false;
    return track_full.Length() + Create(options, create_corrs_data, initial, standard_triangulation, update_calibration);
  }

  return track_full.Length();
  }
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
          // if (CalculateSquaredReprojectionErrorFinal(
          //         test_point2D.XY(), merged_xyz, test_image.Qvec(),
          //         test_image.Tvec(), test_image.GetRawRadii(),test_image.GetFocalLengthParams(),test_camera) > max_squared_reproj_error) {
          //   merge_success = false;
          //   break;
          // }
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
        // if (CalculateSquaredReprojectionErrorFinal(
        //         point2D.XY(), point3D.XYZ(), image.Qvec(), image.Tvec(), image.GetRawRadii(), image.GetFocalLengthParams(), image.GetTheta(),
        //         camera) > max_squared_reproj_error) {
        //   continue;
        // }

        // Success, add observation to point track.
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
