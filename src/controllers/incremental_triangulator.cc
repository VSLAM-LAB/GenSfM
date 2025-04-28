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
    reconstruction_(reconstruction) {
  }

  size_t IncrementalTriangulator::TriangulateImage(const Options& options,
    const image_t image_id, bool initial, bool standard_triangulation) {
    CHECK(options.Check());

    size_t num_tris = 0;

    ClearCaches();
    const Image& image = reconstruction_->Image(image_id);
    if (!image.IsRegistered()) {
      return num_tris;
    }

    const Camera& camera = reconstruction_->Camera(image.CameraId());


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
      if (camera.GetRawRadii().size() < 20) {
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
      }
      else {
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
        if (standard_triangulation) {
          std::vector<double> raw_radii = corr_data.camera->GetRawRadii();
          std::vector<double> focal_lengths = corr_data.camera->GetFocalLengthParams();
          double radius = point_data[i].point_normalized.norm();
          double focal_length = 0;
          if (raw_radii.size() == 0) {
            continue;
          }
          for (int i = 0; i < raw_radii.size() - 1; i++) {
            if (radius >= raw_radii[i] && radius <= raw_radii[i + 1]) {
              // interpolate the focal length
              focal_length = focal_lengths[i] + (focal_lengths[i + 1] - focal_lengths[i]) * (radius - raw_radii[i]) / (raw_radii[i + 1] - raw_radii[i] + 1e-6);
              break;
            }
          }
          // if the radius is smaller than the smallest radius, set the focal length to the smallest focal length
          if (radius < raw_radii[0]) {
            focal_length = focal_lengths[0];
          }
          // if the radius is larger than the largest radius, set the focal length to the largest focal length
          if (radius > raw_radii[raw_radii.size() - 1]) {
            focal_length = focal_lengths[raw_radii.size() - 1];
          }
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
          if (standard_triangulation) {
            num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 2 : 2;
          }
          else {
            num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 1 : 2;
          }

        }
      }

      if (num_constraints < 4) {
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
    // check the number of registered images per camera
    std::map<camera_t, int> num_registrations_per_camera;
    bool standard_triangulation = (num_registrations >= options.min_num_reg_images);
    for (const auto& image : reconstruction_->Images()) {
      if (image.second.IsRegistered()) {
        num_registrations_per_camera[image.second.CameraId()]++;
      }
    }
    for (const auto& pair : num_registrations_per_camera) {
      if (pair.second < options.min_num_reg_images) {
        bool standard_triangulation = false;
        break;
      }
    }


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

      // Find correspondences and perform retriangulation.

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
          const std::vector<CorrData> corrs_data1 = { corr_data1 };
          num_tris += Continue(re_options, corr_data2, corrs_data1);
        }
        else if (!point2D1.HasPoint3D() && point2D2.HasPoint3D()) {
          const std::vector<CorrData> corrs_data2 = { corr_data2 };
          num_tris += Continue(re_options, corr_data1, corrs_data2);
        }
        else if (!point2D1.HasPoint3D() && !point2D2.HasPoint3D()) {
          const std::vector<CorrData> corrs_data = { corr_data1, corr_data2 };
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
      }
      else {
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

  void SortRadiiAndFocalLengths(std::vector<double>& raw_radii, std::vector<double>& focal_lengths) {
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

  std::vector<double> MovingAverage(const std::vector<double>& data, int window_size) {
    std::vector<double> result;
    int n = data.size();
    if (n == 0 || window_size <= 0) return result;

    double window_sum = 0;
    int window_start = 0, window_end = 0;
    for (window_end = 0; window_end < window_size && window_end < n; ++window_end) {
      window_sum += data[window_end];
    }

    result.push_back(window_sum / window_size);
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

    bool full_error = true;
    int debug_counter = 0;
    if (standard_triangulation) {
      std::map<camera_t, std::vector<std::tuple<CameraPose, std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector3d>>>> camera_correspondences;

      for (const CorrData& corr_data : corrs_data) {

        const Image& image = *corr_data.image;
        const Camera& camera = *corr_data.camera;
        camera_t camera_id = camera.CameraId();

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

        CameraPose pose(image.Qvec(), image.Tvec());


        // filter the outliers
        if (update_calibration) {
          PoseRefinementOptions pose_refinement_options;
          CostMatrixOptions cm_options;
          CostMatrix cost_matrix = build_cost_matrix(points2D_xy, cm_options, Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY()));
          Eigen::Vector2d pp = Eigen::Vector2d(camera.PrincipalPointX(), camera.PrincipalPointY());
          CameraPose* pose_refining = &pose;
          Eigen::Vector2d* pp_refining = &pp;
          PoseRefinement1DRadialOptions opt;

          pose_refinement_1D_radial(points2D_xy, points3D_xyz, pose_refining, pp_refining, opt);

          CameraPose pose_refined;

          if (points3D_xyz.size() < 4) {
            pose_refined = pose;
          }
          else { pose_refined = pose_refinement(points2D_xy, points3D_xyz, cost_matrix, pp_refining, pose, pose_refinement_options); }

          filter_result_pose_refinement(points2D_xy, points3D_xyz, pose_refined, pp_refining, pose_refinement_options);
          std::vector<Eigen::Vector2d> points2D_xy_filtered = points2D_xy;
          std::vector<Eigen::Vector3d> points3D_xyz_filtered = points3D_xyz;
          if (camera_correspondences.find(camera_id) == camera_correspondences.end()) {
            camera_correspondences[camera_id] = std::vector<std::tuple<CameraPose, std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector3d>, Eigen::Vector2d>>();
            camera_correspondences[camera_id].push_back(std::make_tuple(pose_refined, points2D_xy, points3D_xyz, pp_refining));
          }
          else {
            camera_correspondences[camera_id].push_back(std::make_tuple(pose_refined, points2D_xy, points3D_xyz, pp_refining));
          }
        }

      }

      debug_counter++;
      if (update_calibration) {

        for (const auto& pair : camera_correspondences) {
          camera_t camera_id = pair.first;
          Camera& camera_this = reconstruction_->Camera(camera_id);

          std::vector<std::tuple<CameraPose, std::vector<Eigen::Vector2d>, std::vector<Eigen::Vector3d>> > correspondences = pair.second;
          std::vector<std::vector<Eigen::Vector2d>> points2D_this;
          std::vector<std::vector<Eigen::Vector3d>> points3D_this;
          std::vector<CameraPose> poses;
          Eigen::Vector2d pp_this;
          for (int i = 0; i < correspondences.size(); i++) {
            CameraPose pose = std::get<0>(correspondences[i]);
            std::vector<Eigen::Vector2d> points2D_i = std::get<1>(correspondences[i]);
            std::vector<Eigen::Vector3d> points3D_i = std::get<2>(correspondences[i]);
            pp_this = std::get<3>(correspondences[i]);
            points2D_this.push_back(points2D_i);
            points3D_this.push_back(points3D_i);
            poses.push_back(pose);
          }

          CostMatrixOptions cm_option;
          CostMatrix cost_matrix_this = build_cost_matrix_multi(points2D_this, cm_option, pp_this);
          IntrinsicCalib intrinsic_calib = calibrate_multi(points2D_this, points3D_this, cost_matrix_this, pp_this, poses);
          debug_counter++;
          if (debug_counter > 50) {
            std::cout << "Breaking out of potentially infinite loop after 100 iterations" << std::endl;
            break;
          }
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
          if (radii.size() > 20 && radii[0] > 0) {

            camera_this.SetRawRadii(radii);
            camera_this.SetTheta(theta);
            camera_this.SetFocalLengthParams(focal_lengths);
            std::cout << "Calibrated camera: " << camera_id << std::endl;
            camera_this.FitSpline(radii, focal_lengths);

            camera_this.FitPieceWiseSpline(radii, focal_lengths);


          }
          else {
            standard_triangulation = false;
          }
        }
      }
    }
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
      // std::cout << "Returning 0" << std::endl;
      return 0;
    }
    else if (options.ignore_two_view_tracks && create_corrs_data.size() == 2) {
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
      if (raw_radii.size() < 20 || (raw_radii.size() > 20 && raw_radii[0] < 0)) {
        standard_triangulation = false;
      }

      if (standard_triangulation) {

        double radius = point_data[i].point_normalized.norm();
        std::vector<double> focal_lengths = corr_data.camera->GetFocalLengthParams();
        std::vector<double> std_points;
        std::string timestamp = currentDateTime();
        std::string filename = "radii_focal_lengths_" + std::to_string(corr_data.image_id) + "_" + timestamp + ".txt";
        std::vector<double> focal_lengths_splined;
        std::vector<double> focal_lengths_splined_piecewise;
        std::vector<double> focal_lengths_splined_grid;
        Eigen::VectorXd points = Eigen::VectorXd::LinSpaced(100, raw_radii[0], raw_radii[raw_radii.size() - 1]);
        std_points = std::vector<double>(points.data(), points.data() + points.size());

        for (int i = 0; i < std_points.size(); i++) {
          focal_lengths_splined.push_back(corr_data.camera->EvalFocalLength(std_points[i]));
          focal_lengths_splined_piecewise.push_back(corr_data.camera->EvalPieceFocalLength(std_points[i]));
          focal_lengths_splined_grid.push_back(corr_data.camera->EvalGridFocalLength(std_points[i]));
        }

        // std::ofstream file(filename);
        // for (int i = 0; i < raw_radii.size(); i++) {
        //   file << raw_radii[i] << " " << focal_lengths[i] << " "<<thetas[i] << std::endl;
        // }
        // file.close();
        double focal_length = 0;
        for (int i = 0; i < raw_radii.size() - 1; i++) {
          if (radius >= raw_radii[i] && radius <= raw_radii[i + 1]) {
            // interpolate the focal length
            focal_length = focal_lengths[i] + (focal_lengths[i + 1] - focal_lengths[i]) * (radius - raw_radii[i]) / (raw_radii[i + 1] - raw_radii[i] + 1e-6);
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
        double focal_length_splined = corr_data.camera->EvalFocalLength(radius);

        // std::ofstream file2(filename, std::ios_base::app);
        // file2 << "Interpolated Focal Length: " << focal_length << std::endl;
        // file2 << "Focal Length Spline: " << focal_length_splined << std::endl;
        // file2 <<"radii: " << radius << std::endl;
        // file2.close();
        // }
        // if(corr_data.image_id == 93 || corr_data.image_id == 107){
        // std::ofstream file2(filename, std::ios_base::app);
        // for (int i = 0; i < std_points.size(); i++) {
        //   file2 <<"spline: "<<std_points[i] << " " << focal_lengths_splined[i] << std::endl;
        // }
        // for (int i = 0; i < std_points.size(); i++) {
        //   file2 <<"piece_spline: "<<std_points[i] << " " << focal_lengths_splined_piecewise[i] << std::endl;
        // }
        // for (int i = 0; i < std_points.size(); i++) {
        //   file2 <<"grid_spline: "<<std_points[i] << " " << focal_lengths_splined_grid[i] << std::endl;
        // }
        // file2.close();
        // }




        // using an interpolated focal length 
        if (focal_length > 0) {
          point_data[i].focal_length = focal_length;
        }
        else {
          point_data[i].focal_length = focal_length;
        }

        point_data[i].point_normalized_standard = point_data[i].point_normalized / point_data[i].focal_length;
        Eigen::Matrix3d K;
        K << focal_length, 0, pose_data[i].camera->PrincipalPointX(),
          0, focal_length, pose_data[i].camera->PrincipalPointY(),
          0, 0, 1;
        pose_data[i].proj_matrix_standard = K * pose_data[i].proj_matrix;
      }
    }
    if (count > 0) {
      full_error = false;
    }
    if (radial_count >= 0.5 * create_corrs_data.size()) {
      standard_triangulation = false;
    }

    // Setup estimation options for radial estimation
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

    // Estimate triangulation.
    Eigen::Vector3d xyz;
    std::vector<char> inlier_mask;
    if (!EstimateTriangulation(tri_options, point_data, pose_data, &inlier_mask,
      &xyz, initial, standard_triangulation)) {
      return 0;
    }
    Eigen::Vector3d xyz_full;
    std::vector<char> inlier_mask_full;

    if (standard_triangulation) {
      EstimateTriangulation(tri_options_full, point_data, pose_data, &inlier_mask_full,
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
        if (standard_triangulation) {
          num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 2 : 2;
        }
        else {
          num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 1 : 2;
        }
      }
    }
    int num_constraints_full = 0;
    Track track_full;

    if (standard_triangulation) {
      track_full.Reserve(create_corrs_data.size());
      for (size_t i = 0; i < inlier_mask_full.size(); ++i) {
        if (inlier_mask_full[i]) {
          const CorrData& corr_data = create_corrs_data[i];
          track_full.AddElement(corr_data.image_id, corr_data.point2D_idx);
          num_constraints += (corr_data.camera->ModelId() == Radial1DCameraModel::model_id || corr_data.camera->ModelId() == ImplicitDistortionModel::model_id) ? 2 : 2;
        }
      }
    }


    // constraint trial

    if ((num_constraints < 4) && (num_constraints_full < 4)) {
      // this is a underconstrained point, we do not add it to the reconstruction since it will get filtered later anyways
      return 0;
    }
    point3D_t point3D_id;
    bool used_full = false;
    // Add estimated point to reconstruction.
    if (track.Length() >= track_full.Length()) {
      point3D_id = reconstruction_->AddPoint3D(xyz, track);
    }
    else {
      point3D_id = reconstruction_->AddPoint3D(xyz_full, track_full);
      used_full = true;
    }
    modified_point3D_ids_.insert(point3D_id);

    const size_t kMinRecursiveTrackLength = 3;
    if (!used_full) {
      if (create_corrs_data.size() - track.Length() >= kMinRecursiveTrackLength) {
        return track.Length() + Create(options, create_corrs_data, initial, standard_triangulation);
      }

      return track.Length();
    }
    else {
      if (create_corrs_data.size() - track_full.Length() >= kMinRecursiveTrackLength) {
        return track_full.Length() + Create(options, create_corrs_data, initial, standard_triangulation);
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
        for (const Track* track : { &point3D.Track(), &corr_point3D.Track() }) {
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
          }
          else {
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

          if (CalculateSquaredReprojectionError(
            point2D.XY(), point3D.XYZ(), image.Qvec(), image.Tvec(),
            camera) > max_squared_reproj_error) {
            continue;
          }
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
    }
    else {
      return it->second;
    }
  }

}  // namespace colmap
