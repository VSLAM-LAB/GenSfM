#ifndef GET_INTRINSICS_CC_
#define GET_INTRINSICS_CC_
#include <iostream>
#include <vector>
#include <cassert>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>

#include <filesystem>
// #include "base/projection.h"
// #include "estimators/triangulation.h"
// #include "estimators/implicit_camera_pose.h"
// #include "estimators/implicit_cost_matrix.h"
// #include "estimators/implicit_pose_refinement.h"
#include "util/math.h"
// #include "util/misc.h"
#include "spline.h"

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;
    if (num <= 0) return result;
    if (num == 1) {
        result.push_back(start);
        return result;
    }

    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + step * i);
    }
    return result;
}

tk::spline<double> ransac_spline(int max_iteration, int degree, double threshold, std::vector<double>radii, std::vector<double>focal_lengths) {

  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, radii.size() - 1);
  tk::spline<double> best_spline;
  std::vector<int> indices;

  for (int i = 0; i < max_iteration; ++i){
    std::vector<std::pair<double, double>> samples;
    indices.clear();
    indices.push_back(0);
    indices.push_back(radii.size()-1);
    for (int j = 0; j < (degree -2); ++j) {
      int idx = dis(gen);
      while(std::find(indices.begin(), indices.end(), idx) != indices.end()){
          idx = dis(gen);
      }
      indices.push_back(idx);
    }
    for (int j = 0; j < indices.size(); j++) {
      samples.emplace_back(radii[indices[j]], focal_lengths[indices[j]]);
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
    for (size_t j = 0; j < radii.size(); ++j) {
      double y_est = s(radii[j]);
      if (fabs(y_est - focal_lengths[j]) < threshold) {
          ++inliers_count;
      }
    }
    if (inliers_count >= best_inliers_count) {
      best_inliers_count = inliers_count;
      best_spline = s;
    }
  }

  return best_spline;
}

void recursiveSplit(const std::vector<double>& radii, const std::vector<double>& focal_lengths,
                        std::vector<std::vector<double>>& radii_segments,
                        std::vector<std::vector<double>>& focal_lengths_segments, 
                        double threshold, double stddev_threshold){
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
        if (max_gap < threshold){
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





std::vector<double> IdentifyCalibratedArea(std::vector<double>& radii, std::vector<double>& focal_lengths){
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

  std::cout << "!!! Original threshold: " << threshold;
  threshold = std::max(threshold, colmap::DegToRad(0.1));
  std::cout << ", New threshold: " << threshold << std::endl;
  recursiveSplit(new_radii, new_focal_lengths, radii_segments, focal_lengths_segments, threshold, std_threshold);
  int longest_segment = 0;
  int longest_segment_size = 0;
  for(int i = 0; i < radii_segments.size(); i++){
    if(radii_segments[i].size() > longest_segment_size){
      longest_segment = i;
      longest_segment_size = radii_segments[i].size();
    }
  }
  std::vector<double> calibrated_area = {radii_segments[longest_segment].front(), radii_segments[longest_segment].back()};
  return calibrated_area;
}


tk::spline<double> FitPieceWiseSpline_binary(std::vector<double>& radii, std::vector<double>& focal_lengths){
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
  threshold = std::max(threshold, colmap::DegToRad(0.1));
  std::cout << ", New threshold: " << threshold << std::endl;
  recursiveSplit(new_radii, new_focal_lengths, radii_segments, focal_lengths_segments, threshold, std_threshold);

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
  int degree = 10;
  double threshold_ransac = 5.0;
  tk::spline<double> best_spline = ransac_spline(max_it, degree, threshold_ransac, radii_calibrated, focal_lengths_calibrated);
return best_spline;

}


int main(){
    

    // the base path to the files
    std::string base_path = "/home/linpan/workspace/implicit_radial_sfm/scripts/babelcalib/eval/";
    std::string output_base = "/home/linpan/workspace/implicit_radial_sfm/scripts/babelcalib/spline/";
    // list of files to read
    for (const auto& entry : std::filesystem::directory_iterator(base_path)){
        std::vector<double> theta{};
        std::vector<double> radii{};
        std::string file_path = entry.path();
        std::ifstream file(file_path);


        std::string line;
        while (std::getline(file, line)){
            std::istringstream iss(line);
            double theta_, radii_;
            if (!(iss >> theta_ >> radii_)) { break; }
            theta.push_back(theta_);
            // std::cout << theta_ << std::endl;
            radii.push_back(radii_);
        }
        std::vector<double> calibrated_area = IdentifyCalibratedArea(theta, radii);
        tk::spline<double> best_spline = FitPieceWiseSpline_binary(theta, radii);
        std::vector<double> sample_x = best_spline.get_x();
        std::vector<double> sample_y = best_spline.get_y();
        // genrate interpolated sample_x to output
       
        // write sample_x and sample_y to a txt file
        std::string base_name = entry.path().stem().string();

        // Remove "_theta_r" if it exists in the filename
        std::string to_remove = "_theta_r";
        size_t pos = base_name.find(to_remove);
        if (pos != std::string::npos) {
            base_name.erase(pos, to_remove.length());
        }
        
        // Create a unique output file path by appending "_spline.txt" to the modified base name
        std::string output_file_path = output_base + base_name + "_spline.txt";
        
        std::ofstream output_file(output_file_path);
        for (int i = 0; i < sample_x.size(); i++){
            output_file << sample_x[i] << " " << sample_y[i] << std::endl;
        }
      }


}

#endif