#ifndef BASE_PIECE_WISE_CAMERA_TEST_CC_
#define BASE_PIECE_WISE_CAMERA_TEST_CC_
#include <iostream>
#include <vector>
#include <cassert>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "spline.h" 
tk::spline ransac_spline(int max_iteration, int degree, double threshold, std::vector<double>radii, std::vector<double>focal_lengths){
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, radii.size() - 1);
  tk::spline best_spline;
  std::vector<int> indices;
  for (int i = 0; i < max_iteration; ++i){
    std::vector<std::pair<double, double>> samples;
    indices.clear();
    for (int j = 0; j < degree; ++j) {
      int idx = dis(gen);
    //   ensure that the samples are unique
        while(std::find(indices.begin(), indices.end(), idx) != indices.end()){
            idx = dis(gen);
        }
        indices.push_back(idx);
      samples.emplace_back(radii[idx], focal_lengths[idx]);
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
tk::spline ransac_spline_binary(int max_iteration, int degree, double threshold, std::vector<double>radii, std::vector<double>focal_lengths){
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, radii.size() - 1);
  tk::spline best_spline;
  std::vector<int> indices;
  for (int i = 0; i < max_iteration; ++i){
    std::vector<std::pair<double, double>> samples;
    indices.clear();
    indices.push_back(0);
    indices.push_back(radii.size()-1);
    for (int j = 0; j < degree-2; ++j) {
      int idx = dis(gen);
    //   ensure that the samples are unique
        while(std::find(indices.begin(), indices.end(), idx) != indices.end()){
            idx = dis(gen);
        }
        indices.push_back(idx);
    }
    for(int j = 0; j < indices.size(); j++){
      samples.emplace_back(radii[indices[j]], focal_lengths[indices[j]]);
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
std::pair<tk::spline,double> ransac_spline_continuous(int max_iteration, int degree, double threshold, std::vector<double>radii, std::vector<double>focal_lengths, double first_derivative){
  int best_inliers_count = 0;
  std::vector<double> best_coeffs;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(2, radii.size() - 3);
  tk::spline best_spline;
  if(first_derivative == 0.0){
    best_spline.set_boundary(tk::spline::second_deriv, 0, tk::spline::second_deriv, 0);
  }else{
  best_spline.set_boundary(tk::spline::first_deriv, first_derivative, tk::spline::second_deriv, 0);}
  std::vector<int> indices;
  for (int i = 0; i < max_iteration; ++i){
    std::vector<std::pair<double, double>> samples;
    indices.clear();
    indices.push_back(0);
    indices.push_back(radii.size()-1);
    for (int j = 0; j < degree-2; ++j) {
      int idx = dis(gen);
    //   ensure that the samples are unique
        while(std::find(indices.begin(), indices.end(), idx) != indices.end()){
            idx = dis(gen);
        }
        indices.push_back(idx);
      // samples.emplace_back(radii[idx], focal_lengths[idx]);
    }
    for(int j = 0; j < indices.size(); j++){
      samples.emplace_back(radii[indices[j]], focal_lengths[indices[j]]);
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
  double first_derivative_best = best_spline.get_m_right_value();
  return std::make_pair(best_spline, first_derivative_best);
}
std::vector<tk::spline> ransac_piece_continuous_splines(int max_iteration, int degree, double threshold, std::vector<std::vector<double>>radii_segments, std::vector<std::vector<double>>focal_segments){
 std::vector<tk::spline> splines;
 for (size_t i = 0; i < radii_segments.size(); ++i) {
      // tk::spline s;
      std::vector<int>piece_indices;
      std::vector<double> radii_data;
      std::vector<double> focal_data;
      double first_derivative = 0.0;
      if (i == 0) {
        for(int j = 0; j < radii_segments[i].size(); j++){
          radii_data.push_back(radii_segments[i][j]);
          focal_data.push_back(focal_segments[i][j]);
        }
      }else{
        if(radii_segments.size()==1){
          continue;
        }
        radii_data.push_back(radii_segments[i-1][radii_segments[i-1].size()-1]);
        focal_data.push_back(focal_segments[i-1][focal_segments[i-1].size()-1]);
        for(int j = 0; j < radii_segments[i].size(); j++){
          radii_data.push_back(radii_segments[i][j]);
          focal_data.push_back(focal_segments[i][j]);
        }
      }

      // s.set_points(new_radii_segments[i], new_focal_lengths_segments[i]);
      // splines.push_back(s);
      if (radii_data.size() > degree) {
          std::pair<tk::spline, double> result = ransac_spline_continuous(max_iteration, degree, threshold, radii_data, focal_data, first_derivative);
          // ransac_spline_continuous(max_iteration, degree, threshold, radii_data, focal_data, first_derivative);
          tk::spline best_piece_spline = result.first;
          first_derivative = result.second;
          splines.push_back(best_piece_spline);
      }else if (radii_data.size() <= degree && radii_data.size() > 2){
        tk::spline s;
        s.set_boundary(tk::spline::first_deriv, first_derivative, tk::spline::second_deriv, 0);
        s.set_points(radii_data, focal_data);
        splines.push_back(s);
      }else{
        // std::cout << "Not enough data points to fit a spline" << std::endl;
        continue;
      }
  }
return splines;
}
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> deciding_segments(tk::spline spline, std::vector<double>radii, std::vector<double>focal_lengths){
//dynamically determine outlier threshold
std::vector<double> errors = {};
  double mean_error=0.0;
  double std_error=0.0;
  for (size_t j = 0; j < radii.size(); ++j) {
    double y_est = spline(radii[j]);
    errors.push_back(fabs(y_est - focal_lengths[j]));
    mean_error += fabs(y_est - focal_lengths[j]);
    std_error += pow(fabs(y_est - focal_lengths[j]), 2);
  }
  mean_error /= radii.size();
  std_error = sqrt(std_error / radii.size() );   
  double outlier_threshold = mean_error + 0.5 * std_error; // assuming a normal distribution of the errors, confidence interval of 95%
  std::vector<double> outlier_xs = {};
  std::vector<double> outliers_indices = {};
  for (size_t j = 0; j < radii.size(); ++j) {
    if(errors[j] > outlier_threshold){
      outlier_xs.push_back(radii[j]);
    }
  }
   for (int i = 0; i < radii.size(); i++) {
    if(abs(spline(radii[i]) - focal_lengths[i]) > outlier_threshold){
      outliers_indices.push_back(i);
    }
  }
  std::cout << "outlier_xs size: " << outlier_xs.size() << std::endl;
  std::vector<std::vector<double>> new_radii_segments = {};
  std::vector<std::vector<double>> new_focal_lengths_segments = {};
  std::vector<double> segment_radii = {} ;
  std::vector<double> segment_focal_lengths = {};
  const double outlier_fraction_threshold = 0.3;
  int outlier_count = 0;
  int segment_count = 0;
  int inlier_count = 0;   
 
  for(int i = 0; i < outliers_indices.size()-1; i++){
    if(i == 0 && outliers_indices[i] != 0){
        for (int j = 0; j < outliers_indices[i]; j++) {
            segment_radii.push_back(radii[j]);
            segment_focal_lengths.push_back(focal_lengths[j]);
        }
        new_radii_segments.push_back(segment_radii);
        new_focal_lengths_segments.push_back(segment_focal_lengths);
        segment_radii.clear();
        segment_focal_lengths.clear();
    }
    if(outliers_indices[i+1] - outliers_indices[i] <5){
        segment_radii.push_back(radii[outliers_indices[i]]);
        segment_focal_lengths.push_back(focal_lengths[outliers_indices[i]]);
    }else{
        new_radii_segments.push_back(segment_radii);
        new_focal_lengths_segments.push_back(segment_focal_lengths);
        segment_radii.clear();
        segment_focal_lengths.clear();
        for(int j = outliers_indices[i]; j < outliers_indices[i+1]; j++){
            segment_radii.push_back(radii[j]);
            segment_focal_lengths.push_back(focal_lengths[j]);
        }
        new_radii_segments.push_back(segment_radii);
        new_focal_lengths_segments.push_back(segment_focal_lengths);
        segment_radii.clear();
        segment_focal_lengths.clear();
    }
}
return std::make_pair(new_radii_segments, new_focal_lengths_segments);
}
void handleUncalibratedAreas(std::vector<double>& radii, double uncalib_threshold, double longest, std::vector<std::vector<double>>& uncalibrated_areas, double& longest_interval) {
    if (radii.empty()) return;

    double current_interval = radii.back() - radii.front();
    longest_interval = std::max(std::max(longest_interval, current_interval), uncalib_threshold * longest);

    double filter_interval = uncalib_threshold * longest;
    if (current_interval < filter_interval) {
        uncalibrated_areas.push_back(radii);
        radii.clear();
    }
}
void recursiveSplit(const std::vector<double>& radii, const std::vector<double>& focal_lengths,
                        std::vector<std::vector<double>>& radii_segments,
                        std::vector<std::vector<double>>& focal_lengths_segments, 
                        double threshold = 0.05, double stddev_threshold = 0.01) {
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
        if (max_gap < threshold || stddev < stddev_threshold) {
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

class Camera {
public:
    inline void FitPieceWiseSpline(std::vector<double>& radii, std::vector<double>& focal_lengths) const;
    inline void FitPieceWiseSpline_2(std::vector<double>& radii, std::vector<double>& focal_lengths) const;
    inline void FitPieceWiseSpline_binary(std::vector<double>& radii, std::vector<double>& focal_lengths) const;
    inline void FitGridSpline(std::vector<double>& radii, std::vector<double>& focal_lengths) const;
    mutable tk::spline spline_;
    mutable std::vector<tk::spline> piece_splines_;
    mutable std::vector<tk::spline> grid_splines_;
    mutable std::vector<std::vector<double>> intervals_;
    mutable std::vector<double> grids_;  
    inline double EvalPieceFocalLength(double radius) const;
    inline double EvalGridFocalLength(double radius) const;
};
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
//   for (int i = 0; i < new_radii.size(); i++) {
//     std::cout << "new_radii: " << new_radii[i] << std::endl;
//   }
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
    tk::spline s;
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
  std::cout << "mean_error: " << mean_error << std::endl;
  std::cout << "std_error before division: " << std_error <<"new_radii size:"<<new_radii.size() << std::endl;
  std_error = sqrt(std_error / new_radii.size() );
  std::cout << "std_error: " << std_error << std::endl;   

  double outlier_threshold = mean_error + 0.5*std_error;
  // Determinig the intervals where a new piece-wise spline should be fitted
  std::vector<double> outlier_xs;
  for (size_t j = 0; j < new_radii.size(); ++j) {
    if(errors[j] > outlier_threshold){
      outlier_xs.push_back(new_radii[j]);
    }
  }
  std::cout << "outlier_xs size: " << outlier_xs.size() << std::endl;
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
std::cout << "new_radii_segments size: " << new_radii_segments.size() << std::endl;

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
    std::cout << "intervals_ set: "  << std::endl; 
    if(intervals_.empty()){
        piece_splines_ = {best_spline};
    }

  // Fit the piece-wise splines
  std::vector<tk::spline> splines;
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
          tk::spline best_piece_spline;
          for (int j = 0; j < piecewise_max_it; ++j) {
              std::vector<std::pair<double, double>> samples_piece;
              piece_indices.clear();
              for (int k = 0; k < degree; ++k) {
                  int idx = dis_piece(gen);
                    while(std::find(piece_indices.begin(), piece_indices.end(), idx) != piece_indices.end()){
                        idx = dis_piece(gen);
                    }
                    piece_indices.push_back(idx);
                    std::cout << "idx: " << idx << std::endl;
                    std::cout << "new_radii_segments[i].size: " << new_radii_segments[i].size() << std::endl;
                  samples_piece.emplace_back(new_radii_segments[i][idx], new_focal_lengths_segments[i][idx]);
              }
              std::sort(samples_piece.begin(), samples_piece.end());
              for (int k = 0; k < piece_indices.size(); k++) {
                  std::cout << "new_radii_segments[i].size: " << new_radii_segments[i].size() << std::endl;
                  std::cout << "piece_indices: " << piece_indices[k] << std::endl;
                  std::cout << "samples_piece: " << samples_piece[k].first << std::endl;

              }
              std::vector<double> sample_x_piece, sample_y_piece;
              for (const auto& pair : samples_piece) {
                  sample_x_piece.push_back(pair.first);
                  sample_y_piece.push_back(pair.second);
              }
              tk::spline s_piece;

              s_piece.set_points(sample_x_piece, sample_y_piece);
              std::cout << "s_piece fitted" << std::endl;
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
        tk::spline s;
        s.set_points(new_radii_segments[i], new_focal_lengths_segments[i]);
        splines.push_back(s);
      }else{
        // std::cout << "Not enough data points to fit a spline" << std::endl;
        continue;
      }
  }

  if(!splines.empty()){
    piece_splines_ = splines;
}
  std::cout << "piece_splines_ size: " << piece_splines_.size() << std::endl;
 
}

inline double Camera::EvalPieceFocalLength(double radius) const {
  if (intervals_.empty() ) {
      return spline_(radius);
    }
    if(radius < intervals_.front()[0]){
        return piece_splines_.front()(radius);
    }else if (radius >= intervals_.back()[1])
    {
        return piece_splines_.back()(radius);
    }else{
    
    for (size_t i = 0; i < intervals_.size(); ++i) {
        if (radius >= intervals_[i][0] && radius < intervals_[i][1]) {
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
  std::cout << "new_radii size: " << new_radii.size() << std::endl;
//   for (int i = 0; i < new_radii.size(); i++) {
//     std::cout << "new_radii: " << new_radii[i] << std::endl;
//   }
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
  std::vector<int> indices;
  
  for (int i = 0; i < max_iterations; ++i){
    std::vector<std::pair<double, double>> samples;
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
    tk::spline s;
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
  double radial_threshold = 0.3*interval;
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
  std::cout << "grid_count: " << gird_count[0] << " " << gird_count[1] << " " << gird_count[2] << " " << gird_count[3] << std::endl;
  std::cout << "grid points:"<< new_radii[0] << " " << new_radii[0] + interval << " " << new_radii[0] + 2*interval << " " << new_radii[0] + 3*interval << " " << new_radii[0] + 4*interval << std::endl;
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
    std::vector<tk::spline> grid_splines;
    for (size_t i = 0; i < grid_points.size() - 1; ++i) {
      std::vector<double> grid_radii;
      std::vector<double> grid_focal_lengths;
      for (size_t j = 0; j < new_radii.size(); ++j) {
        if (new_radii[j] >= grid_points[i] && new_radii[j] < grid_points[i+1]) {
          grid_radii.push_back(new_radii[j]);
          grid_focal_lengths.push_back(new_focal_lengths[j]);
        }
      }
      // using ransac technique to fit the spline
      std::uniform_int_distribution<> dis_grid(0, grid_radii.size() - 1);
      int best_grid_inliers_count = 0;
      tk::spline best_grid_spline;
      std::vector<int> grid_indices;
      for (int j = 0; j < grid_max_it; ++j) {
        std::vector<std::pair<double, double>> samples_grid;
        grid_indices.clear();
        for (int k = 0; k < grid_degree; ++k) {
          int idx = dis_grid(gen);
          while(std::find(grid_indices.begin(), grid_indices.end(), idx) != grid_indices.end()){
            idx = dis_grid(gen);
          }
          grid_indices.push_back(idx);
          samples_grid.emplace_back(grid_radii[idx], grid_focal_lengths[idx]);
        }
        std::sort(samples_grid.begin(), samples_grid.end());
        std::vector<double> sample_x_grid, sample_y_grid;
        for (const auto& pair : samples_grid) {
          sample_x_grid.push_back(pair.first);
          sample_y_grid.push_back(pair.second);
        }
        tk::spline s_grid;
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
      }
      grid_splines.push_back(best_grid_spline);
    }
    grid_splines_ = grid_splines;
  }else{
    grid_splines_.push_back(best_spline); 
  }
  std::cout << "grid_splines_ size: " << grid_splines_.size() << std::endl;  
}

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

inline void Camera::FitPieceWiseSpline_2(std::vector<double>& radii, std::vector<double>& focal_lengths) const{
  assert(radii.size() == focal_lengths.size());

  // Ensure that the input data is sorted in increasing order of radii
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

  for (int i = 0; i < increasing_indices.size(); i++) {
    new_radii.push_back(radii[increasing_indices[i]]);
    new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
  }

  // use RANSAC technique to fit the original whole spline alpha
  const int max_iterations_alpha = 20;
  const double threshold_alpha = 20.0;
  int degree = 10;
  tk::spline best_spline_alpha = ransac_spline(max_iterations_alpha, degree, threshold_alpha, new_radii, new_focal_lengths);
  spline_= best_spline_alpha; // set the class member spline_ as alpha
// Identify outliers to alpha

  // Determinig the intervals where all the outliers should be discarded
  
std::vector<std::vector<double>> new_radii_segments = {};
std::vector<std::vector<double>> new_focal_lengths_segments = {};
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> segments = deciding_segments(best_spline_alpha, new_radii, new_focal_lengths);
new_radii_segments = segments.first;
new_focal_lengths_segments = segments.second;
// std::cout << "new_radii_segments size: " << new_radii_segments.size() << std::endl;
// identify the segments' sizes and lengths
std::vector<double> segment_lengths = {};
std::vector<double> segment_sizes = {};
for (int i = 0; i < new_radii_segments.size(); i++) {
    segment_sizes.push_back(new_radii_segments[i].size());
    // std::cout << "segment_size: " << new_radii_segments[i].size() << std::endl;
    double length = new_radii_segments[i].back() - new_radii_segments[i].front();
    segment_lengths.push_back(length);
    // std::cout << "segment_length: " << length << std::endl;
}
int max_segment_size = *std::max_element(segment_sizes.begin(), segment_sizes.end());
int min_segment_size = *std::min_element(segment_sizes.begin(), segment_sizes.end());
double max_segment_length = *std::max_element(segment_lengths.begin(), segment_lengths.end());
double min_segment_length = *std::min_element(segment_lengths.begin(), segment_lengths.end());
std::cout << "num segment: " << new_radii_segments.size() << std::endl;
std::cout << "total size: " << new_radii.size() << std::endl;
std::cout << "max_segment_size: " << max_segment_size << std::endl;
std::cout << "min_segment_size: " << min_segment_size << std::endl;
std::cout << "total length: "<< new_radii.back() - new_radii.front() << std::endl;
std::cout << "max_segment_length: " << max_segment_length << std::endl;
std::cout << "min_segment_length: " << min_segment_length << std::endl; 
double discard_length_threshold = 0.3*(new_radii.back() - new_radii.front())/new_radii_segments.size();
// discard the segments that are too short
std::vector<double> filtered_radii = {};
std::vector<double> filtered_focal_lengths = {};
for (int i = 0; i < new_radii_segments.size(); i++) {
    if(new_radii_segments[i].size() > 2 && segment_lengths[i] > discard_length_threshold){
        for (int j = 0; j < new_radii_segments[i].size(); j++) {
            filtered_radii.push_back(new_radii_segments[i][j]);
            filtered_focal_lengths.push_back(new_focal_lengths_segments[i][j]);
        }
    }else{
            uncalibrated_areas.push_back(new_radii_segments[i]);}
}

// fit a spline beta for the filtered data
const int max_iterations_beta = 20;
const double threshold_beta = 20.0;
tk::spline best_spline_beta= ransac_spline(max_iterations_beta, degree, threshold_beta, filtered_radii, filtered_focal_lengths);
// re-decide the pieces with beta
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> segments_beta = deciding_segments(best_spline_beta, filtered_radii, filtered_focal_lengths);
std::vector<std::vector<double>> new_radii_segments_beta = segments_beta.first;
std::vector<std::vector<double>> new_focal_lengths_segments_beta = segments_beta.second;
segment_lengths.clear();
segment_sizes.clear();
for (int i = 0; i < new_radii_segments_beta.size(); i++) {
    segment_sizes.push_back(new_radii_segments_beta[i].size());
    std::cout << "segment_size: " << new_radii_segments_beta[i].size() << std::endl;
    double length = new_radii_segments_beta[i].back() - new_radii_segments_beta[i].front();
    std::cout << "segment_length: " << length << std::endl;
    segment_lengths.push_back(length);
}
max_segment_size = *std::max_element(segment_sizes.begin(), segment_sizes.end());
min_segment_size = *std::min_element(segment_sizes.begin(), segment_sizes.end());
max_segment_length = *std::max_element(segment_lengths.begin(), segment_lengths.end());
min_segment_length = *std::min_element(segment_lengths.begin(), segment_lengths.end());
std::cout << "num segment: " << new_radii_segments_beta.size() << std::endl;
std::cout << "total size: " << filtered_radii.size() << std::endl;
std::cout << "max_segment_size: " << max_segment_size << std::endl;
std::cout << "min_segment_size: " << min_segment_size << std::endl;
std::cout << "total length: "<< filtered_radii.back() - filtered_radii.front() << std::endl;
std::cout << "max_segment_length: " << max_segment_length << std::endl;
std::cout << "min_segment_length: " << min_segment_length << std::endl;
// identify the longest segment
int longest_segment_idx = std::distance(segment_lengths.begin(), std::max_element(segment_lengths.begin(), segment_lengths.end()));
bool merge_segments = true;
bool split_segments = false;
// for (int i = 0; i < new_radii_segments_beta.size(); i++) {
//     if(i != longest_segment_idx && segment_lengths[i] > 0.3*segment_lengths[longest_segment_idx]){
//         split_segments = true;
//         merge_segments = false;
//     }
// }
// if(merge_segments){
  //merge all other segments on two sides
  std::vector<double> merged_radii_front= {};
  std::vector<double> merged_focal_lengths_front = {};
  std::vector<double> merged_radii_back = {};
  std::vector<double> merged_focal_lengths_back = {};
  if(longest_segment_idx != 0){
    for (int i = 0; i < longest_segment_idx; i++) {
      for (int j = 0; j < new_radii_segments_beta[i].size(); j++) {
        merged_radii_front.push_back(new_radii_segments_beta[i][j]);
        merged_focal_lengths_front.push_back(new_focal_lengths_segments_beta[i][j]);
      }
    }
  }else{
    for (int i = longest_segment_idx+1; i < new_radii_segments_beta.size(); i++) {
      for (int j = 0; j < new_radii_segments_beta[i].size(); j++) {
        merged_radii_back.push_back(new_radii_segments_beta[i][j]);
        merged_focal_lengths_back.push_back(new_focal_lengths_segments_beta[i][j]);
      }
    }
  }
  // calculate the longest interval of the two sides
  
  double longest_subinterval = 0.0;
  double uncalib_threshold = 0.1;
  if(merged_radii_front.size() > 0){
    std::cout<<"merged_radii_front length: "<<merged_radii_front.back()-merged_radii_front.front()<<std::endl;
  }
  
  handleUncalibratedAreas(merged_radii_front, uncalib_threshold, max_segment_length,uncalibrated_areas, longest_subinterval);
  handleUncalibratedAreas(merged_radii_back, uncalib_threshold, max_segment_length,uncalibrated_areas, longest_subinterval);
  int num_splits = 0;
  if(longest_subinterval > uncalib_threshold * max_segment_length){
    num_splits = max_segment_length/longest_subinterval;
  }
  
  std::cout << "num_splits: " << num_splits << std::endl;
  std::vector<std::vector<double>> radii_segments_post_merge_splitting = {};
  std::vector<std::vector<double>> focal_lengths_segments_post_merge_splitting = {};
  if(num_splits > 1){
    // divide the longest segment into num_splits
    if(merged_radii_front.size() > 0){
      radii_segments_post_merge_splitting.push_back(merged_radii_front);
      focal_lengths_segments_post_merge_splitting.push_back(merged_focal_lengths_front);
    }
    double split_interval = max_segment_length/num_splits;
    double start = new_radii_segments_beta[longest_segment_idx].front();
    for (int i = 0; i < num_splits; i++) {
      std::vector<double> split_radii = {};
      std::vector<double> split_focal_lengths = {};
      for (int j = 0; j < new_radii_segments_beta[longest_segment_idx].size(); j++) {
        if(new_radii_segments_beta[longest_segment_idx][j] >= start && new_radii_segments_beta[longest_segment_idx][j] < start + split_interval){
          split_radii.push_back(new_radii_segments_beta[longest_segment_idx][j]);
          split_focal_lengths.push_back(new_focal_lengths_segments_beta[longest_segment_idx][j]);
        }
      }
      radii_segments_post_merge_splitting.push_back(split_radii);
      focal_lengths_segments_post_merge_splitting.push_back(split_focal_lengths);
      start += split_interval;
    }
    if(merged_radii_back.size() > 0){
      radii_segments_post_merge_splitting.push_back(merged_radii_back);
      focal_lengths_segments_post_merge_splitting.push_back(merged_focal_lengths_back);
    }
  } else{
    if(merged_radii_front.size() > 0){
      radii_segments_post_merge_splitting.push_back(merged_radii_front);
      focal_lengths_segments_post_merge_splitting.push_back(merged_focal_lengths_front);
    }
    focal_lengths_segments_post_merge_splitting.push_back(new_focal_lengths_segments_beta[longest_segment_idx]);
    if(merged_radii_back.size() > 0){
      radii_segments_post_merge_splitting.push_back(merged_radii_back);
      focal_lengths_segments_post_merge_splitting.push_back(merged_focal_lengths_back);
    }
  }
  std::cout << "radii_segments_post_merge_splitting size: " << radii_segments_post_merge_splitting.size() << std::endl;
  for (int i = 0; i < radii_segments_post_merge_splitting.size(); i++) {
    std::cout << "interval size: " << radii_segments_post_merge_splitting[i].size() << std::endl;
    std::cout << "interval length: " << radii_segments_post_merge_splitting[i].back() - radii_segments_post_merge_splitting[i].front() << std::endl;
  }
// }
//   populate the intervals_
    intervals_.clear();
    for (size_t i = 0; i < radii_segments_post_merge_splitting.size(); ++i) {
      if (radii_segments_post_merge_splitting[i].empty()) {
        std::cerr << "Error: new_radii_segments[" << i << "] is empty." << std::endl;
        continue; // Skip empty segments
    }
      if(radii_segments_post_merge_splitting[i].size() > 2){
        std::vector<double> interval;
        interval.push_back(radii_segments_post_merge_splitting[i].front());
        interval.push_back(radii_segments_post_merge_splitting[i].back());
        intervals_.push_back(interval);
        std::cout << "interval: " << interval[0] << " " << interval[1] << std::endl;
      }
    }
    std::cout << "intervals_ set: "  << std::endl; 
    std::cout << "intervals_ size: " << intervals_.size() << std::endl;
    if(intervals_.empty()){
        piece_splines_ = {best_spline_beta};
    }

  // Fit the piece-wise splines
  std::vector<tk::spline> splines;
  splines = ransac_piece_continuous_splines(10,10,threshold_beta, radii_segments_post_merge_splitting,focal_lengths_segments_post_merge_splitting);
  
  if(!splines.empty()){
    piece_splines_ = splines;
}
  std::cout << "piece_splines_ size: " << piece_splines_.size() << std::endl;
  if(!uncalibrated_areas.empty()){
    uncalibrated_areas.clear();}
    double virtual_grid_length = 10; // Length of the virtual grid
    int num_grids = (new_radii.back() - new_radii.front())/virtual_grid_length;
    double uncalib_size_threshold = 0.3 * new_radii.size()/num_grids;
    bool prev_grid_uncalibrated = false;
    std::vector<double> uncalib_area = {};
    for (int i = 0; i < num_grids; i++) {
      int grid_start = new_radii.front() + i*virtual_grid_length;
      int grid_end = new_radii.front() + (i+1)*virtual_grid_length;
      int grid_size = 0;
      for (int j = 0; j < new_radii.size(); j++) {
        if(new_radii[j] >= grid_start && new_radii[j] < grid_end){
          grid_size += 1;
        }
      }
      if(grid_size < uncalib_size_threshold){
        if(prev_grid_uncalibrated){
          uncalib_area.push_back(grid_end);
        }else{
          uncalib_area.push_back(grid_start);
          uncalib_area.push_back(grid_end);
        }
        prev_grid_uncalibrated = true;
      }else{
        if(prev_grid_uncalibrated){
          uncalibrated_areas.push_back(uncalib_area);
          uncalib_area.clear();
        }
        prev_grid_uncalibrated = false;
      }
    }
    std::cout << "uncalibrated_areas size: " << uncalibrated_areas.size() << std::endl;
    for (int i = 0; i < uncalibrated_areas.size(); i++) {
      std::cout << "uncalibrated_area: " << uncalibrated_areas[i][0] << " " << uncalibrated_areas[i].back() << std::endl;
    }
}

inline void Camera::FitPieceWiseSpline_binary(std::vector<double>& radii, std::vector<double>& focal_lengths) const{
  assert(radii.size() == focal_lengths.size());

  // Ensure that the input data is sorted in increasing order of radii
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

  for (int i = 0; i < increasing_indices.size(); i++) {
    new_radii.push_back(radii[increasing_indices[i]]);
    new_focal_lengths.push_back(focal_lengths[increasing_indices[i]]);
  }
  // double mean_interval = (new_radii.back() - new_radii.front())/new_radii.size();
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
  std::cout << "mean_interval: " << mean_interval << std::endl;
  std::cout << "std_interval: " << std_interval << std::endl;
  std::vector<std::vector<double>> radii_segments = {};
  std::vector<std::vector<double>> focal_lengths_segments = {};
  double threshold = mean_interval + 0.3*std_interval;

  double std_threshold = 0.5*std_interval;  
  recursiveSplit(new_radii, new_focal_lengths, radii_segments, focal_lengths_segments, threshold, std_threshold);
  std::cout << "num_segment: " << radii_segments.size() << std::endl;
  for(int i = 0; i < radii_segments.size(); i++){
    std::cout << "segment_size: " << radii_segments[i].size() << std::endl;
  }
  // identify the largest segment and print the front and back
  std::vector<double> segment_lengths = {};
  for (int i = 0; i < radii_segments.size(); i++) {
    double length = radii_segments[i].back() - radii_segments[i].front();
    segment_lengths.push_back(length);
  }
  int longest_segment_idx = std::distance(segment_lengths.begin(), std::max_element(segment_lengths.begin(), segment_lengths.end()));
  // std::cout << "longest_segment_idx: " << longest_segment_idx << std::endl;
  // std::cout << "longest_segment_size: " << radii_segments[longest_segment_idx].size() << std::endl;
  std::cout << "longest_segment_length: " << segment_lengths[longest_segment_idx] << std::endl;
  std::cout << "longest_segment_front: " << radii_segments[longest_segment_idx].front() << std::endl;
  std::cout << "longest_segment_back: " << radii_segments[longest_segment_idx].back() << std::endl;
  std::vector<std::vector<double>> calibrated_segments = {};
  tk::spline best_spline = ransac_spline_binary(40, 10, 10.0, radii_segments[longest_segment_idx], focal_lengths_segments[longest_segment_idx]);
  spline_ = best_spline;
  intervals_ = {{radii_segments[longest_segment_idx].front(), radii_segments[longest_segment_idx].back()}};
}
int main() {
    // Example data for testing
    std::vector<double> radii;
    std::vector<double> focal_lengths;
    // read the radii and focal_lengths from a txt file
    std::ifstream infile("/home/ivonne/radialsfm/focal_lengths/progressive_courtyard/test.txt");
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        double radius, focal_length;
        if (!(iss >> radius >> focal_length)) {
            std::cerr << "Error reading line: " << line << std::endl;
            continue; // Skip invalid lines
        }
        radii.push_back(radius);
        focal_lengths.push_back(focal_length);
    }

    infile.close();
    std::cout << "radii size: " << radii.size() << std::endl;
    
    Camera camera;
    camera.FitPieceWiseSpline_2(radii, focal_lengths);
    camera.FitGridSpline(radii, focal_lengths);
    camera.FitPieceWiseSpline_binary(radii, focal_lengths);
    std::cout << "Piece-wise spline fitted" << std::endl;
    // creating a uniform grid of radii
    std::vector<double> std_points;
    std_points.clear();
    for (double i = radii[0]; i < radii.back(); i += 10) {
        std_points.push_back(i);
    }
    std::cout << "std_points size: " << std_points.size() << std::endl;
    std::vector<double> focal_lengths_est;
    focal_lengths_est.clear();
    std::vector<double> focal_lengths_whole;
    focal_lengths_whole.clear();
    std::vector<double> focal_lengths_grid;
    focal_lengths_grid.clear();
    for (size_t i = 0; i < std_points.size(); ++i) {
        // std::cout << "std_points: " << std_points[i] << std::endl;
        // focal_lengths_est.push_back(camera.EvalPieceFocalLength(std_points[i]));
        focal_lengths_est.push_back(camera.spline_(std_points[i]));
        focal_lengths_whole.push_back(camera.spline_(std_points[i]));
        // std::cout << "focal_lengths_est: " << focal_lengths_est[i] << std::endl;
        focal_lengths_grid.push_back(camera.EvalGridFocalLength(std_points[i]));
    }
    
    std::ofstream outfile("focal_lengths_est.txt");
    for (size_t i = 0; i < std_points.size(); ++i) {
        outfile << std_points[i] << " " << focal_lengths_est[i] << std::endl;
    }
    outfile.close();
    std::ofstream outfile_whole("focal_lengths_whole.txt");
    for (size_t i = 0; i < std_points.size(); ++i) {
        outfile_whole << std_points[i] << " " << focal_lengths_whole[i] << std::endl;
    }
    outfile_whole.close();
    std::ofstream outfile_grid("focal_lengths_grid.txt");
    for (size_t i = 0; i < std_points.size(); ++i) {
        outfile_grid << std_points[i] << " " << focal_lengths_grid[i] << std::endl;
    }
    outfile_grid.close();
}

#endif