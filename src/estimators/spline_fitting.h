
#ifndef COLMAP_SRC_ESTIMATORS_SPLINE_FITTING_H_
#define COLMAP_SRC_ESTIMATORS_SPLINE_FITTING_H_

#include <Eigen/Core>
#include "base/spline.h"


#include <random>
#include <vector>
namespace colmap {

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

} // namespace colmap

#endif  // COLMAP_SRC_ESTIMATORS_SPLINE_FITTING_H_