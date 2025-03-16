#ifndef TEST_SPLINE_SPEED_CC_
#define TEST_SPLINE_SPEED_CC_
#include <iostream>
#include <vector>
#include <cassert>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>

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
int main() {
    std::vector<double> sample_x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> sample_y = {1, 4, 9, 16, 25, 36, 49, 64, 81, 100};
    tk::spline<double> spline_focal_lengths;
    spline_focal_lengths.set_points(sample_x, sample_y);
    std::vector<double> radii = linspace(1, 10, 1000);
    std::vector<double> focal_lengths;
    for (int i = 0; i < radii.size(); i++) {
        focal_lengths.push_back(spline_focal_lengths(radii[i]));
    }
    std::ofstream file("tk_spline_value.txt");
    for (int i = 0; i < radii.size(); i++) {
        file << radii[i] << " " << focal_lengths[i] << std::endl;
    }


}

#endif