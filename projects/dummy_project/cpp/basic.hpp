#ifndef VAKYUME_BASIC_HPP
#define VAKYUME_BASIC_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

std::vector<double> Basic_eqn_1_1__x(double y, double z) {
    std::vector<double> result;
    double x = ((-y) + z);
    result.push_back(x);
    return result;
}

std::vector<double> Basic_eqn_1_1__y(double x, double z) {
    std::vector<double> result;
    double y = ((-x) + z);
    result.push_back(y);
    return result;
}

std::vector<double> Basic_eqn_1_1__z(double x, double y) {
    std::vector<double> result;
    double z = (x + y);
    result.push_back(z);
    return result;
}

#endif // VAKYUME_BASIC_HPP
