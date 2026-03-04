#ifndef VAKYUME_WAVES_HPP
#define VAKYUME_WAVES_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> Waves_eqn_14_11__fn(double n, double v) {
  std::vector<double> result;
  return {((n * v) / (2.0 * std::sqrt(2.0)))};
}

std::vector<double> Waves_eqn_14_11__n(double fn, double v) {
  std::vector<double> result;
  return {std::sqrt(((2.0 * v) / fn))};
}

std::vector<double> Waves_eqn_14_11__v(double fn, double n) {
  std::vector<double> result;
  return {std::sqrt(((2.0 * fn) * n))};
}

#endif // VAKYUME_WAVES_HPP
