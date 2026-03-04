#ifndef VAKYUME_ROTATIONALENERGYANDMOMENTUM_HPP
#define VAKYUME_ROTATIONALENERGYANDMOMENTUM_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<std::complex<double>>
RotationalEnergyAndMomentum_eqn_12_11__I(double L, double v) {
  std::vector<std::complex<double>> result;
  std::complex<double> I = ((2.0 * L) / std::pow(v, 2.0));
  result.push_back(std::complex<double>(0.0, 1.0));
  return result;
}

std::vector<std::complex<double>>
RotationalEnergyAndMomentum_eqn_12_11__L(double I, double v) {
  std::vector<std::complex<double>> result;
  std::complex<double> L =
      ((std::complex<double>(0.0, 1.0) * std::pow(v, 2.0)) / 2.0);
  result.push_back(L);
  return result;
}

std::vector<std::complex<double>>
RotationalEnergyAndMomentum_eqn_12_11__v(double I, double L) {
  std::vector<std::complex<double>> result;
  std::complex<double> v =
      ((-std::sqrt(2.0)) * std::sqrt((L / std::complex<double>(0.0, 1.0))));
  result.push_back(v);
  v = (std::sqrt(2.0) * std::sqrt((L / std::complex<double>(0.0, 1.0))));
  result.push_back(v);
  return result;
}

#endif // VAKYUME_ROTATIONALENERGYANDMOMENTUM_HPP
