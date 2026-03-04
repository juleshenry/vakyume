#ifndef VAKYUME_FLUIDMECHANICS_HPP
#define VAKYUME_FLUIDMECHANICS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> FluidMechanics_eqn_15_8__R(double R1, double R2) {
  std::vector<double> result;
  double R = (R1 + R2);
  result.push_back(R);
  return result;
}

std::vector<double> FluidMechanics_eqn_15_8__R1(double R, double R2) {
  std::vector<double> result;
  double R1 = (R - R2);
  result.push_back(R1);
  return result;
}

std::vector<double> FluidMechanics_eqn_15_8__R2(double R, double R1) {
  std::vector<double> result;
  double R2 = (R - R1);
  result.push_back(R2);
  return result;
}

#endif // VAKYUME_FLUIDMECHANICS_HPP
