#ifndef VAKYUME_ELECTRICCURRENT_HPP
#define VAKYUME_ELECTRICCURRENT_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<std::complex<double>> ElectricCurrent_eqn_19_1__I(double Q,
                                                              double t) {
  std::vector<std::complex<double>> result;
  std::complex<double> I = (Q / t);
  result.push_back(std::complex<double>(0.0, 1.0));
  return result;
}

std::vector<std::complex<double>> ElectricCurrent_eqn_19_1__Q(double I,
                                                              double t) {
  std::vector<std::complex<double>> result;
  std::complex<double> Q = (std::complex<double>(0.0, 1.0) * t);
  result.push_back(Q);
  return result;
}

std::vector<std::complex<double>> ElectricCurrent_eqn_19_1__t(double I,
                                                              double Q) {
  std::vector<std::complex<double>> result;
  std::complex<double> t = (Q / std::complex<double>(0.0, 1.0));
  result.push_back(t);
  return result;
}

std::vector<double> ElectricCurrent_eqn_19_6__R1(double R2, double Reff) {
  throw std::runtime_error("ElectricCurrent_eqn_19_6__R1: requires numerical "
                           "solver (not transpilable)");
}

std::vector<double> ElectricCurrent_eqn_19_6__R2(double R1, double Reff) {
  throw std::runtime_error("ElectricCurrent_eqn_19_6__R2: requires numerical "
                           "solver (not transpilable)");
}

std::vector<double> ElectricCurrent_eqn_19_6__Reff(double R1, double R2) {
  throw std::runtime_error("ElectricCurrent_eqn_19_6__Reff: requires numerical "
                           "solver (not transpilable)");
}

#endif // VAKYUME_ELECTRICCURRENT_HPP
