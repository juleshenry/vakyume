#ifndef VAKYUME_ROTATIONALDYNAMICS_HPP
#define VAKYUME_ROTATIONALDYNAMICS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> RotationalDynamics_eqn_11_10__ICM(double Ih, double M,
                                                      double h) {
  std::vector<double> result;
  double ICM = (Ih - (M * std::pow(h, 2.0)));
  result.push_back(ICM);
  return result;
}

std::vector<double> RotationalDynamics_eqn_11_10__Ih(double ICM, double M,
                                                     double h) {
  std::vector<double> result;
  double Ih = (ICM + (M * std::pow(h, 2.0)));
  result.push_back(Ih);
  return result;
}

std::vector<double> RotationalDynamics_eqn_11_10__M(double ICM, double Ih,
                                                    double h) {
  std::vector<double> result;
  double M = (((-ICM) + Ih) / std::pow(h, 2.0));
  result.push_back(M);
  return result;
}

std::vector<double> RotationalDynamics_eqn_11_10__h(double ICM, double Ih,
                                                    double M) {
  std::vector<double> result;
  double h = std::sqrt((((-ICM) + Ih) / M));
  result.push_back(h);
  h = (-std::sqrt(((-(ICM - Ih)) / M)));
  result.push_back(h);
  return result;
}

std::vector<std::complex<double>>
RotationalDynamics_eqn_11_8__I(double i, double m, double r) {
  std::vector<std::complex<double>> result;
  std::complex<double> I = ((m * std::pow(r, 2.0)) / i);
  result.push_back(std::complex<double>(0.0, 1.0));
  return result;
}

std::vector<std::complex<double>>
RotationalDynamics_eqn_11_8__i(double I, double m, double r) {
  std::vector<std::complex<double>> result;
  std::complex<double> i =
      ((m * std::pow(r, 2.0)) / std::complex<double>(0.0, 1.0));
  result.push_back(i);
  return result;
}

std::vector<std::complex<double>>
RotationalDynamics_eqn_11_8__m(double I, double i, double r) {
  std::vector<std::complex<double>> result;
  std::complex<double> m =
      ((std::complex<double>(0.0, 1.0) * i) / std::pow(r, 2.0));
  result.push_back(m);
  return result;
}

std::vector<std::complex<double>>
RotationalDynamics_eqn_11_8__r(double I, double i, double m) {
  std::vector<std::complex<double>> result;
  std::complex<double> r =
      (-std::sqrt(((std::complex<double>(0.0, 1.0) * i) / m)));
  result.push_back(r);
  r = std::sqrt(((std::complex<double>(0.0, 1.0) * i) / m));
  result.push_back(r);
  return result;
}

#endif // VAKYUME_ROTATIONALDYNAMICS_HPP
