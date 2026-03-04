#ifndef VAKYUME_ELECTROMAGNETICINDUCTION_HPP
#define VAKYUME_ELECTROMAGNETICINDUCTION_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> ElectromagneticInduction_eqn_23_1__B0(double E, double a,
                                                          double r) {
  std::vector<double> result;
  double B0 = (E / (a * std::pow(r, 2.0)));
  result.push_back(B0);
  return result;
}

std::vector<double> ElectromagneticInduction_eqn_23_1__E(double B0, double a,
                                                         double r) {
  std::vector<double> result;
  double E = ((B0 * a) * std::pow(r, 2.0));
  result.push_back(E);
  return result;
}

std::vector<double> ElectromagneticInduction_eqn_23_1__R(double B0, double E,
                                                         double a) {
  std::vector<double> result;
  double R = (-std::sqrt((E / (B0 * a))));
  result.push_back(R);
  R = std::sqrt((E / (B0 * a)));
  result.push_back(R);
  return result;
}

std::vector<double> ElectromagneticInduction_eqn_23_1__V(double dt) {
  throw std::runtime_error("ElectromagneticInduction_eqn_23_1__V: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> ElectromagneticInduction_eqn_23_1__a(double B0, double E,
                                                         double r) {
  std::vector<double> result;
  double a = (E / (B0 * std::pow(r, 2.0)));
  result.push_back(a);
  return result;
}

std::vector<double> ElectromagneticInduction_eqn_23_1__dt(double V) {
  throw std::runtime_error("ElectromagneticInduction_eqn_23_1__dt: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> ElectromagneticInduction_eqn_23_1__r(double B0, double E,
                                                         double a) {
  std::vector<double> result;
  double r = (-std::sqrt((E / (B0 * a))));
  result.push_back(r);
  r = std::sqrt((E / (B0 * a)));
  result.push_back(r);
  return result;
}

#endif // VAKYUME_ELECTROMAGNETICINDUCTION_HPP
