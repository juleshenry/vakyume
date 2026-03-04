#ifndef VAKYUME_SIMPLEHARMONICMOTION_HPP
#define VAKYUME_SIMPLEHARMONICMOTION_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> SimpleHarmonicMotion_eqn_13_1__E(double m, double v,
                                                     double x) {
  std::vector<double> result;
  double E =
      (((0.5 * m) * std::pow(v, 2.0)) + std::pow(2.0, (-std::pow(x, 2.0))));
  result.push_back(E);
  return result;
}

std::vector<double> SimpleHarmonicMotion_eqn_13_1__m(double E, double v,
                                                     double x) {
  std::vector<double> result;
  double m = (((2.0 * E) / std::pow(v, 2.0)) -
              (2.0 / (std::pow(2.0, std::pow(x, 2.0)) * std::pow(v, 2.0))));
  result.push_back(m);
  return result;
}

std::vector<double> SimpleHarmonicMotion_eqn_13_1__v(double E, double m,
                                                     double x) {
  std::vector<double> result;
  double v = ((-1.4142135623731) *
              std::sqrt(((E - (1.0 / std::pow(2.0, std::pow(x, 2.0)))) / m)));
  result.push_back(v);
  v = (1.4142135623731 *
       std::sqrt(((E - (1.0 / std::pow(2.0, std::pow(x, 2.0)))) / m)));
  result.push_back(v);
  return result;
}

std::vector<double> SimpleHarmonicMotion_eqn_13_1__x(double E, double m,
                                                     double v) {
  std::vector<double> result;
  double x =
      ((-1.20112240878645) *
       std::sqrt(std::log((2.0 / ((2.0 * E) - (m * std::pow(v, 2.0)))))));
  result.push_back(x);
  x = (1.20112240878645 *
       std::sqrt(std::log((2.0 / ((2.0 * E) - (m * std::pow(v, 2.0)))))));
  result.push_back(x);
  return result;
}

std::vector<double> SimpleHarmonicMotion_eqn_13_2__d2x(double dt2, double m,
                                                       double x) {
  throw std::runtime_error("SimpleHarmonicMotion_eqn_13_2__d2x: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> SimpleHarmonicMotion_eqn_13_2__dt2(double d2x, double m,
                                                       double x) {
  throw std::runtime_error("SimpleHarmonicMotion_eqn_13_2__dt2: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> SimpleHarmonicMotion_eqn_13_2__m(double d2x, double dt2,
                                                     double x) {
  throw std::runtime_error("SimpleHarmonicMotion_eqn_13_2__m: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> SimpleHarmonicMotion_eqn_13_2__x(double d2x, double dt2,
                                                     double m) {
  throw std::runtime_error("SimpleHarmonicMotion_eqn_13_2__x: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> SimpleHarmonicMotion_eqn_13_2__x0(double x1, double x2) {
  std::vector<double> result;
  double x0 = ((x1 / 3.0) + ((2.0 * x2) / 3.0));
  result.push_back(x0);
  return result;
}

std::vector<double> SimpleHarmonicMotion_eqn_13_2__x1(double x0, double x2) {
  std::vector<double> result;
  double x1 = ((3.0 * x0) - (2.0 * x2));
  result.push_back(x1);
  return result;
}

std::vector<double> SimpleHarmonicMotion_eqn_13_2__x2(double x0, double x1) {
  std::vector<double> result;
  double x2 = (((3.0 * x0) / 2.0) - (x1 / 2.0));
  result.push_back(x2);
  return result;
}

#endif // VAKYUME_SIMPLEHARMONICMOTION_HPP
