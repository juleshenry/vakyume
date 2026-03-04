#ifndef VAKYUME_THEMAGNETICFORCE_HPP
#define VAKYUME_THEMAGNETICFORCE_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> TheMagneticForce_eqn_21_4__B(double T, double m, double q,
                                                 double v) {
  std::vector<double> result;
  double B = ((q * v) / (2.0 * T));
  result.push_back(B);
  return result;
}

std::vector<double> TheMagneticForce_eqn_21_4__T(double B, double m, double q,
                                                 double v) {
  std::vector<double> result;
  double T = ((q * v) / (2.0 * B));
  result.push_back(T);
  return result;
}

std::vector<double> TheMagneticForce_eqn_21_4__m(double B, double T, double q,
                                                 double v) {
  throw std::runtime_error("TheMagneticForce_eqn_21_4__m: requires numerical "
                           "solver (not transpilable)");
}

std::vector<double> TheMagneticForce_eqn_21_4__q(double B, double T, double m,
                                                 double v) {
  std::vector<double> result;
  double q = (((2.0 * B) * T) / v);
  result.push_back(q);
  return result;
}

std::vector<double> TheMagneticForce_eqn_21_4__v(double B, double T, double m,
                                                 double q) {
  std::vector<double> result;
  double v = (((2.0 * B) * T) / q);
  result.push_back(v);
  return result;
}

#endif // VAKYUME_THEMAGNETICFORCE_HPP
