#ifndef VAKYUME_LINEARMOMENTUMANDTHECENTREOFMASS_HPP
#define VAKYUME_LINEARMOMENTUMANDTHECENTREOFMASS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> LinearMomentumAndTheCentreOfMass_eqn_10_1__m(double p,
                                                                 double v) {
  std::vector<double> result;
  double m = (p / v);
  result.push_back(m);
  return result;
}

std::vector<double> LinearMomentumAndTheCentreOfMass_eqn_10_1__p(double m,
                                                                 double v) {
  std::vector<double> result;
  double p = (m * v);
  result.push_back(p);
  return result;
}

std::vector<double> LinearMomentumAndTheCentreOfMass_eqn_10_1__v(double m,
                                                                 double p) {
  std::vector<double> result;
  double v = (p / m);
  result.push_back(v);
  return result;
}

std::vector<double>
LinearMomentumAndTheCentreOfMass_eqn_10_11__F_ext(double M, double a_CM) {
  std::vector<double> result;
  double F_ext = (M * a_CM);
  result.push_back(F_ext);
  return result;
}

std::vector<double> LinearMomentumAndTheCentreOfMass_eqn_10_11__M(double F_ext,
                                                                  double a_CM) {
  std::vector<double> result;
  double M = (F_ext / a_CM);
  result.push_back(M);
  return result;
}

std::vector<double>
LinearMomentumAndTheCentreOfMass_eqn_10_11__a_CM(double F_ext, double M) {
  std::vector<double> result;
  double a_CM = (F_ext / M);
  result.push_back(a_CM);
  return result;
}

#endif // VAKYUME_LINEARMOMENTUMANDTHECENTREOFMASS_HPP
