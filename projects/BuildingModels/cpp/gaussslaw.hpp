#ifndef VAKYUME_GAUSSSLAW_HPP
#define VAKYUME_GAUSSSLAW_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> GaussSLaw_eqn_17_2__E(double R) {
  std::vector<double> result;
  double E = 0.0;
  result.push_back(E);
  return result;
}

std::vector<double> GaussSLaw_eqn_17_2__R(double E) {
  std::vector<double> result;
  double R = ((-1.0) / 2.0);
  result.push_back(R);
  R = (1.0 / 2.0);
  result.push_back(R);
  return result;
}

std::vector<double> GaussSLaw_eqn_17_3__E(double Q, double R) {
  std::vector<double> result;
  double E = (Q / (4.0 * std::pow(R, 2.0)));
  result.push_back(E);
  return result;
}

std::vector<double> GaussSLaw_eqn_17_3__Q(double E, double R) {
  std::vector<double> result;
  double Q = ((4.0 * E) * std::pow(R, 2.0));
  result.push_back(Q);
  return result;
}

std::vector<double> GaussSLaw_eqn_17_3__R(double E, double Q) {
  std::vector<double> result;
  double R = ((-std::sqrt((Q / E))) / 2.0);
  result.push_back(R);
  R = (std::sqrt((Q / E)) / 2.0);
  result.push_back(R);
  return result;
}

#endif // VAKYUME_GAUSSSLAW_HPP
