#ifndef VAKYUME_APPLYINGNEWTONSLAWS_HPP
#define VAKYUME_APPLYINGNEWTONSLAWS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> ApplyingNewtonSLaws_eqn_6_1__F(double F_N) {
  std::vector<double> result;
  double F = F_N;
  result.push_back(F);
  return result;
}

std::vector<double> ApplyingNewtonSLaws_eqn_6_1__F_N(double F) {
  std::vector<double> result;
  double F_N = F;
  result.push_back(F_N);
  return result;
}

std::vector<double> ApplyingNewtonSLaws_eqn_6_2__a_2(double m, double x) {
  std::vector<double> result;
  double a_2 = ((-x) / m);
  result.push_back(a_2);
  return result;
}

std::vector<double> ApplyingNewtonSLaws_eqn_6_2__m(double a_2, double x) {
  std::vector<double> result;
  double m = ((-x) / a_2);
  result.push_back(m);
  return result;
}

std::vector<double> ApplyingNewtonSLaws_eqn_6_2__x(double a_2, double m) {
  std::vector<double> result;
  double x = ((-a_2) * m);
  result.push_back(x);
  return result;
}

#endif // VAKYUME_APPLYINGNEWTONSLAWS_HPP
