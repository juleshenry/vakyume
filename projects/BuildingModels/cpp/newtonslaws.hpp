#ifndef VAKYUME_NEWTONSLAWS_HPP
#define VAKYUME_NEWTONSLAWS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> NewtonSLaws_eqn_5_1__F(double a, double m) {
  std::vector<double> result;
  double F = (a * m);
  result.push_back(F);
  return result;
}

std::vector<double> NewtonSLaws_eqn_5_1__a(double F, double m) {
  std::vector<double> result;
  double a = (F / m);
  result.push_back(a);
  return result;
}

std::vector<double> NewtonSLaws_eqn_5_1__m(double F, double a) {
  std::vector<double> result;
  double m = (F / a);
  result.push_back(m);
  return result;
}

std::vector<double> NewtonSLaws_eqn_5_2__F(double x) {
  std::vector<double> result;
  return {std::sqrt(x)};
}

std::vector<double> NewtonSLaws_eqn_5_2__F_action(double F_reaction) {
  std::vector<double> result;
  double F_action = F_reaction;
  result.push_back(F_action);
  return result;
}

std::vector<double> NewtonSLaws_eqn_5_2__F_reaction(double F_action) {
  std::vector<double> result;
  double F_reaction = F_action;
  result.push_back(F_reaction);
  return result;
}

std::vector<double> NewtonSLaws_eqn_5_2__x(double F) {
  std::vector<double> result;
  return {std::sqrt(F)};
}

#endif // VAKYUME_NEWTONSLAWS_HPP
