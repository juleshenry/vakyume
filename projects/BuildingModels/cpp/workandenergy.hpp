#ifndef VAKYUME_WORKANDENERGY_HPP
#define VAKYUME_WORKANDENERGY_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> WorkAndEnergy_eqn_7_1__F(double W, double d) {
  std::vector<double> result;
  double F = (W / d);
  result.push_back(F);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_1__W(double F, double d) {
  std::vector<double> result;
  double W = (F * d);
  result.push_back(W);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_1__d(double F, double W) {
  std::vector<double> result;
  double d = (W / F);
  result.push_back(d);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_2__F1(double F2, double F3,
                                              double W_tot, double x) {
  std::vector<double> result;
  double F1 = (((-F2) - F3) + (W_tot / x));
  result.push_back(F1);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_2__F2(double F1, double F3,
                                              double W_tot, double x) {
  std::vector<double> result;
  double F2 = (((-F1) - F3) + (W_tot / x));
  result.push_back(F2);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_2__F3(double F1, double F2,
                                              double W_tot, double x) {
  std::vector<double> result;
  double F3 = (((-F1) - F2) + (W_tot / x));
  result.push_back(F3);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_2__W_tot(double F1, double F2,
                                                 double F3, double x) {
  std::vector<double> result;
  double W_tot = (x * ((F1 + F2) + F3));
  result.push_back(W_tot);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_2__x(double F1, double F2, double F3,
                                             double W_tot) {
  std::vector<double> result;
  double x = (W_tot / ((F1 + F2) + F3));
  result.push_back(x);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_6__F_(double W_f, double d) {
  std::vector<double> result;
  double F_ = (W_f / d);
  result.push_back(F_);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_6__W_f(double F_, double d) {
  std::vector<double> result;
  double W_f = (F_ * d);
  result.push_back(W_f);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_6__d(double F_, double W_f) {
  std::vector<double> result;
  double d = (W_f / F_);
  result.push_back(d);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_7__W_F(double W_f, double W_g,
                                               double W_net) {
  std::vector<double> result;
  double W_F = (((-W_f) - W_g) + W_net);
  result.push_back(W_F);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_7__W_f(double W_F, double W_g,
                                               double W_net) {
  std::vector<double> result;
  double W_f = (((-W_F) - W_g) + W_net);
  result.push_back(W_f);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_7__W_g(double W_F, double W_f,
                                               double W_net) {
  std::vector<double> result;
  double W_g = (((-W_F) - W_f) + W_net);
  result.push_back(W_g);
  return result;
}

std::vector<double> WorkAndEnergy_eqn_7_7__W_net(double W_F, double W_f,
                                                 double W_g) {
  std::vector<double> result;
  double W_net = ((W_F + W_f) + W_g);
  result.push_back(W_net);
  return result;
}

#endif // VAKYUME_WORKANDENERGY_HPP
