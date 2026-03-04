#ifndef VAKYUME_AIRLEAK_HPP
#define VAKYUME_AIRLEAK_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> AirLeak_eqn_4_10__T(double V, double del_P, double leakage,
                                        double t) {
  std::vector<double> result;
  double T = (((3.127 * V) * del_P) / (leakage * t));
  result.push_back(T);
  return result;
}

std::vector<double> AirLeak_eqn_4_10__V(double T, double del_P, double leakage,
                                        double t) {
  std::vector<double> result;
  double V = ((((0.319795330988168 * T) * leakage) * t) / del_P);
  result.push_back(V);
  return result;
}

std::vector<double> AirLeak_eqn_4_10__del_P(double T, double V, double leakage,
                                            double t) {
  std::vector<double> result;
  double del_P = ((((0.319795330988168 * T) * leakage) * t) / V);
  result.push_back(del_P);
  return result;
}

std::vector<double> AirLeak_eqn_4_10__leakage(double T, double V, double del_P,
                                              double t) {
  std::vector<double> result;
  double leakage = (((3.127 * V) * del_P) / (T * t));
  result.push_back(leakage);
  return result;
}

std::vector<double> AirLeak_eqn_4_7__W(double W_T,
                                       double sum_individual_leak_rates) {
  std::vector<double> result;
  double W = (W_T - sum_individual_leak_rates);
  result.push_back(W);
  return result;
}

std::vector<double> AirLeak_eqn_4_7__W_T(double W,
                                         double sum_individual_leak_rates) {
  std::vector<double> result;
  double W_T = (W + sum_individual_leak_rates);
  result.push_back(W_T);
  return result;
}

std::vector<double> AirLeak_eqn_4_7__sum_individual_leak_rates(double W,
                                                               double W_T) {
  std::vector<double> result;
  double sum_individual_leak_rates = ((-W) + W_T);
  result.push_back(sum_individual_leak_rates);
  return result;
}

#endif // VAKYUME_AIRLEAK_HPP
