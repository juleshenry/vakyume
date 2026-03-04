#ifndef VAKYUME_ELECTRICCHARGESANDFIELDS_HPP
#define VAKYUME_ELECTRICCHARGESANDFIELDS_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> ElectricChargesAndFields_eqn_16_3__F_g(double G, double m_e,
                                                           double m_p,
                                                           double r) {
  std::vector<double> result;
  double F_g = (((G * m_e) * m_p) / std::pow(r, 2.0));
  result.push_back(F_g);
  return result;
}

std::vector<double> ElectricChargesAndFields_eqn_16_3__G(double F_g, double m_e,
                                                         double m_p, double r) {
  std::vector<double> result;
  double G = ((F_g * std::pow(r, 2.0)) / (m_e * m_p));
  result.push_back(G);
  return result;
}

std::vector<double> ElectricChargesAndFields_eqn_16_3__m_e(double F_g, double G,
                                                           double m_p,
                                                           double r) {
  std::vector<double> result;
  double m_e = ((F_g * std::pow(r, 2.0)) / (G * m_p));
  result.push_back(m_e);
  return result;
}

std::vector<double> ElectricChargesAndFields_eqn_16_3__m_p(double F_g, double G,
                                                           double m_e,
                                                           double r) {
  std::vector<double> result;
  double m_p = ((F_g * std::pow(r, 2.0)) / (G * m_e));
  result.push_back(m_p);
  return result;
}

std::vector<double> ElectricChargesAndFields_eqn_16_3__r(double F_g, double G,
                                                         double m_e,
                                                         double m_p) {
  std::vector<double> result;
  double r = (-std::sqrt((((G * m_e) * m_p) / F_g)));
  result.push_back(r);
  r = std::sqrt((((G * m_e) * m_p) / F_g));
  result.push_back(r);
  return result;
}

#endif // VAKYUME_ELECTRICCHARGESANDFIELDS_HPP
