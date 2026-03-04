#ifndef VAKYUME_THETHEORYOFSPECIALRELATIVITY_HPP
#define VAKYUME_THETHEORYOFSPECIALRELATIVITY_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_11__E(double c, double m_0, double p) {
  std::vector<double> result;
  double E = (-std::sqrt(
      ((std::pow(c, 2.0) * std::pow(p, 2.0)) + std::pow(m_0, 2.0))));
  result.push_back(E);
  E = std::sqrt(((std::pow(c, 2.0) * std::pow(p, 2.0)) + std::pow(m_0, 2.0)));
  result.push_back(E);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_11__c(double E, double m_0, double p) {
  std::vector<double> result;
  double c = ((-std::sqrt(((E - m_0) * (E + m_0)))) / p);
  result.push_back(c);
  c = (std::sqrt(((E - m_0) * (E + m_0))) / p);
  result.push_back(c);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_11__m_0(double E, double c, double p) {
  std::vector<double> result;
  double m_0 = (-std::sqrt(((E - (c * p)) * (E + (c * p)))));
  result.push_back(m_0);
  m_0 = std::sqrt(((E - (c * p)) * (E + (c * p))));
  result.push_back(m_0);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_11__p(double E, double c, double m_0) {
  std::vector<double> result;
  double p = ((-std::sqrt(((E - m_0) * (E + m_0)))) / c);
  result.push_back(p);
  p = (std::sqrt(((E - m_0) * (E + m_0))) / c);
  result.push_back(p);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__FE(double l,
                                                              double r) {
  std::vector<double> result;
  double FE = (l / r);
  result.push_back(FE);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__c(double t, double v,
                                                             double x) {
  throw std::runtime_error("TheTheoryOfSpecialRelativity_eqn_24_2__c: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__l(double FE,
                                                             double r) {
  std::vector<double> result;
  double l = (FE * r);
  result.push_back(l);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__r(double FE,
                                                             double l) {
  std::vector<double> result;
  double r = (l / FE);
  result.push_back(r);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__t(double c, double v,
                                                             double x) {
  throw std::runtime_error("TheTheoryOfSpecialRelativity_eqn_24_2__t: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__u_0(double c,
                                                               double u_x,
                                                               double v,
                                                               double x) {
  std::vector<double> result;
  double u_0 =
      ((u_x * (std::pow(c, 2.0) - (u_x * v))) / (std::pow(c, 2.0) * x));
  result.push_back(u_0);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__u_x(double c,
                                                               double u_0,
                                                               double v,
                                                               double x) {
  std::vector<double> result;
  double u_x =
      ((c * (c - std::sqrt((std::pow(c, 2.0) - (((4.0 * u_0) * v) * x))))) /
       (2.0 * v));
  result.push_back(u_x);
  u_x = ((c * (c + std::sqrt((std::pow(c, 2.0) - (((4.0 * u_0) * v) * x))))) /
         (2.0 * v));
  result.push_back(u_x);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__v(double c, double t,
                                                             double x) {
  throw std::runtime_error("TheTheoryOfSpecialRelativity_eqn_24_2__v: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__vx(double c,
                                                              double t,
                                                              double v,
                                                              double x) {
  std::vector<double> result;
  double vx = (std::pow(c, 2.0) * x);
  result.push_back(vx);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_2__x(double c, double t,
                                                             double v) {
  std::vector<double> result;
  double x = 0.0;
  result.push_back(x);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_6__c(double u, double ux, double v) {
  std::vector<double> result;
  double c = std::sqrt((((u * v) * ((-u) + ux)) / ((u - ux) + v)));
  result.push_back(c);
  c = (-std::sqrt(((((-u) * v) * (u - ux)) / ((u - ux) + v))));
  result.push_back(c);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_6__u(double c, double ux, double v) {
  std::vector<double> result;
  double u = ((((-std::pow(c, 2.0)) + (ux * v)) -
               std::sqrt((((std::pow(c, 2.0) - ((2.0 * c) * v)) + (ux * v)) *
                          ((std::pow(c, 2.0) + ((2.0 * c) * v)) + (ux * v))))) /
              (2.0 * v));
  result.push_back(u);
  u = ((((-std::pow(c, 2.0)) + (ux * v)) +
        std::sqrt((((std::pow(c, 2.0) - ((2.0 * c) * v)) + (ux * v)) *
                   ((std::pow(c, 2.0) + ((2.0 * c) * v)) + (ux * v))))) /
       (2.0 * v));
  result.push_back(u);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_6__ux(double c, double u, double v) {
  std::vector<double> result;
  double ux = ((((std::pow(c, 2.0) * u) + (std::pow(c, 2.0) * v)) +
                (std::pow(u, 2.0) * v)) /
               (std::pow(c, 2.0) + (u * v)));
  result.push_back(ux);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_6__v(double c, double u,
                                                             double ux) {
  std::vector<double> result;
  double v = ((std::pow(c, 2.0) * ((-u) + ux)) /
              ((std::pow(c, 2.0) + std::pow(u, 2.0)) - (u * ux)));
  result.push_back(v);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_9__K(double c, double m_0, double u) {
  std::vector<double> result;
  double K = (m_0 * (std::pow(c, 2.0) - std::pow(u, 2.0)));
  result.push_back(K);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_9__c(double K, double m_0, double u) {
  std::vector<double> result;
  double c = (-std::sqrt(((K + (m_0 * std::pow(u, 2.0))) / m_0)));
  result.push_back(c);
  c = std::sqrt(((K + (m_0 * std::pow(u, 2.0))) / m_0));
  result.push_back(c);
  return result;
}

std::vector<double>
TheTheoryOfSpecialRelativity_eqn_24_9__m_0(double K, double c, double u) {
  std::vector<double> result;
  double m_0 = (K / (std::pow(c, 2.0) - std::pow(u, 2.0)));
  result.push_back(m_0);
  return result;
}

std::vector<double> TheTheoryOfSpecialRelativity_eqn_24_9__u(double K, double c,
                                                             double m_0) {
  std::vector<double> result;
  double u = (-std::sqrt((((-K) / m_0) + std::pow(c, 2.0))));
  result.push_back(u);
  u = std::sqrt((((-K) / m_0) + std::pow(c, 2.0)));
  result.push_back(u);
  return result;
}

#endif // VAKYUME_THETHEORYOFSPECIALRELATIVITY_HPP
