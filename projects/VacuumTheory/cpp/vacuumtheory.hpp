#ifndef VAKYUME_VACUUMTHEORY_HPP
#define VAKYUME_VACUUMTHEORY_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> VacuumTheory_eqn_1_10__P_1(double P_2, double T_1,
                                               double T_2, double V_1,
                                               double V_2) {
  std::vector<double> result;
  double P_1 = (((P_2 * T_1) * V_2) / (T_2 * V_1));
  result.push_back(P_1);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_10__P_2(double P_1, double T_1,
                                               double T_2, double V_1,
                                               double V_2) {
  std::vector<double> result;
  double P_2 = (((P_1 * T_2) * V_1) / (T_1 * V_2));
  result.push_back(P_2);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_10__T_1(double P_1, double P_2,
                                               double T_2, double V_1,
                                               double V_2) {
  std::vector<double> result;
  double T_1 = (((P_1 * T_2) * V_1) / (P_2 * V_2));
  result.push_back(T_1);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_10__T_2(double P_1, double P_2,
                                               double T_1, double V_1,
                                               double V_2) {
  std::vector<double> result;
  double T_2 = (((P_2 * T_1) * V_2) / (P_1 * V_1));
  result.push_back(T_2);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_10__V_1(double P_1, double P_2,
                                               double T_1, double T_2,
                                               double V_2) {
  std::vector<double> result;
  double V_1 = (((P_2 * T_1) * V_2) / (P_1 * T_2));
  result.push_back(V_1);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_10__V_2(double P_1, double P_2,
                                               double T_1, double T_2,
                                               double V_1) {
  std::vector<double> result;
  double V_2 = (((P_1 * T_2) * V_1) / (P_2 * T_1));
  result.push_back(V_2);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_11__M(double P, double T, double W,
                                             double q) {
  std::vector<double> result;
  double M = (((6821.0 * T) * W) / ((738.0 * P) * q));
  result.push_back(M);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_11__P(double M, double T, double W,
                                             double q) {
  std::vector<double> result;
  double P = (((6821.0 * T) * W) / ((738.0 * M) * q));
  result.push_back(P);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_11__T(double M, double P, double W,
                                             double q) {
  std::vector<double> result;
  double T = ((((738.0 * M) * P) * q) / (6821.0 * W));
  result.push_back(T);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_11__W(double M, double P, double T,
                                             double q) {
  std::vector<double> result;
  double W = ((((738.0 * M) * P) * q) / (6821.0 * T));
  result.push_back(W);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_11__q(double M, double P, double T,
                                             double W) {
  std::vector<double> result;
  double q = (((6821.0 * T) * W) / ((738.0 * M) * P));
  result.push_back(q);
  return result;
}

std::vector<double>
VacuumTheory_eqn_1_12__Total_P(double sum_partial_pressures) {
  std::vector<double> result;
  double Total_P = sum_partial_pressures;
  result.push_back(Total_P);
  return result;
}

std::vector<double>
VacuumTheory_eqn_1_12__sum_partial_pressures(double Total_P) {
  std::vector<double> result;
  double sum_partial_pressures = Total_P;
  result.push_back(sum_partial_pressures);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_13a__n(double n_a, double y_a) {
  std::vector<double> result;
  double n = (n_a / y_a);
  result.push_back(n);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_13a__n_a(double n, double y_a) {
  std::vector<double> result;
  double n_a = (n * y_a);
  result.push_back(n_a);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_13a__y_a(double n, double n_a) {
  std::vector<double> result;
  double y_a = (n_a / n);
  result.push_back(y_a);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_13b__P(double p_a, double y_a) {
  std::vector<double> result;
  double P = (p_a / y_a);
  result.push_back(P);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_13b__p_a(double P, double y_a) {
  std::vector<double> result;
  double p_a = (P * y_a);
  result.push_back(p_a);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_13b__y_a(double P, double p_a) {
  std::vector<double> result;
  double y_a = (p_a / P);
  result.push_back(y_a);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_3__T(double k, double m, double v) {
  std::vector<double> result;
  double T = (((0.333333333333333 * m) * std::pow(v, 2.0)) / k);
  result.push_back(T);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_3__k(double T, double m, double v) {
  std::vector<double> result;
  double k = (((0.333333333333333 * m) * std::pow(v, 2.0)) / T);
  result.push_back(k);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_3__m(double T, double k, double v) {
  std::vector<double> result;
  double m = (((3.0 * T) * k) / std::pow(v, 2.0));
  result.push_back(m);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_3__v(double T, double k, double m) {
  std::vector<double> result;
  double v = ((-1.73205080756888) * std::sqrt(((T * k) / m)));
  result.push_back(v);
  v = (1.73205080756888 * std::sqrt(((T * k) / m)));
  result.push_back(v);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_7__R(double T, double V, double n,
                                            double p) {
  std::vector<double> result;
  double R = ((V * p) / (T * n));
  result.push_back(R);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_7__T(double R, double V, double n,
                                            double p) {
  std::vector<double> result;
  double T = ((V * p) / (R * n));
  result.push_back(T);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_7__V(double R, double T, double n,
                                            double p) {
  std::vector<double> result;
  double V = (((R * T) * n) / p);
  result.push_back(V);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_7__n(double R, double T, double V,
                                            double p) {
  std::vector<double> result;
  double n = ((V * p) / (R * T));
  result.push_back(n);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_7__p(double R, double T, double V,
                                            double n) {
  std::vector<double> result;
  double p = (((R * T) * n) / V);
  result.push_back(p);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_8__M(double P, double R, double T,
                                            double V, double m) {
  std::vector<double> result;
  double M = (((R * T) * m) / (P * V));
  result.push_back(M);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_8__P(double M, double R, double T,
                                            double V, double m) {
  std::vector<double> result;
  double P = (((R * T) * m) / (M * V));
  result.push_back(P);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_8__R(double M, double P, double T,
                                            double V, double m) {
  std::vector<double> result;
  double R = (((M * P) * V) / (T * m));
  result.push_back(R);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_8__T(double M, double P, double R,
                                            double V, double m) {
  std::vector<double> result;
  double T = (((M * P) * V) / (R * m));
  result.push_back(T);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_8__V(double M, double P, double R,
                                            double T, double m) {
  std::vector<double> result;
  double V = (((R * T) * m) / (M * P));
  result.push_back(V);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_9__M(double P, double R, double T,
                                            double rho) {
  std::vector<double> result;
  double M = (((R * T) * rho) / P);
  result.push_back(M);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_9__P(double M, double R, double T,
                                            double rho) {
  std::vector<double> result;
  double P = (((R * T) * rho) / M);
  result.push_back(P);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_9__R(double M, double P, double T,
                                            double rho) {
  std::vector<double> result;
  double R = ((M * P) / (T * rho));
  result.push_back(R);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_9__T(double M, double P, double R,
                                            double rho) {
  std::vector<double> result;
  double T = ((M * P) / (R * rho));
  result.push_back(T);
  return result;
}

std::vector<double> VacuumTheory_eqn_1_9__rho(double M, double P, double R,
                                              double T) {
  std::vector<double> result;
  double rho = ((M * P) / (R * T));
  result.push_back(rho);
  return result;
}

#endif // VAKYUME_VACUUMTHEORY_HPP
