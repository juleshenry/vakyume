#ifndef VAKYUME_PROCESSAPP1_HPP
#define VAKYUME_PROCESSAPP1_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> ProcessApp1_eqn_5_1__K_i(double x_i, double y_i) {
  std::vector<double> result;
  double K_i = (y_i / x_i);
  result.push_back(K_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_1__x_i(double K_i, double y_i) {
  std::vector<double> result;
  double x_i = (y_i / K_i);
  result.push_back(x_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_1__y_i(double K_i, double x_i) {
  std::vector<double> result;
  double y_i = (K_i * x_i);
  result.push_back(y_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10a__D(double L_0, double V_1) {
  std::vector<double> result;
  double D = ((-L_0) + V_1);
  result.push_back(D);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10a__L_0(double D, double V_1) {
  std::vector<double> result;
  double L_0 = 0.0;
  result.push_back(L_0);
  L_0 = ((-D) + V_1);
  result.push_back(L_0);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10a__V_1(double D, double L_0) {
  std::vector<double> result;
  double V_1 = (D + L_0);
  result.push_back(V_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10b__L_0(double R, double V_1) {
  std::vector<double> result;
  double L_0 = ((R * V_1) / (R + 1.0));
  result.push_back(L_0);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10b__R(double L_0, double V_1) {
  std::vector<double> result;
  double R = ((-L_0) / (L_0 - V_1));
  result.push_back(R);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10b__V_1(double L_0, double R) {
  std::vector<double> result;
  double V_1 = (L_0 + (L_0 / R));
  result.push_back(V_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10c__D(double L_0, double R) {
  std::vector<double> result;
  double D = (L_0 / R);
  result.push_back(D);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10c__L_0(double D, double R) {
  std::vector<double> result;
  double L_0 = (D * R);
  result.push_back(L_0);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_10c__R(double D, double L_0) {
  std::vector<double> result;
  double R = (L_0 / D);
  result.push_back(R);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_11__B(double L_N, double V_0) {
  std::vector<double> result;
  double B = (L_N - V_0);
  result.push_back(B);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_11__L_N(double B, double V_0) {
  std::vector<double> result;
  double L_N = (B + V_0);
  result.push_back(L_N);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_11__V_0(double B, double L_N) {
  std::vector<double> result;
  double V_0 = ((-B) + L_N);
  result.push_back(V_0);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_12__Eff(double N_ES, double N_t,
                                              double T) {
  std::vector<double> result;
  double Eff = std::pow((N_ES / N_t), (1.0 / T));
  result.push_back(Eff);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_12__N_ES(double Eff, double N_t,
                                               double T) {
  std::vector<double> result;
  double N_ES = (std::pow(Eff, T) * N_t);
  result.push_back(N_ES);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_12__N_t(double Eff, double N_ES,
                                              double T) {
  std::vector<double> result;
  double N_t = (N_ES / std::pow(Eff, T));
  result.push_back(N_t);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_12__T(double Eff, double N_ES,
                                            double N_t) {
  std::vector<double> result;
  double T = (std::log((N_ES / N_t)) / std::log(Eff));
  result.push_back(T);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_13__HETP(double H_p, double N_ES) {
  std::vector<double> result;
  double HETP = (H_p / N_ES);
  result.push_back(HETP);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_13__H_p(double HETP, double N_ES) {
  std::vector<double> result;
  double H_p = (HETP * N_ES);
  result.push_back(H_p);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_13__N_ES(double HETP, double H_p) {
  std::vector<double> result;
  double N_ES = (H_p / HETP);
  result.push_back(N_ES);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_14__M(double P_0, double T, double W_E) {
  std::vector<double> result;
  double M =
      (((294.213699178261 * T) * std::pow(W_E, 2.0)) / std::pow(P_0, 2.0));
  result.push_back(M);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_14__P_0(double M, double T, double W_E) {
  std::vector<double> result;
  double P_0 = ((17.1526586620926 * W_E) / std::sqrt((M / T)));
  result.push_back(P_0);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_14__T(double M, double P_0, double W_E) {
  std::vector<double> result;
  double T = (((0.00339889 * M) * std::pow(P_0, 2.0)) / std::pow(W_E, 2.0));
  result.push_back(T);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_14__W_E(double M, double P_0, double T) {
  std::vector<double> result;
  double W_E = ((0.0583 * P_0) * std::sqrt((M / T)));
  result.push_back(W_E);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_15__M_1(double M_2, double P_0_1,
                                              double P_0_2, double a_M_12) {
  std::vector<double> result;
  double M_1 = ((-M_2) / std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
  result.push_back(M_1);
  M_1 = (M_2 / std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
  result.push_back(M_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_15__M_2(double M_1, double P_0_1,
                                              double P_0_2, double a_M_12) {
  std::vector<double> result;
  double M_2 = ((-M_1) * std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
  result.push_back(M_2);
  M_2 = (M_1 * std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
  result.push_back(M_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_15__P_0_1(double M_1, double M_2,
                                                double P_0_2, double a_M_12) {
  std::vector<double> result;
  double P_0_1 = ((P_0_2 * a_M_12) / std::pow((M_2 / M_1), (2.0 / 5.0)));
  result.push_back(P_0_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_15__P_0_2(double M_1, double M_2,
                                                double P_0_1, double a_M_12) {
  std::vector<double> result;
  double P_0_2 = ((P_0_1 * std::pow((M_2 / M_1), (2.0 / 5.0))) / a_M_12);
  result.push_back(P_0_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_15__a_M_12(double M_1, double M_2,
                                                 double P_0_1, double P_0_2) {
  std::vector<double> result;
  double a_M_12 = ((P_0_1 * std::pow((M_2 / M_1), (2.0 / 5.0))) / P_0_2);
  result.push_back(a_M_12);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_16__H_i(double p_i, double x_i) {
  std::vector<double> result;
  double H_i = (p_i / x_i);
  result.push_back(H_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_16__p_i(double H_i, double x_i) {
  std::vector<double> result;
  double p_i = (H_i * x_i);
  result.push_back(p_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_16__x_i(double H_i, double p_i) {
  std::vector<double> result;
  double x_i = (p_i / H_i);
  result.push_back(x_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_17__H_2_1(double H_2_3, double H_2_mi,
                                                double x_1, double x_3) {
  std::vector<double> result;
  double H_2_1 =
      std::exp(((((-x_3) * std::log(H_2_3)) + std::log(H_2_mi)) / x_1));
  result.push_back(H_2_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_17__H_2_3(double H_2_1, double H_2_mi,
                                                double x_1, double x_3) {
  std::vector<double> result;
  double H_2_3 =
      std::exp(((((-x_1) * std::log(H_2_1)) + std::log(H_2_mi)) / x_3));
  result.push_back(H_2_3);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_17__H_2_mi(double H_2_1, double H_2_3,
                                                 double x_1, double x_3) {
  std::vector<double> result;
  double H_2_mi = std::exp(((x_1 * std::log(H_2_1)) + (x_3 * std::log(H_2_3))));
  result.push_back(H_2_mi);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_17__x_1(double H_2_1, double H_2_3,
                                              double H_2_mi, double x_3) {
  std::vector<double> result;
  double x_1 =
      ((((-x_3) * std::log(H_2_3)) + std::log(H_2_mi)) / std::log(H_2_1));
  result.push_back(x_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_17__x_3(double H_2_1, double H_2_3,
                                              double H_2_mi, double x_1) {
  std::vector<double> result;
  double x_3 =
      ((((-x_1) * std::log(H_2_1)) + std::log(H_2_mi)) / std::log(H_2_3));
  result.push_back(x_3);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2a__K_1(double K_2, double alpha_1_2) {
  std::vector<double> result;
  double K_1 = (K_2 * alpha_1_2);
  result.push_back(K_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2a__K_2(double K_1, double alpha_1_2) {
  std::vector<double> result;
  double K_2 = (K_1 / alpha_1_2);
  result.push_back(K_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2a__alpha_1_2(double K_1, double K_2) {
  std::vector<double> result;
  double alpha_1_2 = (K_1 / K_2);
  result.push_back(alpha_1_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2b__K_1(double K_2, double x_1,
                                              double x_2, double y_1,
                                              double y_2) {
  std::vector<double> result;
  double K_1 = (((K_2 * x_2) * y_1) / (x_1 * y_2));
  result.push_back(K_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2b__K_2(double K_1, double x_1,
                                              double x_2, double y_1,
                                              double y_2) {
  std::vector<double> result;
  double K_2 = (((K_1 * x_1) * y_2) / (x_2 * y_1));
  result.push_back(K_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2b__x_1(double K_1, double K_2,
                                              double x_2, double y_1,
                                              double y_2) {
  std::vector<double> result;
  double x_1 = (((K_2 * x_2) * y_1) / (K_1 * y_2));
  result.push_back(x_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2b__x_2(double K_1, double K_2,
                                              double x_1, double y_1,
                                              double y_2) {
  std::vector<double> result;
  double x_2 = (((K_1 * x_1) * y_2) / (K_2 * y_1));
  result.push_back(x_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2b__y_1(double K_1, double K_2,
                                              double x_1, double x_2,
                                              double y_2) {
  std::vector<double> result;
  double y_1 = (((K_1 * x_1) * y_2) / (K_2 * x_2));
  result.push_back(y_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_2b__y_2(double K_1, double K_2,
                                              double x_1, double x_2,
                                              double y_1) {
  std::vector<double> result;
  double y_2 = (((K_2 * x_2) * y_1) / (K_1 * x_1));
  result.push_back(y_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_3__P_0_i(double p_i, double x_i) {
  std::vector<double> result;
  double P_0_i = (p_i / x_i);
  result.push_back(P_0_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_3__p_i(double P_0_i, double x_i) {
  std::vector<double> result;
  double p_i = (P_0_i * x_i);
  result.push_back(p_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_3__x_i(double P_0_i, double p_i) {
  std::vector<double> result;
  double x_i = (p_i / P_0_i);
  result.push_back(x_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_4__P(double P_0_i, double x_i,
                                           double y_i) {
  std::vector<double> result;
  double P = ((P_0_i * x_i) / y_i);
  result.push_back(P);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_4__P_0_i(double P, double x_i,
                                               double y_i) {
  std::vector<double> result;
  double P_0_i = ((P * y_i) / x_i);
  result.push_back(P_0_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_4__x_i(double P, double P_0_i,
                                             double y_i) {
  std::vector<double> result;
  double x_i = ((P * y_i) / P_0_i);
  result.push_back(x_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_4__y_i(double P, double P_0_i,
                                             double x_i) {
  std::vector<double> result;
  double y_i = ((P_0_i * x_i) / P);
  result.push_back(y_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_5__P_0_1(double P_0_2, double alpha_12) {
  std::vector<double> result;
  double P_0_1 = (P_0_2 * alpha_12);
  result.push_back(P_0_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_5__P_0_2(double P_0_1, double alpha_12) {
  std::vector<double> result;
  double P_0_2 = (P_0_1 / alpha_12);
  result.push_back(P_0_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_5__alpha_12(double P_0_1, double P_0_2) {
  std::vector<double> result;
  double alpha_12 = (P_0_1 / P_0_2);
  result.push_back(alpha_12);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_6__P_0_i(double gamma_i, double p_i,
                                               double x_i) {
  std::vector<double> result;
  double P_0_i = (p_i / (gamma_i * x_i));
  result.push_back(P_0_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_6__gamma_i(double P_0_i, double p_i,
                                                 double x_i) {
  std::vector<double> result;
  double gamma_i = (p_i / (P_0_i * x_i));
  result.push_back(gamma_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_6__p_i(double P_0_i, double gamma_i,
                                             double x_i) {
  std::vector<double> result;
  double p_i = ((P_0_i * gamma_i) * x_i);
  result.push_back(p_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_6__x_i(double P_0_i, double gamma_i,
                                             double p_i) {
  std::vector<double> result;
  double x_i = (p_i / (P_0_i * gamma_i));
  result.push_back(x_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_7__P(double P_0_i, double gamma_i,
                                           double x_i, double y_i) {
  std::vector<double> result;
  double P = (((P_0_i * gamma_i) * x_i) / y_i);
  result.push_back(P);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_7__P_0_i(double P, double gamma_i,
                                               double x_i, double y_i) {
  std::vector<double> result;
  double P_0_i = ((P * y_i) / (gamma_i * x_i));
  result.push_back(P_0_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_7__gamma_i(double P, double P_0_i,
                                                 double x_i, double y_i) {
  std::vector<double> result;
  double gamma_i = ((P * y_i) / (P_0_i * x_i));
  result.push_back(gamma_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_7__x_i(double P, double P_0_i,
                                             double gamma_i, double y_i) {
  std::vector<double> result;
  double x_i = ((P * y_i) / (P_0_i * gamma_i));
  result.push_back(x_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_7__y_i(double P, double P_0_i,
                                             double gamma_i, double x_i) {
  std::vector<double> result;
  double y_i = (((P_0_i * gamma_i) * x_i) / P);
  result.push_back(y_i);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_8__P_0_1(double P_0_2, double alpha_12,
                                               double gamma_1, double gamma_2) {
  std::vector<double> result;
  double P_0_1 = (((P_0_2 * alpha_12) * gamma_2) / gamma_1);
  result.push_back(P_0_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_8__P_0_2(double P_0_1, double alpha_12,
                                               double gamma_1, double gamma_2) {
  std::vector<double> result;
  double P_0_2 = ((P_0_1 * gamma_1) / (alpha_12 * gamma_2));
  result.push_back(P_0_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_8__alpha_12(double P_0_1, double P_0_2,
                                                  double gamma_1,
                                                  double gamma_2) {
  std::vector<double> result;
  double alpha_12 = ((P_0_1 * gamma_1) / (P_0_2 * gamma_2));
  result.push_back(alpha_12);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_8__gamma_1(double P_0_1, double P_0_2,
                                                 double alpha_12,
                                                 double gamma_2) {
  std::vector<double> result;
  double gamma_1 = (((P_0_2 * alpha_12) * gamma_2) / P_0_1);
  result.push_back(gamma_1);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_8__gamma_2(double P_0_1, double P_0_2,
                                                 double alpha_12,
                                                 double gamma_1) {
  std::vector<double> result;
  double gamma_2 = ((P_0_1 * gamma_1) / (P_0_2 * alpha_12));
  result.push_back(gamma_2);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_9__D(double L_0, double V_1) {
  std::vector<double> result;
  double D = ((-L_0) + V_1);
  result.push_back(D);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_9__L_0(double D, double V_1) {
  std::vector<double> result;
  double L_0 = 0.0;
  result.push_back(L_0);
  L_0 = ((-D) + V_1);
  result.push_back(L_0);
  return result;
}

std::vector<double> ProcessApp1_eqn_5_9__V_1(double D, double L_0) {
  std::vector<double> result;
  double V_1 = (D + L_0);
  result.push_back(V_1);
  return result;
}

#endif // VAKYUME_PROCESSAPP1_HPP
