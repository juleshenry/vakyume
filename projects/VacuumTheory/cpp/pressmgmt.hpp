#ifndef VAKYUME_PRESSMGMT_HPP
#define VAKYUME_PRESSMGMT_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> PressMgmt_eqn_3_1__Abs_Pressure(double BarometricPressure,
                                                    double Vacuum) {
  std::vector<double> result;
  double Abs_Pressure = (BarometricPressure - Vacuum);
  result.push_back(Abs_Pressure);
  return result;
}

std::vector<double> PressMgmt_eqn_3_1__BarometricPressure(double Abs_Pressure,
                                                          double Vacuum) {
  std::vector<double> result;
  double BarometricPressure = (Abs_Pressure + Vacuum);
  result.push_back(BarometricPressure);
  return result;
}

std::vector<double> PressMgmt_eqn_3_1__Vacuum(double Abs_Pressure,
                                              double BarometricPressure) {
  std::vector<double> result;
  double Vacuum = ((-Abs_Pressure) + BarometricPressure);
  result.push_back(Vacuum);
  return result;
}

std::vector<double> PressMgmt_eqn_3_11__A_C(double H_2, double P, double V) {
  std::vector<double> result;
  double A_C = ((P * V) / std::pow(H_2, 2.0));
  result.push_back(A_C);
  return result;
}

std::vector<double> PressMgmt_eqn_3_11__H_2(double A_C, double P, double V) {
  std::vector<double> result;
  double H_2 = (-std::sqrt(((P * V) / A_C)));
  result.push_back(H_2);
  H_2 = std::sqrt(((P * V) / A_C));
  result.push_back(H_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_11__P(double A_C, double H_2, double V) {
  std::vector<double> result;
  double P = ((A_C * std::pow(H_2, 2.0)) / V);
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_11__V(double A_C, double H_2, double P) {
  std::vector<double> result;
  double V = ((A_C * std::pow(H_2, 2.0)) / P);
  result.push_back(V);
  return result;
}

std::vector<double> PressMgmt_eqn_3_12__H_2(double KAPPA_1, double P) {
  std::vector<double> result;
  double H_2 = (-std::sqrt((P / KAPPA_1)));
  result.push_back(H_2);
  H_2 = std::sqrt((P / KAPPA_1));
  result.push_back(H_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_12__KAPPA_1(double H_2, double P) {
  std::vector<double> result;
  double KAPPA_1 = (P / std::pow(H_2, 2.0));
  result.push_back(KAPPA_1);
  return result;
}

std::vector<double> PressMgmt_eqn_3_12__P(double H_2, double KAPPA_1) {
  std::vector<double> result;
  double P = (std::pow(H_2, 2.0) * KAPPA_1);
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_13__H_1(double H_2, double KAPPA_2,
                                            double P) {
  std::vector<double> result;
  double H_1 = (H_2 - (P / KAPPA_2));
  result.push_back(H_1);
  return result;
}

std::vector<double> PressMgmt_eqn_3_13__H_2(double H_1, double KAPPA_2,
                                            double P) {
  std::vector<double> result;
  double H_2 = (H_1 + (P / KAPPA_2));
  result.push_back(H_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_13__KAPPA_2(double H_1, double H_2,
                                                double P) {
  std::vector<double> result;
  double KAPPA_2 = ((-P) / (H_1 - H_2));
  result.push_back(KAPPA_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_13__P(double H_1, double H_2,
                                          double KAPPA_2) {
  std::vector<double> result;
  double P = (KAPPA_2 * ((-H_1) + H_2));
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_15__V_PMIN() {
  std::vector<double> result;
  double V_PMIN = 0.785398163397448;
  result.push_back(V_PMIN);
  return result;
}

std::vector<double> PressMgmt_eqn_3_16__V_div_V_P_MAX() {
  std::vector<double> result;
  double V_div_V_P_MAX = 254647.908947033;
  result.push_back(V_div_V_P_MAX);
  return result;
}

std::vector<double> PressMgmt_eqn_3_17__P_MIN() {
  std::vector<double> result;
  double P_MIN = 3.92699081698724e-06;
  result.push_back(P_MIN);
  return result;
}

std::vector<double> PressMgmt_eqn_3_2__G(double G_C, double H, double P,
                                         double rho) {
  std::vector<double> result;
  double G = (((G_C * H) * P) * rho);
  result.push_back(G);
  return result;
}

std::vector<double> PressMgmt_eqn_3_2__G_C(double G, double H, double P,
                                           double rho) {
  std::vector<double> result;
  double G_C = (G / ((H * P) * rho));
  result.push_back(G_C);
  return result;
}

std::vector<double> PressMgmt_eqn_3_2__H(double G, double G_C, double P,
                                         double rho) {
  std::vector<double> result;
  double H = (G / ((G_C * P) * rho));
  result.push_back(H);
  return result;
}

std::vector<double> PressMgmt_eqn_3_2__P(double G, double G_C, double H,
                                         double rho) {
  std::vector<double> result;
  double P = (G / ((G_C * H) * rho));
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_2__rho(double G, double G_C, double H,
                                           double P) {
  std::vector<double> result;
  double rho = (G / ((G_C * H) * P));
  result.push_back(rho);
  return result;
}

std::vector<double> PressMgmt_eqn_3_3__H_1(double H_2, double P, double P_P) {
  std::vector<double> result;
  double H_1 = ((H_2 + P) - P_P);
  result.push_back(H_1);
  return result;
}

std::vector<double> PressMgmt_eqn_3_3__H_2(double H_1, double P, double P_P) {
  std::vector<double> result;
  double H_2 = ((H_1 - P) + P_P);
  result.push_back(H_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_3__P(double H_1, double H_2, double P_P) {
  std::vector<double> result;
  double P = ((H_1 - H_2) + P_P);
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_3__P_P(double H_1, double H_2, double P) {
  std::vector<double> result;
  double P_P = (((-H_1) + H_2) + P);
  result.push_back(P_P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_4__KAPPA(double P, double V) {
  std::vector<double> result;
  double KAPPA = (P * V);
  result.push_back(KAPPA);
  return result;
}

std::vector<double> PressMgmt_eqn_3_4__P(double KAPPA, double V) {
  std::vector<double> result;
  double P = (KAPPA / V);
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_4__V(double KAPPA, double P) {
  std::vector<double> result;
  double V = (KAPPA / P);
  result.push_back(V);
  return result;
}

std::vector<double> PressMgmt_eqn_3_5__P(double P_P, double V, double V_P) {
  std::vector<double> result;
  double P = ((P_P * V_P) / V);
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_5__P_P(double P, double V, double V_P) {
  std::vector<double> result;
  double P_P = ((P * V) / V_P);
  result.push_back(P_P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_5__V(double P, double P_P, double V_P) {
  std::vector<double> result;
  double V = ((P_P * V_P) / P);
  result.push_back(V);
  return result;
}

std::vector<double> PressMgmt_eqn_3_5__V_P(double P, double P_P, double V) {
  std::vector<double> result;
  double V_P = ((P * V) / P_P);
  result.push_back(V_P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_6__H_1(double H_2, double P, double V,
                                           double V_P) {
  std::vector<double> result;
  double H_1 = ((H_2 - ((P * V) / V_P)) + P);
  result.push_back(H_1);
  return result;
}

std::vector<double> PressMgmt_eqn_3_6__H_2(double H_1, double P, double V,
                                           double V_P) {
  std::vector<double> result;
  double H_2 = ((H_1 + ((P * V) / V_P)) - P);
  result.push_back(H_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_6__P(double H_1, double H_2, double V,
                                         double V_P) {
  std::vector<double> result;
  double P = ((V_P * ((-H_1) + H_2)) / (V - V_P));
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_6__V(double H_1, double H_2, double P,
                                         double V_P) {
  std::vector<double> result;
  double V = ((V_P * (((-H_1) + H_2) + P)) / P);
  result.push_back(V);
  return result;
}

std::vector<double> PressMgmt_eqn_3_6__V_P(double H_1, double H_2, double P,
                                           double V) {
  std::vector<double> result;
  double V_P = ((P * V) / (((-H_1) + H_2) + P));
  result.push_back(V_P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_8__A_C(double H_2, double V_P) {
  std::vector<double> result;
  double A_C = (V_P / H_2);
  result.push_back(A_C);
  return result;
}

std::vector<double> PressMgmt_eqn_3_8__H_2(double A_C, double V_P) {
  std::vector<double> result;
  double H_2 = (V_P / A_C);
  result.push_back(H_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_8__V_P(double A_C, double H_2) {
  std::vector<double> result;
  double V_P = (A_C * H_2);
  result.push_back(V_P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_9__A_C(double H_1, double H_2, double P,
                                           double V) {
  std::vector<double> result;
  double A_C = ((P * V) / (H_2 * (((-H_1) + H_2) + P)));
  result.push_back(A_C);
  return result;
}

std::vector<double> PressMgmt_eqn_3_9__H_1(double A_C, double H_2, double P,
                                           double V) {
  std::vector<double> result;
  double H_1 = ((H_2 + P) - ((P * V) / (A_C * H_2)));
  result.push_back(H_1);
  return result;
}

std::vector<double> PressMgmt_eqn_3_9__H_2(double A_C, double H_1, double P,
                                           double V) {
  std::vector<double> result;
  double H_2 =
      (((A_C * (H_1 - P)) - std::sqrt((A_C * ((((A_C * std::pow(H_1, 2.0)) -
                                                (((2.0 * A_C) * H_1) * P)) +
                                               (A_C * std::pow(P, 2.0))) +
                                              ((4.0 * P) * V))))) /
       (2.0 * A_C));
  result.push_back(H_2);
  H_2 = (((A_C * (H_1 - P)) + std::sqrt((A_C * ((((A_C * std::pow(H_1, 2.0)) -
                                                  (((2.0 * A_C) * H_1) * P)) +
                                                 (A_C * std::pow(P, 2.0))) +
                                                ((4.0 * P) * V))))) /
         (2.0 * A_C));
  result.push_back(H_2);
  return result;
}

std::vector<double> PressMgmt_eqn_3_9__P(double A_C, double H_1, double H_2,
                                         double V) {
  std::vector<double> result;
  double P = (((A_C * H_2) * (H_1 - H_2)) / ((A_C * H_2) - V));
  result.push_back(P);
  return result;
}

std::vector<double> PressMgmt_eqn_3_9__V(double A_C, double H_1, double H_2,
                                         double P) {
  std::vector<double> result;
  double V = (((A_C * H_2) * (((-H_1) + H_2) + P)) / P);
  result.push_back(V);
  return result;
}

#endif // VAKYUME_PRESSMGMT_HPP
