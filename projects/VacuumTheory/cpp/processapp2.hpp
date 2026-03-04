#ifndef VAKYUME_PROCESSAPP2_HPP
#define VAKYUME_PROCESSAPP2_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> ProcessApp2_eqn_6_1__T_1(double T_2, double T_R, double c_p,
                                             double del_h_v, double w_1,
                                             double w_2, double w_v) {
  std::vector<double> result;
  double T_1 = (((((T_R * c_p) * w_1) + ((c_p * w_2) * ((-T_2) + T_R))) +
                 (del_h_v * w_v)) /
                (c_p * w_1));
  result.push_back(T_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_1__T_2(double T_1, double T_R, double c_p,
                                             double del_h_v, double w_1,
                                             double w_2, double w_v) {
  std::vector<double> result;
  double T_2 = (((((T_R * c_p) * w_2) + ((c_p * w_1) * ((-T_1) + T_R))) +
                 (del_h_v * w_v)) /
                (c_p * w_2));
  result.push_back(T_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_1__T_R(double T_1, double T_2, double c_p,
                                             double del_h_v, double w_1,
                                             double w_2, double w_v) {
  std::vector<double> result;
  double T_R =
      (((((T_1 * c_p) * w_1) + ((T_2 * c_p) * w_2)) - (del_h_v * w_v)) /
       (c_p * (w_1 + w_2)));
  result.push_back(T_R);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_1__c_p(double T_1, double T_2, double T_R,
                                             double del_h_v, double w_1,
                                             double w_2, double w_v) {
  std::vector<double> result;
  double c_p = ((del_h_v * w_v) /
                ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2)));
  result.push_back(c_p);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_1__del_h_v(double T_1, double T_2,
                                                 double T_R, double c_p,
                                                 double w_1, double w_2,
                                                 double w_v) {
  std::vector<double> result;
  double del_h_v =
      ((c_p * ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2))) /
       w_v);
  result.push_back(del_h_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_1__w_1(double T_1, double T_2, double T_R,
                                             double c_p, double del_h_v,
                                             double w_2, double w_v) {
  std::vector<double> result;
  double w_1 =
      ((((((-T_2) * c_p) * w_2) + ((T_R * c_p) * w_2)) + (del_h_v * w_v)) /
       (c_p * (T_1 - T_R)));
  result.push_back(w_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_1__w_2(double T_1, double T_2, double T_R,
                                             double c_p, double del_h_v,
                                             double w_1, double w_v) {
  std::vector<double> result;
  double w_2 =
      ((((((-T_1) * c_p) * w_1) + ((T_R * c_p) * w_1)) + (del_h_v * w_v)) /
       (c_p * (T_2 - T_R)));
  result.push_back(w_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_1__w_v(double T_1, double T_2, double T_R,
                                             double c_p, double del_h_v,
                                             double w_1, double w_2) {
  std::vector<double> result;
  double w_v =
      ((c_p * ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2))) /
       del_h_v);
  result.push_back(w_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_10__A(double dV_dt, double delta_P,
                                            double mu, double r_c, double s,
                                            double tau) {
  std::vector<double> result;
  double A = ((((dV_dt * std::pow(delta_P, (s - 1.0))) * mu) * r_c) * tau);
  result.push_back(A);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_10__dV_dt(double A, double delta_P,
                                                double mu, double r_c, double s,
                                                double tau) {
  std::vector<double> result;
  double dV_dt = ((A * std::pow(delta_P, (1.0 - s))) / ((mu * r_c) * tau));
  result.push_back(dV_dt);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_10__delta_P(double A, double dV_dt,
                                                  double mu, double r_c,
                                                  double s, double tau) {
  std::vector<double> result;
  double delta_P =
      std::pow(((((dV_dt * mu) * r_c) * tau) / A), ((-1.0) / (s - 1.0)));
  result.push_back(delta_P);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_10__mu(double A, double dV_dt,
                                             double delta_P, double r_c,
                                             double s, double tau) {
  std::vector<double> result;
  double mu = ((A * std::pow(delta_P, (1.0 - s))) / ((dV_dt * r_c) * tau));
  result.push_back(mu);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_10__r_c(double A, double dV_dt,
                                              double delta_P, double mu,
                                              double s, double tau) {
  std::vector<double> result;
  double r_c = ((A * std::pow(delta_P, (1.0 - s))) / ((dV_dt * mu) * tau));
  result.push_back(r_c);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_10__s(double A, double dV_dt,
                                            double delta_P, double mu,
                                            double r_c, double tau) {
  std::vector<double> result;
  double s = (std::log(((A * delta_P) / (((dV_dt * mu) * r_c) * tau))) /
              std::log(delta_P));
  result.push_back(s);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_10__tau(double A, double dV_dt,
                                              double delta_P, double mu,
                                              double r_c, double s) {
  std::vector<double> result;
  double tau = ((A * std::pow(delta_P, (1.0 - s))) / ((dV_dt * mu) * r_c));
  result.push_back(tau);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_11a__A_d(double delta_T, double delta_h_i,
                                               double delta_m, double h_d,
                                               double m_b, double t_R) {
  std::vector<double> result;
  double A_d = (((delta_h_i * delta_m) * m_b) / ((delta_T * h_d) * t_R));
  result.push_back(A_d);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_11a__delta_T(double A_d, double delta_h_i,
                                                   double delta_m, double h_d,
                                                   double m_b, double t_R) {
  std::vector<double> result;
  double delta_T = (((delta_h_i * delta_m) * m_b) / ((A_d * h_d) * t_R));
  result.push_back(delta_T);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_11a__delta_h_i(double A_d, double delta_T,
                                                     double delta_m, double h_d,
                                                     double m_b, double t_R) {
  std::vector<double> result;
  double delta_h_i = ((((A_d * delta_T) * h_d) * t_R) / (delta_m * m_b));
  result.push_back(delta_h_i);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_11a__delta_m(double A_d, double delta_T,
                                                   double delta_h_i, double h_d,
                                                   double m_b, double t_R) {
  std::vector<double> result;
  double delta_m = ((((A_d * delta_T) * h_d) * t_R) / (delta_h_i * m_b));
  result.push_back(delta_m);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_11a__h_d(double A_d, double delta_T,
                                               double delta_h_i, double delta_m,
                                               double m_b, double t_R) {
  std::vector<double> result;
  double h_d = (((delta_h_i * delta_m) * m_b) / ((A_d * delta_T) * t_R));
  result.push_back(h_d);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_11a__m_b(double A_d, double delta_T,
                                               double delta_h_i, double delta_m,
                                               double h_d, double t_R) {
  std::vector<double> result;
  double m_b = ((((A_d * delta_T) * h_d) * t_R) / (delta_h_i * delta_m));
  result.push_back(m_b);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_11a__t_R(double A_d, double delta_T,
                                               double delta_h_i, double delta_m,
                                               double h_d, double m_b) {
  std::vector<double> result;
  double t_R = (((delta_h_i * delta_m) * m_b) / ((A_d * delta_T) * h_d));
  result.push_back(t_R);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_2__Q_v(double T_1, double T_2, double T_R,
                                             double c_p, double w_1,
                                             double w_2) {
  std::vector<double> result;
  double Q_v =
      ((c_p * ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2))) /
       12000.0);
  result.push_back(Q_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_2__T_1(double Q_v, double T_2, double T_R,
                                             double c_p, double w_1,
                                             double w_2) {
  std::vector<double> result;
  double T_1 = ((((12000.0 * Q_v) + ((T_R * c_p) * w_1)) +
                 ((c_p * w_2) * ((-T_2) + T_R))) /
                (c_p * w_1));
  result.push_back(T_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_2__T_2(double Q_v, double T_1, double T_R,
                                             double c_p, double w_1,
                                             double w_2) {
  std::vector<double> result;
  double T_2 = ((((12000.0 * Q_v) + ((T_R * c_p) * w_2)) +
                 ((c_p * w_1) * ((-T_1) + T_R))) /
                (c_p * w_2));
  result.push_back(T_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_2__T_R(double Q_v, double T_1, double T_2,
                                             double c_p, double w_1,
                                             double w_2) {
  std::vector<double> result;
  double T_R =
      (((((-12000.0) * Q_v) + ((T_1 * c_p) * w_1)) + ((T_2 * c_p) * w_2)) /
       (c_p * (w_1 + w_2)));
  result.push_back(T_R);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_2__c_p(double Q_v, double T_1, double T_2,
                                             double T_R, double w_1,
                                             double w_2) {
  std::vector<double> result;
  double c_p = ((12000.0 * Q_v) /
                ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2)));
  result.push_back(c_p);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_2__w_1(double Q_v, double T_1, double T_2,
                                             double T_R, double c_p,
                                             double w_2) {
  std::vector<double> result;
  double w_1 =
      ((((12000.0 * Q_v) - ((T_2 * c_p) * w_2)) + ((T_R * c_p) * w_2)) /
       (c_p * (T_1 - T_R)));
  result.push_back(w_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_2__w_2(double Q_v, double T_1, double T_2,
                                             double T_R, double c_p,
                                             double w_1) {
  std::vector<double> result;
  double w_2 =
      ((((12000.0 * Q_v) - ((T_1 * c_p) * w_1)) + ((T_R * c_p) * w_1)) /
       (c_p * (T_2 - T_R)));
  result.push_back(w_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_4__Q_v(double delta_h_v, double w_v) {
  std::vector<double> result;
  double Q_v = ((delta_h_v * w_v) / 12000.0);
  result.push_back(Q_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_4__delta_h_v(double Q_v, double w_v) {
  std::vector<double> result;
  double delta_h_v = ((12000.0 * Q_v) / w_v);
  result.push_back(delta_h_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_4__w_v(double Q_v, double delta_h_v) {
  std::vector<double> result;
  double w_v = ((12000.0 * Q_v) / delta_h_v);
  result.push_back(w_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__C_1(double C_2, double T_1, double T_2,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double m_b,
                                             double m_v) {
  std::vector<double> result;
  double C_1 = (((((C_2 * delta_h_c) * m_b) + ((c_p * m_b) * ((-T_1) + T_2))) +
                 (delta_h_v * m_v)) /
                (delta_h_c * m_b));
  result.push_back(C_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__C_2(double C_1, double T_1, double T_2,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double m_b,
                                             double m_v) {
  std::vector<double> result;
  double C_2 = (((((C_1 * delta_h_c) * m_b) + ((c_p * m_b) * (T_1 - T_2))) -
                 (delta_h_v * m_v)) /
                (delta_h_c * m_b));
  result.push_back(C_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__T_1(double C_1, double C_2, double T_2,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double m_b,
                                             double m_v) {
  std::vector<double> result;
  double T_1 = (((((T_2 * c_p) * m_b) + ((delta_h_c * m_b) * ((-C_1) + C_2))) +
                 (delta_h_v * m_v)) /
                (c_p * m_b));
  result.push_back(T_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__T_2(double C_1, double C_2, double T_1,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double m_b,
                                             double m_v) {
  std::vector<double> result;
  double T_2 = (((((T_1 * c_p) * m_b) + ((delta_h_c * m_b) * (C_1 - C_2))) -
                 (delta_h_v * m_v)) /
                (c_p * m_b));
  result.push_back(T_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__c_p(double C_1, double C_2, double T_1,
                                             double T_2, double delta_h_c,
                                             double delta_h_v, double m_b,
                                             double m_v) {
  std::vector<double> result;
  double c_p = ((((((-C_1) * delta_h_c) * m_b) + ((C_2 * delta_h_c) * m_b)) +
                 (delta_h_v * m_v)) /
                (m_b * (T_1 - T_2)));
  result.push_back(c_p);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__delta_h_c(double C_1, double C_2,
                                                   double T_1, double T_2,
                                                   double c_p, double delta_h_v,
                                                   double m_b, double m_v) {
  std::vector<double> result;
  double delta_h_c =
      ((((((-T_1) * c_p) * m_b) + ((T_2 * c_p) * m_b)) + (delta_h_v * m_v)) /
       (m_b * (C_1 - C_2)));
  result.push_back(delta_h_c);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__delta_h_v(double C_1, double C_2,
                                                   double T_1, double T_2,
                                                   double c_p, double delta_h_c,
                                                   double m_b, double m_v) {
  std::vector<double> result;
  double delta_h_v =
      ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) -
               (T_2 * c_p))) /
       m_v);
  result.push_back(delta_h_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__m_b(double C_1, double C_2, double T_1,
                                             double T_2, double c_p,
                                             double delta_h_c, double delta_h_v,
                                             double m_v) {
  std::vector<double> result;
  double m_b =
      ((delta_h_v * m_v) /
       ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p)));
  result.push_back(m_b);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_7__m_v(double C_1, double C_2, double T_1,
                                             double T_2, double c_p,
                                             double delta_h_c, double delta_h_v,
                                             double m_b) {
  std::vector<double> result;
  double m_v =
      ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) -
               (T_2 * c_p))) /
       delta_h_v);
  result.push_back(m_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__C_1(double C_2, double T_1, double T_2,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double delta_t,
                                             double m_b, double w_v) {
  std::vector<double> result;
  double C_1 = (((((C_2 * delta_h_c) * m_b) + ((c_p * m_b) * ((-T_1) + T_2))) +
                 ((delta_h_v * delta_t) * w_v)) /
                (delta_h_c * m_b));
  result.push_back(C_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__C_2(double C_1, double T_1, double T_2,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double delta_t,
                                             double m_b, double w_v) {
  std::vector<double> result;
  double C_2 = (((((C_1 * delta_h_c) * m_b) + ((c_p * m_b) * (T_1 - T_2))) -
                 ((delta_h_v * delta_t) * w_v)) /
                (delta_h_c * m_b));
  result.push_back(C_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__T_1(double C_1, double C_2, double T_2,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double delta_t,
                                             double m_b, double w_v) {
  std::vector<double> result;
  double T_1 = (((((T_2 * c_p) * m_b) + ((delta_h_c * m_b) * ((-C_1) + C_2))) +
                 ((delta_h_v * delta_t) * w_v)) /
                (c_p * m_b));
  result.push_back(T_1);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__T_2(double C_1, double C_2, double T_1,
                                             double c_p, double delta_h_c,
                                             double delta_h_v, double delta_t,
                                             double m_b, double w_v) {
  std::vector<double> result;
  double T_2 = (((((T_1 * c_p) * m_b) + ((delta_h_c * m_b) * (C_1 - C_2))) -
                 ((delta_h_v * delta_t) * w_v)) /
                (c_p * m_b));
  result.push_back(T_2);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__c_p(double C_1, double C_2, double T_1,
                                             double T_2, double delta_h_c,
                                             double delta_h_v, double delta_t,
                                             double m_b, double w_v) {
  std::vector<double> result;
  double c_p = ((((((-C_1) * delta_h_c) * m_b) + ((C_2 * delta_h_c) * m_b)) +
                 ((delta_h_v * delta_t) * w_v)) /
                (m_b * (T_1 - T_2)));
  result.push_back(c_p);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__delta_h_c(double C_1, double C_2,
                                                   double T_1, double T_2,
                                                   double c_p, double delta_h_v,
                                                   double delta_t, double m_b,
                                                   double w_v) {
  std::vector<double> result;
  double delta_h_c = ((((((-T_1) * c_p) * m_b) + ((T_2 * c_p) * m_b)) +
                       ((delta_h_v * delta_t) * w_v)) /
                      (m_b * (C_1 - C_2)));
  result.push_back(delta_h_c);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__delta_h_v(double C_1, double C_2,
                                                   double T_1, double T_2,
                                                   double c_p, double delta_h_c,
                                                   double delta_t, double m_b,
                                                   double w_v) {
  std::vector<double> result;
  double delta_h_v =
      ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) -
               (T_2 * c_p))) /
       (delta_t * w_v));
  result.push_back(delta_h_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__delta_t(double C_1, double C_2,
                                                 double T_1, double T_2,
                                                 double c_p, double delta_h_c,
                                                 double delta_h_v, double m_b,
                                                 double w_v) {
  std::vector<double> result;
  double delta_t =
      ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) -
               (T_2 * c_p))) /
       (delta_h_v * w_v));
  result.push_back(delta_t);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__m_b(double C_1, double C_2, double T_1,
                                             double T_2, double c_p,
                                             double delta_h_c, double delta_h_v,
                                             double delta_t, double w_v) {
  std::vector<double> result;
  double m_b =
      (((delta_h_v * delta_t) * w_v) /
       ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p)));
  result.push_back(m_b);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_8__w_v(double C_1, double C_2, double T_1,
                                             double T_2, double c_p,
                                             double delta_h_c, double delta_h_v,
                                             double delta_t, double m_b) {
  std::vector<double> result;
  double w_v =
      ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) -
               (T_2 * c_p))) /
       (delta_h_v * delta_t));
  result.push_back(w_v);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_9__A(double dV_dt, double delta_P,
                                           double m, double mu, double r,
                                           double r_M) {
  std::vector<double> result;
  double A =
      (((dV_dt * r_M) -
        std::sqrt(
            (dV_dt * ((dV_dt * std::pow(r_M, 2.0)) +
                      ((((4.0 * std::pow(delta_P, 2.0)) * m) * mu) * r))))) /
       (2.0 * delta_P));
  result.push_back(A);
  A = (((dV_dt * r_M) +
        std::sqrt(
            (dV_dt * ((dV_dt * std::pow(r_M, 2.0)) +
                      ((((4.0 * std::pow(delta_P, 2.0)) * m) * mu) * r))))) /
       (2.0 * delta_P));
  result.push_back(A);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_9__dV_dt(double A, double delta_P,
                                               double m, double mu, double r,
                                               double r_M) {
  std::vector<double> result;
  double dV_dt =
      ((std::pow(A, 2.0) * delta_P) / ((A * r_M) + (((delta_P * m) * mu) * r)));
  result.push_back(dV_dt);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_9__delta_P(double A, double dV_dt,
                                                 double m, double mu, double r,
                                                 double r_M) {
  std::vector<double> result;
  double delta_P =
      (((A * dV_dt) * r_M) / (std::pow(A, 2.0) - (((dV_dt * m) * mu) * r)));
  result.push_back(delta_P);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_9__m(double A, double dV_dt,
                                           double delta_P, double mu, double r,
                                           double r_M) {
  std::vector<double> result;
  double m =
      ((A * ((A * delta_P) - (dV_dt * r_M))) / (((dV_dt * delta_P) * mu) * r));
  result.push_back(m);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_9__mu(double A, double dV_dt,
                                            double delta_P, double m, double r,
                                            double r_M) {
  std::vector<double> result;
  double mu =
      ((A * ((A * delta_P) - (dV_dt * r_M))) / (((dV_dt * delta_P) * m) * r));
  result.push_back(mu);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_9__r(double A, double dV_dt,
                                           double delta_P, double m, double mu,
                                           double r_M) {
  std::vector<double> result;
  double r =
      ((A * ((A * delta_P) - (dV_dt * r_M))) / (((dV_dt * delta_P) * m) * mu));
  result.push_back(r);
  return result;
}

std::vector<double> ProcessApp2_eqn_6_9__r_M(double A, double dV_dt,
                                             double delta_P, double m,
                                             double mu, double r) {
  std::vector<double> result;
  double r_M = (((A * delta_P) / dV_dt) - ((((delta_P * m) * mu) * r) / A));
  result.push_back(r_M);
  return result;
}

#endif // VAKYUME_PROCESSAPP2_HPP
