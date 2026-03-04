#ifndef VAKYUME_FLUIDFLOWVACUUMLINES_HPP
#define VAKYUME_FLUIDFLOWVACUUMLINES_HPP

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std::complex_literals;

std::vector<double> FluidFlowVacuumLines_eqn_2_1__D(double Re, double mu,
                                                    double rho, double v) {
  std::vector<double> result;
  double D = ((Re * mu) / (rho * v));
  result.push_back(D);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_1__Re(double D, double mu,
                                                     double rho, double v) {
  std::vector<double> result;
  double Re = (((D * rho) * v) / mu);
  result.push_back(Re);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_1__mu(double D, double Re,
                                                     double rho, double v) {
  std::vector<double> result;
  double mu = (((D * rho) * v) / Re);
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_1__rho(double D, double Re,
                                                      double mu, double v) {
  std::vector<double> result;
  double rho = ((Re * mu) / (D * v));
  result.push_back(rho);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_1__v(double D, double Re,
                                                    double mu, double rho) {
  std::vector<double> result;
  double v = ((Re * mu) / (D * rho));
  result.push_back(v);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_10__Suc_Pres(double delta_P,
                                                            double oper_press) {
  std::vector<double> result;
  double Suc_Pres = ((-delta_P) + oper_press);
  result.push_back(Suc_Pres);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_10__delta_P(double Suc_Pres,
                                                           double oper_press) {
  std::vector<double> result;
  double delta_P = ((-Suc_Pres) + oper_press);
  result.push_back(delta_P);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_10__oper_press(double Suc_Pres,
                                                              double delta_P) {
  std::vector<double> result;
  double oper_press = (Suc_Pres + delta_P);
  result.push_back(oper_press);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_11__D(double L, double f,
                                                     double g_c, double h_r,
                                                     double v) {
  std::vector<double> result;
  double D = (((L * f) * std::pow(v, 2.0)) / ((2.0 * g_c) * h_r));
  result.push_back(D);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_11__L(double D, double f,
                                                     double g_c, double h_r,
                                                     double v) {
  std::vector<double> result;
  double L = ((((2.0 * D) * g_c) * h_r) / (f * std::pow(v, 2.0)));
  result.push_back(L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_11__f(double D, double L,
                                                     double g_c, double h_r,
                                                     double v) {
  std::vector<double> result;
  double f = ((((2.0 * D) * g_c) * h_r) / (L * std::pow(v, 2.0)));
  result.push_back(f);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_11__g_c(double D, double L,
                                                       double f, double h_r,
                                                       double v) {
  std::vector<double> result;
  double g_c = (((L * f) * std::pow(v, 2.0)) / ((2.0 * D) * h_r));
  result.push_back(g_c);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_11__h_r(double D, double L,
                                                       double f, double g_c,
                                                       double v) {
  std::vector<double> result;
  double h_r = (((L * f) * std::pow(v, 2.0)) / ((2.0 * D) * g_c));
  result.push_back(h_r);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_11__v(double D, double L,
                                                     double f, double g_c,
                                                     double h_r) {
  std::vector<double> result;
  double v = ((-std::sqrt(2.0)) * std::sqrt((((D * g_c) * h_r) / (L * f))));
  result.push_back(v);
  v = (std::sqrt(2.0) * std::sqrt((((D * g_c) * h_r) / (L * f))));
  result.push_back(v);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_12__L(double d, double delta_P,
                                                     double f, double g,
                                                     double rho, double v) {
  std::vector<double> result;
  double L = ((((0.464037122969838 * d) * delta_P) * g) /
              ((f * rho) * std::pow(v, 2.0)));
  result.push_back(L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_12__d(double L, double delta_P,
                                                     double f, double g,
                                                     double rho, double v) {
  std::vector<double> result;
  double d = (((((2.155 * L) * f) * rho) * std::pow(v, 2.0)) / (delta_P * g));
  result.push_back(d);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_12__delta_P(double L, double d,
                                                           double f, double g,
                                                           double rho,
                                                           double v) {
  std::vector<double> result;
  double delta_P = (((((2.155 * L) * f) * rho) * std::pow(v, 2.0)) / (d * g));
  result.push_back(delta_P);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_12__f(double L, double d,
                                                     double delta_P, double g,
                                                     double rho, double v) {
  std::vector<double> result;
  double f = ((((0.464037122969838 * d) * delta_P) * g) /
              ((L * rho) * std::pow(v, 2.0)));
  result.push_back(f);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_12__g(double L, double d,
                                                     double delta_P, double f,
                                                     double rho, double v) {
  std::vector<double> result;
  double g = (((((2.155 * L) * f) * rho) * std::pow(v, 2.0)) / (d * delta_P));
  result.push_back(g);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_12__rho(double L, double d,
                                                       double delta_P, double f,
                                                       double g, double v) {
  std::vector<double> result;
  double rho = ((((0.464037122969838 * d) * delta_P) * g) /
                ((L * f) * std::pow(v, 2.0)));
  result.push_back(rho);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_12__v(double L, double d,
                                                     double delta_P, double f,
                                                     double g, double rho) {
  std::vector<double> result;
  double v = ((-0.681202703290172) *
              std::sqrt((((d * delta_P) * g) / ((L * f) * rho))));
  result.push_back(v);
  v = (0.681202703290172 * std::sqrt((((d * delta_P) * g) / ((L * f) * rho))));
  result.push_back(v);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_13__L(double d, double delta_P,
                                                     double f, double q,
                                                     double rho) {
  std::vector<double> result;
  double L = (((0.465116279069767 * std::pow(d, 5.0)) * delta_P) /
              ((f * std::pow(q, 2.0)) * rho));
  result.push_back(L);
  return result;
}

std::vector<std::complex<double>>
FluidFlowVacuumLines_eqn_2_13__d(double L, double delta_P, double f, double q,
                                 double rho) {
  std::vector<std::complex<double>> result;
  std::complex<double> d =
      (1.16543402167043 *
       std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0)));
  result.push_back(d);
  d = (((-0.942855929354115) *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))) -
       ((0.685024930457783 * std::complex<double>(0.0, 1.0)) *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))));
  result.push_back(d);
  d = (((-0.942855929354115) *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))) +
       ((0.685024930457783 * std::complex<double>(0.0, 1.0)) *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))));
  result.push_back(d);
  d = ((0.360138918518902 *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))) -
       ((1.10839362062173 * std::complex<double>(0.0, 1.0)) *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))));
  result.push_back(d);
  d = ((0.360138918518902 *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))) +
       ((1.10839362062173 * std::complex<double>(0.0, 1.0)) *
        std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P),
                 (1.0 / 5.0))));
  result.push_back(d);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_13__delta_P(double L, double d,
                                                           double f, double q,
                                                           double rho) {
  std::vector<double> result;
  double delta_P =
      (((((2.15 * L) * f) * std::pow(q, 2.0)) * rho) / std::pow(d, 5.0));
  result.push_back(delta_P);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_13__f(double L, double d,
                                                     double delta_P, double q,
                                                     double rho) {
  std::vector<double> result;
  double f = (((0.465116279069767 * std::pow(d, 5.0)) * delta_P) /
              ((L * std::pow(q, 2.0)) * rho));
  result.push_back(f);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_13__q(double L, double d,
                                                     double delta_P, double f,
                                                     double rho) {
  std::vector<double> result;
  double q = ((-0.681994339470473) *
              std::sqrt(((std::pow(d, 5.0) * delta_P) / ((L * f) * rho))));
  result.push_back(q);
  q = (0.681994339470473 *
       std::sqrt(((std::pow(d, 5.0) * delta_P) / ((L * f) * rho))));
  result.push_back(q);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_13__rho(double L, double d,
                                                       double delta_P, double f,
                                                       double q) {
  std::vector<double> result;
  double rho = (((0.465116279069767 * std::pow(d, 5.0)) * delta_P) /
                ((L * f) * std::pow(q, 2.0)));
  result.push_back(rho);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_14__M(double R, double T,
                                                     double g_c, double k,
                                                     double v_s) {
  std::vector<double> result;
  double M = ((((R * T) * g_c) * k) / std::pow(v_s, 2.0));
  result.push_back(M);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_14__R(double M, double T,
                                                     double g_c, double k,
                                                     double v_s) {
  std::vector<double> result;
  double R = ((M * std::pow(v_s, 2.0)) / ((T * g_c) * k));
  result.push_back(R);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_14__T(double M, double R,
                                                     double g_c, double k,
                                                     double v_s) {
  std::vector<double> result;
  double T = ((M * std::pow(v_s, 2.0)) / ((R * g_c) * k));
  result.push_back(T);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_14__g_c(double M, double R,
                                                       double T, double k,
                                                       double v_s) {
  std::vector<double> result;
  double g_c = ((M * std::pow(v_s, 2.0)) / ((R * T) * k));
  result.push_back(g_c);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_14__k(double M, double R,
                                                     double T, double g_c,
                                                     double v_s) {
  std::vector<double> result;
  double k = ((M * std::pow(v_s, 2.0)) / ((R * T) * g_c));
  result.push_back(k);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_14__v_s(double M, double R,
                                                       double T, double g_c,
                                                       double k) {
  std::vector<double> result;
  double v_s = std::sqrt(((((R * T) * g_c) * k) / M));
  result.push_back(v_s);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_15__Re(double f) {
  std::vector<double> result;
  double Re = (0.009971220736 / std::pow(f, 4.0));
  result.push_back(Re);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_15__f(double Re) {
  std::vector<double> result;
  double f = (0.316 / std::pow(Re, (1.0 / 4.0)));
  result.push_back(f);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_16__Re(double f) {
  std::vector<double> result;
  double Re = (64.0 / f);
  result.push_back(Re);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_16__f(double Re) {
  std::vector<double> result;
  double f = (64.0 / Re);
  result.push_back(f);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17a__L(double d, double delta_P,
                                                      double mu, double v) {
  std::vector<double> result;
  double L = (((28.9855072463768 * std::pow(d, 2.0)) * delta_P) / (mu * v));
  result.push_back(L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17a__d(double L, double delta_P,
                                                      double mu, double v) {
  std::vector<double> result;
  double d = ((-0.185741756210067) * std::sqrt((((L * mu) * v) / delta_P)));
  result.push_back(d);
  d = (0.185741756210067 * std::sqrt((((L * mu) * v) / delta_P)));
  result.push_back(d);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17a__delta_P(double L, double d,
                                                            double mu,
                                                            double v) {
  std::vector<double> result;
  double delta_P = ((((0.0345 * L) * mu) * v) / std::pow(d, 2.0));
  result.push_back(delta_P);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17a__mu(double L, double d,
                                                       double delta_P,
                                                       double v) {
  std::vector<double> result;
  double mu = (((28.9855072463768 * std::pow(d, 2.0)) * delta_P) / (L * v));
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17a__v(double L, double d,
                                                      double delta_P,
                                                      double mu) {
  std::vector<double> result;
  double v = (((28.9855072463768 * std::pow(d, 2.0)) * delta_P) / (L * mu));
  result.push_back(v);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17b__L(double d, double delta_P,
                                                      double mu, double q) {
  std::vector<double> result;
  double L = (((9.52380952380952 * std::pow(d, 4.0)) * delta_P) / (mu * q));
  result.push_back(L);
  return result;
}

std::vector<std::complex<double>>
FluidFlowVacuumLines_eqn_2_17b__d(double L, double delta_P, double mu,
                                  double q) {
  std::vector<std::complex<double>> result;
  std::complex<double> d =
      (((-0.569242509762222) * std::complex<double>(0.0, 1.0)) *
       std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(d);
  d = ((0.569242509762222 * std::complex<double>(0.0, 1.0)) *
       std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(d);
  d = ((-0.569242509762222) *
       std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(d);
  d = (0.569242509762222 * std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(d);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17b__delta_P(double L, double d,
                                                            double mu,
                                                            double q) {
  std::vector<double> result;
  double delta_P = ((((0.105 * L) * mu) * q) / std::pow(d, 4.0));
  result.push_back(delta_P);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17b__mu(double L, double d,
                                                       double delta_P,
                                                       double q) {
  std::vector<double> result;
  double mu = (((9.52380952380952 * std::pow(d, 4.0)) * delta_P) / (L * q));
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_17b__q(double L, double d,
                                                      double delta_P,
                                                      double mu) {
  std::vector<double> result;
  double q = (((9.52380952380952 * std::pow(d, 4.0)) * delta_P) / (L * mu));
  result.push_back(q);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_18a__D_eq(double R_ll) {
  std::vector<double> result;
  double D_eq = (4.0 * R_ll);
  result.push_back(D_eq);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_18a__R_ll(double D_eq) {
  std::vector<double> result;
  double R_ll = (D_eq / 4.0);
  result.push_back(R_ll);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_18b__R_ll(double h, double w) {
  std::vector<double> result;
  double R_ll = ((h * w) / (2.0 * (h + w)));
  result.push_back(R_ll);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_18b__h(double R_ll, double w) {
  std::vector<double> result;
  double h = (((2.0 * R_ll) * w) / (((-2.0) * R_ll) + w));
  result.push_back(h);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_18b__w(double R_ll, double h) {
  std::vector<double> result;
  double w = (((2.0 * R_ll) * h) / (((-2.0) * R_ll) + h));
  result.push_back(w);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19a__R_ll(double Re, double mu,
                                                         double rho, double v) {
  std::vector<double> result;
  double R_ll = ((Re * mu) / ((4.0 * rho) * v));
  result.push_back(R_ll);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19a__Re(double R_ll, double mu,
                                                       double rho, double v) {
  std::vector<double> result;
  double Re = ((((4.0 * R_ll) * rho) * v) / mu);
  result.push_back(Re);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19a__mu(double R_ll, double Re,
                                                       double rho, double v) {
  std::vector<double> result;
  double mu = ((((4.0 * R_ll) * rho) * v) / Re);
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19a__rho(double R_ll, double Re,
                                                        double mu, double v) {
  std::vector<double> result;
  double rho = ((Re * mu) / ((4.0 * R_ll) * v));
  result.push_back(rho);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19a__v(double R_ll, double Re,
                                                      double mu, double rho) {
  std::vector<double> result;
  double v = ((Re * mu) / ((4.0 * R_ll) * rho));
  result.push_back(v);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19b__Re(double h, double mu,
                                                       double rho, double v,
                                                       double w) {
  std::vector<double> result;
  double Re = (((((2.0 * h) * rho) * v) * w) / (mu * (h + w)));
  result.push_back(Re);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19b__h(double Re, double mu,
                                                      double rho, double v,
                                                      double w) {
  std::vector<double> result;
  double h = (((Re * mu) * w) / (((-Re) * mu) + (((2.0 * rho) * v) * w)));
  result.push_back(h);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19b__mu(double Re, double h,
                                                       double rho, double v,
                                                       double w) {
  std::vector<double> result;
  double mu = (((((2.0 * h) * rho) * v) * w) / (Re * (h + w)));
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19b__rho(double Re, double h,
                                                        double mu, double v,
                                                        double w) {
  std::vector<double> result;
  double rho = (((Re * mu) * (h + w)) / (((2.0 * h) * v) * w));
  result.push_back(rho);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19b__v(double Re, double h,
                                                      double mu, double rho,
                                                      double w) {
  std::vector<double> result;
  double v = (((Re * mu) * (h + w)) / (((2.0 * h) * rho) * w));
  result.push_back(v);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_19b__w(double Re, double h,
                                                      double mu, double rho,
                                                      double v) {
  std::vector<double> result;
  double w = (((Re * h) * mu) / (((-Re) * mu) + (((2.0 * h) * rho) * v)));
  result.push_back(w);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_2__delta(double lambd,
                                                        double psi) {
  std::vector<double> result;
  double delta = ((-0.474424998328794) * std::sqrt((lambd / psi)));
  result.push_back(delta);
  delta = (0.474424998328794 * std::sqrt((lambd / psi)));
  result.push_back(delta);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_2__lambd(double delta,
                                                        double psi) {
  std::vector<double> result;
  double lambd = ((4.44288293815837 * std::pow(delta, 2.0)) * psi);
  result.push_back(lambd);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_2__psi(double delta,
                                                      double lambd) {
  std::vector<double> result;
  double psi = ((0.225079079039277 * lambd) / std::pow(delta, 2.0));
  result.push_back(psi);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_20__L(double sum_equivalent_length,
                                 double sum_pipe) {
  std::vector<double> result;
  double L = (sum_equivalent_length + sum_pipe);
  result.push_back(L);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_20__sum_equivalent_length(double L,
                                                     double sum_pipe) {
  std::vector<double> result;
  double sum_equivalent_length = (L - sum_pipe);
  result.push_back(sum_equivalent_length);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_20__sum_pipe(double L,
                                        double sum_equivalent_length) {
  std::vector<double> result;
  double sum_pipe = (L - sum_equivalent_length);
  result.push_back(sum_pipe);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_22__P_s(double Q_throughput,
                                                       double S_p) {
  std::vector<double> result;
  double P_s = (Q_throughput / S_p);
  result.push_back(P_s);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_22__Q_throughput(double P_s,
                                                                double S_p) {
  std::vector<double> result;
  double Q_throughput = (P_s * S_p);
  result.push_back(Q_throughput);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_22__S_p(double P_s,
                                                       double Q_throughput) {
  std::vector<double> result;
  double S_p = (Q_throughput / P_s);
  result.push_back(S_p);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_25__C(double P_1, double P_2,
                                                     double Q_throughput) {
  std::vector<double> result;
  double C = (Q_throughput / (P_1 - P_2));
  result.push_back(C);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_25__P_1(double C, double P_2,
                                                       double Q_throughput) {
  std::vector<double> result;
  double P_1 = (P_2 + (Q_throughput / C));
  result.push_back(P_1);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_25__P_2(double C, double P_1,
                                                       double Q_throughput) {
  std::vector<double> result;
  double P_2 = (P_1 - (Q_throughput / C));
  result.push_back(P_2);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_25__Q_throughput(double C, double P_1, double P_2) {
  std::vector<double> result;
  double Q_throughput = (C * (P_1 - P_2));
  result.push_back(Q_throughput);
  return result;
}

std::vector<std::complex<double>>
FluidFlowVacuumLines_eqn_2_26__D(double L, double P_downstream, double P_p,
                                 double P_upstream, double mu, double q) {
  std::vector<std::complex<double>> result;
  std::complex<double> D =
      (((-14953.4878122122) * std::complex<double>(0.0, 1.0)) *
       std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) -
                                      (1227184630308510.0 * P_upstream))),
                (1.0 / 4.0)));
  result.push_back(D);
  D = ((14953.4878122122 * std::complex<double>(0.0, 1.0)) *
       std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) -
                                      (1227184630308510.0 * P_upstream))),
                (1.0 / 4.0)));
  result.push_back(D);
  D = ((-14953.4878122122) *
       std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) -
                                      (1227184630308510.0 * P_upstream))),
                (1.0 / 4.0)));
  result.push_back(D);
  D = (14953.4878122122 *
       std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) -
                                      (1227184630308510.0 * P_upstream))),
                (1.0 / 4.0)));
  result.push_back(D);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_26__L(double D, double P_downstream, double P_p,
                                 double P_upstream, double mu, double q) {
  std::vector<double> result;
  double L = (((0.0245436926061703 * std::pow(D, 4.0)) *
               ((-P_downstream) + P_upstream)) /
              (mu * q));
  result.push_back(L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_26__P_downstream(
    double D, double L, double P_p, double P_upstream, double mu, double q) {
  std::vector<double> result;
  double P_downstream =
      (P_upstream - ((((40.7436654315252 * L) * mu) * q) / std::pow(D, 4.0)));
  result.push_back(P_downstream);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_26__P_p(double D, double L,
                                                       double P_downstream,
                                                       double P_upstream,
                                                       double mu, double q) {
  std::vector<double> result;
  double P_p = 0.0;
  result.push_back(P_p);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_26__P_upstream(
    double D, double L, double P_downstream, double P_p, double mu, double q) {
  std::vector<double> result;
  double P_upstream =
      (P_downstream + ((((40.7436654315252 * L) * mu) * q) / std::pow(D, 4.0)));
  result.push_back(P_upstream);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_26__mu(double D, double L, double P_downstream,
                                  double P_p, double P_upstream, double q) {
  std::vector<double> result;
  double mu = (((0.0245436926061703 * std::pow(D, 4.0)) *
                ((-P_downstream) + P_upstream)) /
               (L * q));
  result.push_back(mu);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_26__q(double D, double L, double P_downstream,
                                 double P_p, double P_upstream, double mu) {
  std::vector<double> result;
  double q = (((0.0245436926061703 * std::pow(D, 4.0)) *
               ((-P_downstream) + P_upstream)) /
              (L * mu));
  result.push_back(q);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_28__C(double D, double L,
                                                     double P_p, double mu) {
  std::vector<double> result;
  double C = (((0.0245436926061703 * std::pow(D, 4.0)) * P_p) / (L * mu));
  result.push_back(C);
  return result;
}

std::vector<std::complex<double>>
FluidFlowVacuumLines_eqn_2_28__D(double C, double L, double P_p, double mu) {
  std::vector<std::complex<double>> result;
  std::complex<double> D =
      (((-2.52647511098426) * std::complex<double>(0.0, 1.0)) *
       std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
  result.push_back(D);
  D = ((2.52647511098426 * std::complex<double>(0.0, 1.0)) *
       std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
  result.push_back(D);
  D = ((-2.52647511098426) * std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
  result.push_back(D);
  D = (2.52647511098426 * std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
  result.push_back(D);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_28__L(double C, double D,
                                                     double P_p, double mu) {
  std::vector<double> result;
  double L = (((0.0245436926061703 * std::pow(D, 4.0)) * P_p) / (C * mu));
  result.push_back(L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_28__P_p(double C, double D,
                                                       double L, double mu) {
  std::vector<double> result;
  double P_p = ((((40.7436654315252 * C) * L) * mu) / std::pow(D, 4.0));
  result.push_back(P_p);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_28__mu(double C, double D,
                                                      double L, double P_p) {
  std::vector<double> result;
  double mu = (((0.0245436926061703 * std::pow(D, 4.0)) * P_p) / (C * L));
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_29__C(double S_1, double S_2) {
  std::vector<double> result;
  double C = (((-S_1) * S_2) / (S_1 - S_2));
  result.push_back(C);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_29__S_1(double C, double S_2) {
  std::vector<double> result;
  double S_1 = ((C * S_2) / (C + S_2));
  result.push_back(S_1);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_29__S_2(double C, double S_1) {
  std::vector<double> result;
  double S_2 = ((C * S_1) / (C - S_1));
  result.push_back(S_2);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_3__D(double kn, double lambd) {
  std::vector<double> result;
  double D = (lambd / kn);
  result.push_back(D);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_3__kn(double D, double lambd) {
  std::vector<double> result;
  double kn = (lambd / D);
  result.push_back(kn);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_3__lambd(double D, double kn) {
  std::vector<double> result;
  double lambd = (D * kn);
  result.push_back(lambd);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_31__C(double S_p,
                                                     double S_pump_speed) {
  std::vector<double> result;
  double C = ((S_p * S_pump_speed) / (S_p - S_pump_speed));
  result.push_back(C);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_31__S_p(double C,
                                                       double S_pump_speed) {
  std::vector<double> result;
  double S_p = ((C * S_pump_speed) / (C - S_pump_speed));
  result.push_back(S_p);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_31__S_pump_speed(double C,
                                                                double S_p) {
  std::vector<double> result;
  double S_pump_speed = ((C * S_p) / (C + S_p));
  result.push_back(S_pump_speed);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_32__C_series(double geometric_sum_C) {
  std::vector<double> result;
  double C_series = (1.0 / geometric_sum_C);
  result.push_back(C_series);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_32__geometric_sum_C(double C_series) {
  std::vector<double> result;
  double geometric_sum_C = (1.0 / C_series);
  result.push_back(geometric_sum_C);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_33__C_paralell(double arithmetic_sum_C) {
  std::vector<double> result;
  double C_paralell = (1.0 / arithmetic_sum_C);
  result.push_back(C_paralell);
  return result;
}

std::vector<double>
FluidFlowVacuumLines_eqn_2_33__arithmetic_sum_C(double C_paralell) {
  std::vector<double> result;
  double arithmetic_sum_C = (1.0 / C_paralell);
  result.push_back(arithmetic_sum_C);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_34__C(double C_1, double C_2,
                                                     double D, double L,
                                                     double P_p, double mu) {
  std::vector<double> result;
  double C = ((std::pow(D, 3.0) * (((C_1 * D) * P_p) + (C_2 * mu))) / (L * mu));
  result.push_back(C);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_34__C_1(double C, double C_2,
                                                       double D, double L,
                                                       double P_p, double mu) {
  std::vector<double> result;
  double C_1 =
      ((mu * ((C * L) - (C_2 * std::pow(D, 3.0)))) / (std::pow(D, 4.0) * P_p));
  result.push_back(C_1);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_34__C_2(double C, double C_1,
                                                       double D, double L,
                                                       double P_p, double mu) {
  std::vector<double> result;
  double C_2 = (((C * L) / std::pow(D, 3.0)) - (((C_1 * D) * P_p) / mu));
  result.push_back(C_2);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_34__D(double C, double C_1,
                                                     double C_2, double L,
                                                     double P_p, double mu) {
  throw std::runtime_error("FluidFlowVacuumLines_eqn_2_34__D: requires "
                           "numerical solver (not transpilable)");
}

std::vector<double> FluidFlowVacuumLines_eqn_2_34__L(double C, double C_1,
                                                     double C_2, double D,
                                                     double P_p, double mu) {
  std::vector<double> result;
  double L = ((std::pow(D, 3.0) * (((C_1 * D) * P_p) + (C_2 * mu))) / (C * mu));
  result.push_back(L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_34__P_p(double C, double C_1,
                                                       double C_2, double D,
                                                       double L, double mu) {
  std::vector<double> result;
  double P_p =
      ((mu * ((C * L) - (C_2 * std::pow(D, 3.0)))) / (C_1 * std::pow(D, 4.0)));
  result.push_back(P_p);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_34__mu(double C, double C_1,
                                                      double C_2, double D,
                                                      double L, double P_p) {
  std::vector<double> result;
  double mu =
      (((C_1 * std::pow(D, 4.0)) * P_p) / ((C * L) - (C_2 * std::pow(D, 3.0))));
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_35__C_L(double C_T, double F_p) {
  std::vector<double> result;
  double C_L = (C_T / F_p);
  result.push_back(C_L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_35__C_T(double C_L, double F_p) {
  std::vector<double> result;
  double C_T = (C_L * F_p);
  result.push_back(C_T);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_35__F_p(double C_L, double C_T) {
  std::vector<double> result;
  double F_p = (C_T / C_L);
  result.push_back(F_p);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_36__C(double C_0, double F_t) {
  std::vector<double> result;
  double C = (C_0 * F_t);
  result.push_back(C);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_36__C_0(double C, double F_t) {
  std::vector<double> result;
  double C_0 = (C / F_t);
  result.push_back(C_0);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_36__F_t(double C, double C_0) {
  std::vector<double> result;
  double F_t = (C / C_0);
  result.push_back(F_t);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_37__A(double C, double F_t,
                                                     double M, double T) {
  std::vector<double> result;
  double A = (((0.000681714375311032 * std::pow(C, 2.0)) * M) / (F_t * T));
  result.push_back(A);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_37__C(double A, double F_t,
                                                     double M, double T) {
  std::vector<double> result;
  double C = (38.3 * std::sqrt((((A * F_t) * T) / M)));
  result.push_back(C);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_37__F_t(double A, double C,
                                                       double M, double T) {
  std::vector<double> result;
  double F_t = (((0.000681714375311032 * std::pow(C, 2.0)) * M) / (A * T));
  result.push_back(F_t);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_37__M(double A, double C,
                                                     double F_t, double T) {
  std::vector<double> result;
  double M = ((((1466.89 * A) * F_t) * T) / std::pow(C, 2.0));
  result.push_back(M);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_37__T(double A, double C,
                                                     double F_t, double M) {
  std::vector<double> result;
  double T = (((0.000681714375311032 * std::pow(C, 2.0)) * M) / (A * F_t));
  result.push_back(T);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_4___beta(double mu,
                                                        double vel_grad) {
  std::vector<double> result;
  double _beta = (mu * vel_grad);
  result.push_back(_beta);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_4__mu(double _beta,
                                                     double vel_grad) {
  std::vector<double> result;
  double mu = (_beta / vel_grad);
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_4__vel_grad(double _beta,
                                                           double mu) {
  std::vector<double> result;
  double vel_grad = (_beta / mu);
  result.push_back(vel_grad);
  return result;
}

std::vector<std::complex<double>>
FluidFlowVacuumLines_eqn_2_5__D(double L, double delta_P, double mu, double q) {
  std::vector<std::complex<double>> result;
  std::complex<double> D =
      (((-2.52647511098426) * std::complex<double>(0.0, 1.0)) *
       std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(D);
  D = ((2.52647511098426 * std::complex<double>(0.0, 1.0)) *
       std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(D);
  D = ((-2.52647511098426) * std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(D);
  D = (2.52647511098426 * std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
  result.push_back(D);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_5__L(double D, double delta_P,
                                                    double mu, double q) {
  std::vector<double> result;
  double L = (((0.0245436926061703 * std::pow(D, 4.0)) * delta_P) / (mu * q));
  result.push_back(L);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_5__delta_P(double D, double L,
                                                          double mu, double q) {
  std::vector<double> result;
  double delta_P = ((((40.7436654315252 * L) * mu) * q) / std::pow(D, 4.0));
  result.push_back(delta_P);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_5__mu(double D, double L,
                                                     double delta_P, double q) {
  std::vector<double> result;
  double mu = (((0.0245436926061703 * std::pow(D, 4.0)) * delta_P) / (L * q));
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_5__q(double D, double L,
                                                    double delta_P, double mu) {
  std::vector<double> result;
  double q = (((0.0245436926061703 * std::pow(D, 4.0)) * delta_P) / (L * mu));
  result.push_back(q);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_6__lambd(double mu, double rho,
                                                        double v_a) {
  std::vector<double> result;
  double lambd = ((2.85714285714286 * mu) / (rho * v_a));
  result.push_back(lambd);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_6__mu(double lambd, double rho,
                                                     double v_a) {
  std::vector<double> result;
  double mu = (((0.35 * lambd) * rho) * v_a);
  result.push_back(mu);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_6__rho(double lambd, double mu,
                                                      double v_a) {
  std::vector<double> result;
  double rho = ((2.85714285714286 * mu) / (lambd * v_a));
  result.push_back(rho);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_6__v_a(double lambd, double mu,
                                                      double rho) {
  std::vector<double> result;
  double v_a = ((2.85714285714286 * mu) / (lambd * rho));
  result.push_back(v_a);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_7__T(double k, double m,
                                                    double v_a) {
  std::vector<double> result;
  double T = (((0.392699081698724 * m) * std::pow(v_a, 2.0)) / k);
  result.push_back(T);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_7__k(double T, double m,
                                                    double v_a) {
  std::vector<double> result;
  double k = (((0.392699081698724 * m) * std::pow(v_a, 2.0)) / T);
  result.push_back(k);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_7__m(double T, double k,
                                                    double v_a) {
  std::vector<double> result;
  double m = (((2.54647908947033 * T) * k) / std::pow(v_a, 2.0));
  result.push_back(m);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_7__v_a(double T, double k,
                                                      double m) {
  std::vector<double> result;
  double v_a = (1.59576912160573 * std::sqrt(((T * k) / m)));
  result.push_back(v_a);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_8__M(double P_c, double T_c,
                                                    double mu_c) {
  std::vector<double> result;
  double M = (((0.0168662506324844 * std::pow(T_c, (1.0 / 3.0))) *
               std::pow(mu_c, 2.0)) /
              std::pow(P_c, (4.0 / 3.0)));
  result.push_back(M);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_8__P_c(double M, double T_c,
                                                      double mu_c) {
  std::vector<double> result;
  double P_c =
      ((-0.046801946114055) *
       std::pow(((std::pow(T_c, 0.166666666666667) * mu_c) / std::pow(M, 0.5)),
                (3.0 / 2.0)));
  result.push_back(P_c);
  P_c =
      (0.046801946114055 *
       std::pow(((std::pow(T_c, 0.166666666666667) * mu_c) / std::pow(M, 0.5)),
                (3.0 / 2.0)));
  result.push_back(P_c);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_8__T_c(double M, double P_c,
                                                      double mu_c) {
  std::vector<double> result;
  double T_c = (((208422.380089 * std::pow(M, 3.0)) * std::pow(P_c, 4.0)) /
                std::pow(mu_c, 6.0));
  result.push_back(T_c);
  return result;
}

std::vector<double> FluidFlowVacuumLines_eqn_2_8__mu_c(double M, double P_c,
                                                       double T_c) {
  std::vector<double> result;
  double mu_c = (((7.7 * std::sqrt(M)) * std::pow(P_c, (2.0 / 3.0))) /
                 std::pow(T_c, (1.0 / 6.0)));
  result.push_back(mu_c);
  return result;
}

#endif // VAKYUME_FLUIDFLOWVACUUMLINES_HPP
