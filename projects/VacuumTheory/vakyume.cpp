#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

// Lambert W₀ (principal branch) via Halley iteration
inline double _lambertW(double x) {
    if (x == 0.0) return 0.0;
    double w = (x > 1.0) ? std::log(x) - std::log(std::log(x)) : 0.5;
    for (int i = 0; i < 50; ++i) {
        double ew = std::exp(w);
        double wew = w * ew;
        double f  = wew - x;
        double fp = ew * (w + 1.0);
        double d  = f / (fp - (w + 2.0) * f / (2.0 * (w + 1.0)));
        w -= d;
        if (std::abs(d) < 1e-15) break;
    }
    return w;
}

std::vector<double> AirLeak_eqn_4_10__T(double V, double del_P, double leakage, double t) {
    std::vector<double> result;
    double T = (((3.127 * V) * del_P) / (leakage * t));
    result.push_back(T);
    return result;
}
std::vector<double> AirLeak_eqn_4_10__V(double T, double del_P, double leakage, double t) {
    std::vector<double> result;
    double V = ((((0.319795330988168 * T) * leakage) * t) / del_P);
    result.push_back(V);
    return result;
}
std::vector<double> AirLeak_eqn_4_10__del_P(double T, double V, double leakage, double t) {
    std::vector<double> result;
    double del_P = ((((0.319795330988168 * T) * leakage) * t) / V);
    result.push_back(del_P);
    return result;
}
std::vector<double> AirLeak_eqn_4_10__leakage(double T, double V, double del_P, double t) {
    std::vector<double> result;
    double leakage = (((3.127 * V) * del_P) / (T * t));
    result.push_back(leakage);
    return result;
}
std::vector<double> AirLeak_eqn_4_10__t(double T, double V, double del_P, double leakage) {
    std::vector<double> result;
    double t = (((3.127 * V) * del_P) / (T * leakage));
    result.push_back(t);
    return result;
}
std::vector<double> AirLeak_eqn_4_7__W(double W_T, double sum_individual_leak_rates) {
    std::vector<double> result;
    double W = (W_T - sum_individual_leak_rates);
    result.push_back(W);
    return result;
}
std::vector<double> AirLeak_eqn_4_7__W_T(double W, double sum_individual_leak_rates) {
    std::vector<double> result;
    double W_T = (W + sum_individual_leak_rates);
    result.push_back(W_T);
    return result;
}
std::vector<double> AirLeak_eqn_4_7__sum_individual_leak_rates(double W, double W_T) {
    std::vector<double> result;
    double sum_individual_leak_rates = ((-W) + W_T);
    result.push_back(sum_individual_leak_rates);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_1__D(double Re, double mu, double rho, double v) {
    std::vector<double> result;
    double D = ((Re * mu) / (rho * v));
    result.push_back(D);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_1__Re(double D, double mu, double rho, double v) {
    std::vector<double> result;
    double Re = (((D * rho) * v) / mu);
    result.push_back(Re);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_1__mu(double D, double Re, double rho, double v) {
    std::vector<double> result;
    double mu = (((D * rho) * v) / Re);
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_1__rho(double D, double Re, double mu, double v) {
    std::vector<double> result;
    double rho = ((Re * mu) / (D * v));
    result.push_back(rho);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_1__v(double D, double Re, double mu, double rho) {
    std::vector<double> result;
    double v = ((Re * mu) / (D * rho));
    result.push_back(v);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_10__Suc_Pres(double delta_P, double oper_press) {
    std::vector<double> result;
    double Suc_Pres = ((-delta_P) + oper_press);
    result.push_back(Suc_Pres);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_10__delta_P(double Suc_Pres, double oper_press) {
    std::vector<double> result;
    double delta_P = ((-Suc_Pres) + oper_press);
    result.push_back(delta_P);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_10__oper_press(double Suc_Pres, double delta_P) {
    std::vector<double> result;
    double oper_press = (Suc_Pres + delta_P);
    result.push_back(oper_press);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_11__D(double L, double f, double g_c, double h_r, double v) {
    std::vector<double> result;
    double D = (((L * f) * std::pow(v, 2.0)) / ((2.0 * g_c) * h_r));
    result.push_back(D);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_11__L(double D, double f, double g_c, double h_r, double v) {
    std::vector<double> result;
    double L = ((((2.0 * D) * g_c) * h_r) / (f * std::pow(v, 2.0)));
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_11__f(double D, double L, double g_c, double h_r, double v) {
    std::vector<double> result;
    double f = ((((2.0 * D) * g_c) * h_r) / (L * std::pow(v, 2.0)));
    result.push_back(f);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_11__g_c(double D, double L, double f, double h_r, double v) {
    std::vector<double> result;
    double g_c = (((L * f) * std::pow(v, 2.0)) / ((2.0 * D) * h_r));
    result.push_back(g_c);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_11__h_r(double D, double L, double f, double g_c, double v) {
    std::vector<double> result;
    double h_r = (((L * f) * std::pow(v, 2.0)) / ((2.0 * D) * g_c));
    result.push_back(h_r);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_11__v(double D, double L, double f, double g_c, double h_r) {
    std::vector<double> result;
    double v = ((-std::sqrt(2.0)) * std::sqrt((((D * g_c) * h_r) / (L * f))));
    result.push_back(v);
    v = (std::sqrt(2.0) * std::sqrt((((D * g_c) * h_r) / (L * f))));
    result.push_back(v);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_12__L(double d, double delta_P, double f, double g, double rho, double v) {
    std::vector<double> result;
    double L = ((((0.464037122969838 * d) * delta_P) * g) / ((f * rho) * std::pow(v, 2.0)));
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_12__d(double L, double delta_P, double f, double g, double rho, double v) {
    std::vector<double> result;
    double d = (((((2.155 * L) * f) * rho) * std::pow(v, 2.0)) / (delta_P * g));
    result.push_back(d);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_12__delta_P(double L, double d, double f, double g, double rho, double v) {
    std::vector<double> result;
    double delta_P = (((((2.155 * L) * f) * rho) * std::pow(v, 2.0)) / (d * g));
    result.push_back(delta_P);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_12__f(double L, double d, double delta_P, double g, double rho, double v) {
    std::vector<double> result;
    double f = ((((0.464037122969838 * d) * delta_P) * g) / ((L * rho) * std::pow(v, 2.0)));
    result.push_back(f);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_12__g(double L, double d, double delta_P, double f, double rho, double v) {
    std::vector<double> result;
    double g = (((((2.155 * L) * f) * rho) * std::pow(v, 2.0)) / (d * delta_P));
    result.push_back(g);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_12__rho(double L, double d, double delta_P, double f, double g, double v) {
    std::vector<double> result;
    double rho = ((((0.464037122969838 * d) * delta_P) * g) / ((L * f) * std::pow(v, 2.0)));
    result.push_back(rho);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_12__v(double L, double d, double delta_P, double f, double g, double rho) {
    std::vector<double> result;
    double v = ((-0.681202703290172) * std::sqrt((((d * delta_P) * g) / ((L * f) * rho))));
    result.push_back(v);
    v = (0.681202703290172 * std::sqrt((((d * delta_P) * g) / ((L * f) * rho))));
    result.push_back(v);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_13__L(double d, double delta_P, double f, double q, double rho) {
    std::vector<double> result;
    double L = (((0.465116279069767 * std::pow(d, 5.0)) * delta_P) / ((f * std::pow(q, 2.0)) * rho));
    result.push_back(L);
    return result;
}
std::vector<std::complex<double>> FluidFlowVacuumLines_eqn_2_13__d(double L, double delta_P, double f, double q, double rho) {
    std::vector<std::complex<double>> result;
    std::complex<double> d = (1.16543402167043 * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0)));
    result.push_back(d);
    d = (((-0.942855929354115) * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))) - ((0.685024930457783 * std::complex<double>(0.0, 1.0)) * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))));
    result.push_back(d);
    d = (((-0.942855929354115) * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))) + ((0.685024930457783 * std::complex<double>(0.0, 1.0)) * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))));
    result.push_back(d);
    d = ((0.360138918518902 * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))) - ((1.10839362062173 * std::complex<double>(0.0, 1.0)) * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))));
    result.push_back(d);
    d = ((0.360138918518902 * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))) + ((1.10839362062173 * std::complex<double>(0.0, 1.0)) * std::pow(((((L * f) * std::pow(q, 2.0)) * rho) / delta_P), (1.0 / 5.0))));
    result.push_back(d);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_13__delta_P(double L, double d, double f, double q, double rho) {
    std::vector<double> result;
    double delta_P = (((((2.15 * L) * f) * std::pow(q, 2.0)) * rho) / std::pow(d, 5.0));
    result.push_back(delta_P);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_13__f(double L, double d, double delta_P, double q, double rho) {
    std::vector<double> result;
    double f = (((0.465116279069767 * std::pow(d, 5.0)) * delta_P) / ((L * std::pow(q, 2.0)) * rho));
    result.push_back(f);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_13__q(double L, double d, double delta_P, double f, double rho) {
    std::vector<double> result;
    double q = ((-0.681994339470473) * std::sqrt(((std::pow(d, 5.0) * delta_P) / ((L * f) * rho))));
    result.push_back(q);
    q = (0.681994339470473 * std::sqrt(((std::pow(d, 5.0) * delta_P) / ((L * f) * rho))));
    result.push_back(q);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_13__rho(double L, double d, double delta_P, double f, double q) {
    std::vector<double> result;
    double rho = (((0.465116279069767 * std::pow(d, 5.0)) * delta_P) / ((L * f) * std::pow(q, 2.0)));
    result.push_back(rho);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_14__M(double R, double T, double g_c, double k, double v_s) {
    std::vector<double> result;
    double M = ((((R * T) * g_c) * k) / std::pow(v_s, 2.0));
    result.push_back(M);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_14__R(double M, double T, double g_c, double k, double v_s) {
    std::vector<double> result;
    double R = ((M * std::pow(v_s, 2.0)) / ((T * g_c) * k));
    result.push_back(R);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_14__T(double M, double R, double g_c, double k, double v_s) {
    std::vector<double> result;
    double T = ((M * std::pow(v_s, 2.0)) / ((R * g_c) * k));
    result.push_back(T);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_14__g_c(double M, double R, double T, double k, double v_s) {
    std::vector<double> result;
    double g_c = ((M * std::pow(v_s, 2.0)) / ((R * T) * k));
    result.push_back(g_c);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_14__k(double M, double R, double T, double g_c, double v_s) {
    std::vector<double> result;
    double k = ((M * std::pow(v_s, 2.0)) / ((R * T) * g_c));
    result.push_back(k);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_14__v_s(double M, double R, double T, double g_c, double k) {
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
std::vector<double> FluidFlowVacuumLines_eqn_2_17__L(double d, double delta_P, double mu, double v) {
    std::vector<double> result;
    double L = (((28.9855072463768 * std::pow(d, 2.0)) * delta_P) / (mu * v));
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_17__d(double L, double delta_P, double mu, double v) {
    std::vector<double> result;
    double d = ((-0.185741756210067) * std::sqrt((((L * mu) * v) / delta_P)));
    result.push_back(d);
    d = (0.185741756210067 * std::sqrt((((L * mu) * v) / delta_P)));
    result.push_back(d);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_17__delta_P(double L, double d, double mu, double v) {
    std::vector<double> result;
    double delta_P = ((((0.0345 * L) * mu) * v) / std::pow(d, 2.0));
    result.push_back(delta_P);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_17__mu(double L, double d, double delta_P, double v) {
    std::vector<double> result;
    double mu = (((28.9855072463768 * std::pow(d, 2.0)) * delta_P) / (L * v));
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_17__q(double L, double d, double delta_P, double mu) {
    std::vector<double> result;
    double q = (((9.52380952380952 * std::pow(d, 4.0)) * delta_P) / (L * mu));
    result.push_back(q);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_17__v(double L, double d, double delta_P, double mu) {
    std::vector<double> result;
    double v = (((28.9855072463768 * std::pow(d, 2.0)) * delta_P) / (L * mu));
    result.push_back(v);
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
std::vector<double> FluidFlowVacuumLines_eqn_2_19a__R_ll(double Re, double mu, double rho, double v) {
    std::vector<double> result;
    double R_ll = ((Re * mu) / ((4.0 * rho) * v));
    result.push_back(R_ll);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19a__Re(double R_ll, double mu, double rho, double v) {
    std::vector<double> result;
    double Re = ((((4.0 * R_ll) * rho) * v) / mu);
    result.push_back(Re);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19a__mu(double R_ll, double Re, double rho, double v) {
    std::vector<double> result;
    double mu = ((((4.0 * R_ll) * rho) * v) / Re);
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19a__rho(double R_ll, double Re, double mu, double v) {
    std::vector<double> result;
    double rho = ((Re * mu) / ((4.0 * R_ll) * v));
    result.push_back(rho);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19a__v(double R_ll, double Re, double mu, double rho) {
    std::vector<double> result;
    double v = ((Re * mu) / ((4.0 * R_ll) * rho));
    result.push_back(v);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19b__Re(double h, double mu, double rho, double v, double w) {
    std::vector<double> result;
    double Re = (((((2.0 * h) * rho) * v) * w) / (mu * (h + w)));
    result.push_back(Re);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19b__h(double Re, double mu, double rho, double v, double w) {
    std::vector<double> result;
    double h = (((Re * mu) * w) / (((-Re) * mu) + (((2.0 * rho) * v) * w)));
    result.push_back(h);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19b__mu(double Re, double h, double rho, double v, double w) {
    std::vector<double> result;
    double mu = (((((2.0 * h) * rho) * v) * w) / (Re * (h + w)));
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19b__rho(double Re, double h, double mu, double v, double w) {
    std::vector<double> result;
    double rho = (((Re * mu) * (h + w)) / (((2.0 * h) * v) * w));
    result.push_back(rho);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19b__v(double Re, double h, double mu, double rho, double w) {
    std::vector<double> result;
    double v = (((Re * mu) * (h + w)) / (((2.0 * h) * rho) * w));
    result.push_back(v);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_19b__w(double Re, double h, double mu, double rho, double v) {
    std::vector<double> result;
    double w = (((Re * h) * mu) / (((-Re) * mu) + (((2.0 * h) * rho) * v)));
    result.push_back(w);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_2__delta(double lambd, double psi) {
    std::vector<double> result;
    double delta = ((-0.474424998328794) * std::sqrt((lambd / psi)));
    result.push_back(delta);
    delta = (0.474424998328794 * std::sqrt((lambd / psi)));
    result.push_back(delta);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_2__lambd(double delta, double psi) {
    std::vector<double> result;
    double lambd = ((4.44288293815837 * std::pow(delta, 2.0)) * psi);
    result.push_back(lambd);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_2__psi(double delta, double lambd) {
    std::vector<double> result;
    double psi = ((0.225079079039277 * lambd) / std::pow(delta, 2.0));
    result.push_back(psi);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_20__L(double sum_equivalent_length, double sum_pipe) {
    std::vector<double> result;
    double L = (sum_equivalent_length + sum_pipe);
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_20__sum_equivalent_length(double L, double sum_pipe) {
    std::vector<double> result;
    double sum_equivalent_length = (L - sum_pipe);
    result.push_back(sum_equivalent_length);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_20__sum_pipe(double L, double sum_equivalent_length) {
    std::vector<double> result;
    double sum_pipe = (L - sum_equivalent_length);
    result.push_back(sum_pipe);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_22__P_s(double Q_throughput, double S_p) {
    std::vector<double> result;
    double P_s = (Q_throughput / S_p);
    result.push_back(P_s);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_22__Q_throughput(double P_s, double S_p) {
    std::vector<double> result;
    double Q_throughput = (P_s * S_p);
    result.push_back(Q_throughput);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_22__S_p(double P_s, double Q_throughput) {
    std::vector<double> result;
    double S_p = (Q_throughput / P_s);
    result.push_back(S_p);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_25__C(double P_1, double P_2, double Q_throughput) {
    std::vector<double> result;
    double C = (Q_throughput / (P_1 - P_2));
    result.push_back(C);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_25__P_1(double C, double P_2, double Q_throughput) {
    std::vector<double> result;
    double P_1 = (P_2 + (Q_throughput / C));
    result.push_back(P_1);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_25__P_2(double C, double P_1, double Q_throughput) {
    std::vector<double> result;
    double P_2 = (P_1 - (Q_throughput / C));
    result.push_back(P_2);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_25__Q_throughput(double C, double P_1, double P_2) {
    std::vector<double> result;
    double Q_throughput = (C * (P_1 - P_2));
    result.push_back(Q_throughput);
    return result;
}
std::vector<std::complex<double>> FluidFlowVacuumLines_eqn_2_26__D(double L, double P_downstream, double P_p, double P_upstream, double mu, double q) {
    std::vector<std::complex<double>> result;
    std::complex<double> D = (((-14953.4878122122) * std::complex<double>(0.0, 1.0)) * std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) - (1227184630308510.0 * P_upstream))), (1.0 / 4.0)));
    result.push_back(D);
    D = ((14953.4878122122 * std::complex<double>(0.0, 1.0)) * std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) - (1227184630308510.0 * P_upstream))), (1.0 / 4.0)));
    result.push_back(D);
    D = ((-14953.4878122122) * std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) - (1227184630308510.0 * P_upstream))), (1.0 / 4.0)));
    result.push_back(D);
    D = (14953.4878122122 * std::pow(((((-L) * mu) * q) / ((1227184630308510.0 * P_downstream) - (1227184630308510.0 * P_upstream))), (1.0 / 4.0)));
    result.push_back(D);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_26__L(double D, double P_downstream, double P_p, double P_upstream, double mu, double q) {
    std::vector<double> result;
    double L = (((0.0245436926061703 * std::pow(D, 4.0)) * ((-P_downstream) + P_upstream)) / (mu * q));
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_26__P_downstream(double D, double L, double P_p, double P_upstream, double mu, double q) {
    std::vector<double> result;
    double P_downstream = (P_upstream - ((((40.7436654315252 * L) * mu) * q) / std::pow(D, 4.0)));
    result.push_back(P_downstream);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_26__P_p(double D, double L, double P_downstream, double P_upstream, double mu, double q) {
    std::vector<double> result;
    double P_p = 0.0;
    result.push_back(P_p);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_26__P_upstream(double D, double L, double P_downstream, double P_p, double mu, double q) {
    std::vector<double> result;
    double P_upstream = (P_downstream + ((((40.7436654315252 * L) * mu) * q) / std::pow(D, 4.0)));
    result.push_back(P_upstream);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_26__mu(double D, double L, double P_downstream, double P_p, double P_upstream, double q) {
    std::vector<double> result;
    double mu = (((0.0245436926061703 * std::pow(D, 4.0)) * ((-P_downstream) + P_upstream)) / (L * q));
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_26__q(double D, double L, double P_downstream, double P_p, double P_upstream, double mu) {
    std::vector<double> result;
    double q = (((0.0245436926061703 * std::pow(D, 4.0)) * ((-P_downstream) + P_upstream)) / (L * mu));
    result.push_back(q);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_28__C(double D, double L, double P_p, double mu) {
    std::vector<double> result;
    double C = (((0.0245436926061703 * std::pow(D, 4.0)) * P_p) / (L * mu));
    result.push_back(C);
    return result;
}
std::vector<std::complex<double>> FluidFlowVacuumLines_eqn_2_28__D(double C, double L, double P_p, double mu) {
    std::vector<std::complex<double>> result;
    std::complex<double> D = (((-2.52647511098426) * std::complex<double>(0.0, 1.0)) * std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
    result.push_back(D);
    D = ((2.52647511098426 * std::complex<double>(0.0, 1.0)) * std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
    result.push_back(D);
    D = ((-2.52647511098426) * std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
    result.push_back(D);
    D = (2.52647511098426 * std::pow((((C * L) * mu) / P_p), (1.0 / 4.0)));
    result.push_back(D);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_28__L(double C, double D, double P_p, double mu) {
    std::vector<double> result;
    double L = (((0.0245436926061703 * std::pow(D, 4.0)) * P_p) / (C * mu));
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_28__P_p(double C, double D, double L, double mu) {
    std::vector<double> result;
    double P_p = ((((40.7436654315252 * C) * L) * mu) / std::pow(D, 4.0));
    result.push_back(P_p);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_28__mu(double C, double D, double L, double P_p) {
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
std::vector<double> FluidFlowVacuumLines_eqn_2_31__C(double S_p, double S_pump_speed) {
    std::vector<double> result;
    double C = ((S_p * S_pump_speed) / (S_p - S_pump_speed));
    result.push_back(C);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_31__S_p(double C, double S_pump_speed) {
    std::vector<double> result;
    double S_p = ((C * S_pump_speed) / (C - S_pump_speed));
    result.push_back(S_p);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_31__S_pump_speed(double C, double S_p) {
    std::vector<double> result;
    double S_pump_speed = ((C * S_p) / (C + S_p));
    result.push_back(S_pump_speed);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_32__C_series(double geometric_sum_C) {
    std::vector<double> result;
    double C_series = (1.0 / geometric_sum_C);
    result.push_back(C_series);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_32__geometric_sum_C(double C_series) {
    std::vector<double> result;
    double geometric_sum_C = (1.0 / C_series);
    result.push_back(geometric_sum_C);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_33__C_paralell(double arithmetic_sum_C) {
    std::vector<double> result;
    double C_paralell = (1.0 / arithmetic_sum_C);
    result.push_back(C_paralell);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_33__arithmetic_sum_C(double C_paralell) {
    std::vector<double> result;
    double arithmetic_sum_C = (1.0 / C_paralell);
    result.push_back(arithmetic_sum_C);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_34__C(double C_1, double C_2, double D, double L, double P_p, double mu) {
    std::vector<double> result;
    double C = ((std::pow(D, 3.0) * (((C_1 * D) * P_p) + (C_2 * mu))) / (L * mu));
    result.push_back(C);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_34__C_1(double C, double C_2, double D, double L, double P_p, double mu) {
    std::vector<double> result;
    double C_1 = ((mu * ((C * L) - (C_2 * std::pow(D, 3.0)))) / (std::pow(D, 4.0) * P_p));
    result.push_back(C_1);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_34__C_2(double C, double C_1, double D, double L, double P_p, double mu) {
    std::vector<double> result;
    double C_2 = (((C * L) / std::pow(D, 3.0)) - (((C_1 * D) * P_p) / mu));
    result.push_back(C_2);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_34__D(double C, double C_1, double C_2, double L, double P_p, double mu) {
    throw std::runtime_error("FluidFlowVacuumLines_eqn_2_34__D: requires numerical solver (not transpilable)");
}
std::vector<double> FluidFlowVacuumLines_eqn_2_34__L(double C, double C_1, double C_2, double D, double P_p, double mu) {
    std::vector<double> result;
    double L = ((std::pow(D, 3.0) * (((C_1 * D) * P_p) + (C_2 * mu))) / (C * mu));
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_34__P_p(double C, double C_1, double C_2, double D, double L, double mu) {
    std::vector<double> result;
    double P_p = ((mu * ((C * L) - (C_2 * std::pow(D, 3.0)))) / (C_1 * std::pow(D, 4.0)));
    result.push_back(P_p);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_34__mu(double C, double C_1, double C_2, double D, double L, double P_p) {
    std::vector<double> result;
    double mu = (((C_1 * std::pow(D, 4.0)) * P_p) / ((C * L) - (C_2 * std::pow(D, 3.0))));
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
std::vector<double> FluidFlowVacuumLines_eqn_2_37__A(double C, double F_t, double M, double T) {
    std::vector<double> result;
    double A = (((0.000681714375311032 * std::pow(C, 2.0)) * M) / (F_t * T));
    result.push_back(A);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_37__C(double A, double F_t, double M, double T) {
    std::vector<double> result;
    double C = (38.3 * std::sqrt((((A * F_t) * T) / M)));
    result.push_back(C);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_37__F_t(double A, double C, double M, double T) {
    std::vector<double> result;
    double F_t = (((0.000681714375311032 * std::pow(C, 2.0)) * M) / (A * T));
    result.push_back(F_t);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_37__M(double A, double C, double F_t, double T) {
    std::vector<double> result;
    double M = ((((1466.89 * A) * F_t) * T) / std::pow(C, 2.0));
    result.push_back(M);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_37__T(double A, double C, double F_t, double M) {
    std::vector<double> result;
    double T = (((0.000681714375311032 * std::pow(C, 2.0)) * M) / (A * F_t));
    result.push_back(T);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_4___beta(double mu, double vel_grad) {
    std::vector<double> result;
    double _beta = (mu * vel_grad);
    result.push_back(_beta);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_4__mu(double _beta, double vel_grad) {
    std::vector<double> result;
    double mu = (_beta / vel_grad);
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_4__vel_grad(double _beta, double mu) {
    std::vector<double> result;
    double vel_grad = (_beta / mu);
    result.push_back(vel_grad);
    return result;
}
std::vector<std::complex<double>> FluidFlowVacuumLines_eqn_2_5__D(double L, double delta_P, double mu, double q) {
    std::vector<std::complex<double>> result;
    std::complex<double> D = (((-2.52647511098426) * std::complex<double>(0.0, 1.0)) * std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
    result.push_back(D);
    D = ((2.52647511098426 * std::complex<double>(0.0, 1.0)) * std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
    result.push_back(D);
    D = ((-2.52647511098426) * std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
    result.push_back(D);
    D = (2.52647511098426 * std::pow((((L * mu) * q) / delta_P), (1.0 / 4.0)));
    result.push_back(D);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_5__L(double D, double delta_P, double mu, double q) {
    std::vector<double> result;
    double L = (((0.0245436926061703 * std::pow(D, 4.0)) * delta_P) / (mu * q));
    result.push_back(L);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_5__delta_P(double D, double L, double mu, double q) {
    std::vector<double> result;
    double delta_P = ((((40.7436654315252 * L) * mu) * q) / std::pow(D, 4.0));
    result.push_back(delta_P);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_5__mu(double D, double L, double delta_P, double q) {
    std::vector<double> result;
    double mu = (((0.0245436926061703 * std::pow(D, 4.0)) * delta_P) / (L * q));
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_5__q(double D, double L, double delta_P, double mu) {
    std::vector<double> result;
    double q = (((0.0245436926061703 * std::pow(D, 4.0)) * delta_P) / (L * mu));
    result.push_back(q);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_6__lambd(double mu, double rho, double v_a) {
    std::vector<double> result;
    double lambd = ((2.85714285714286 * mu) / (rho * v_a));
    result.push_back(lambd);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_6__mu(double lambd, double rho, double v_a) {
    std::vector<double> result;
    double mu = (((0.35 * lambd) * rho) * v_a);
    result.push_back(mu);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_6__rho(double lambd, double mu, double v_a) {
    std::vector<double> result;
    double rho = ((2.85714285714286 * mu) / (lambd * v_a));
    result.push_back(rho);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_6__v_a(double lambd, double mu, double rho) {
    std::vector<double> result;
    double v_a = ((2.85714285714286 * mu) / (lambd * rho));
    result.push_back(v_a);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_7__T(double k, double m, double v_a) {
    std::vector<double> result;
    double T = (((0.392699081698724 * m) * std::pow(v_a, 2.0)) / k);
    result.push_back(T);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_7__k(double T, double m, double v_a) {
    std::vector<double> result;
    double k = (((0.392699081698724 * m) * std::pow(v_a, 2.0)) / T);
    result.push_back(k);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_7__m(double T, double k, double v_a) {
    std::vector<double> result;
    double m = (((2.54647908947033 * T) * k) / std::pow(v_a, 2.0));
    result.push_back(m);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_7__v_a(double T, double k, double m) {
    std::vector<double> result;
    double v_a = (1.59576912160573 * std::sqrt(((T * k) / m)));
    result.push_back(v_a);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_8__M(double P_c, double T_c, double mu_c) {
    std::vector<double> result;
    double M = (((0.0168662506324844 * std::pow(T_c, (1.0 / 3.0))) * std::pow(mu_c, 2.0)) / std::pow(P_c, (4.0 / 3.0)));
    result.push_back(M);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_8__P_c(double M, double T_c, double mu_c) {
    std::vector<double> result;
    double P_c = ((-0.046801946114055) * std::pow(((std::pow(T_c, 0.166666666666667) * mu_c) / std::pow(M, 0.5)), (3.0 / 2.0)));
    result.push_back(P_c);
    P_c = (0.046801946114055 * std::pow(((std::pow(T_c, 0.166666666666667) * mu_c) / std::pow(M, 0.5)), (3.0 / 2.0)));
    result.push_back(P_c);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_8__T_c(double M, double P_c, double mu_c) {
    std::vector<double> result;
    double T_c = (((208422.380089 * std::pow(M, 3.0)) * std::pow(P_c, 4.0)) / std::pow(mu_c, 6.0));
    result.push_back(T_c);
    return result;
}
std::vector<double> FluidFlowVacuumLines_eqn_2_8__mu_c(double M, double P_c, double T_c) {
    std::vector<double> result;
    double mu_c = (((7.7 * std::sqrt(M)) * std::pow(P_c, (2.0 / 3.0))) / std::pow(T_c, (1.0 / 6.0)));
    result.push_back(mu_c);
    return result;
}
std::vector<double> LiquidRing_eqn_10_1__D_r(double sig_R, double w) {
    std::vector<double> result;
    double D_r = ((229.357798165138 * sig_R) / w);
    result.push_back(D_r);
    return result;
}
std::vector<double> LiquidRing_eqn_10_1__sig_R(double D_r, double w) {
    std::vector<double> result;
    double sig_R = ((0.00436 * D_r) * w);
    result.push_back(sig_R);
    return result;
}
std::vector<double> LiquidRing_eqn_10_1__w(double D_r, double sig_R) {
    std::vector<double> result;
    double w = ((229.357798165138 * sig_R) / D_r);
    result.push_back(w);
    return result;
}
std::vector<double> LiquidRing_eqn_10_10__bhp(double bhp_0, double mu, double rho) {
    std::vector<double> result;
    double bhp = ((0.0005 * bhp_0) * (((31.0 * std::pow(mu, (4.0 / 25.0))) * std::pow(rho, (21.0 / 25.0))) + 1000.0));
    result.push_back(bhp);
    return result;
}
std::vector<double> LiquidRing_eqn_10_10__bhp_0(double bhp, double mu, double rho) {
    std::vector<double> result;
    double bhp_0 = ((2000.0 * bhp) / (((31.0 * std::pow(mu, 0.16)) * std::pow(rho, 0.84)) + 1000.0));
    result.push_back(bhp_0);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_10__mu(double bhp, double bhp_0, double rho) {
    std::vector<std::complex<double>> result;
    std::complex<double> mu = (((-4.7751763343393e-10) * std::complex<double>(0.0, 1.0)) * std::pow((((2000.0 * bhp) / (bhp_0 * std::pow(rho, 0.84))) - (1000.0 / std::pow(rho, 0.84))), (25.0 / 4.0)));
    result.push_back(mu);
    mu = ((4.7751763343393e-10 * std::complex<double>(0.0, 1.0)) * std::pow((((2000.0 * bhp) / (bhp_0 * std::pow(rho, 0.84))) - (1000.0 / std::pow(rho, 0.84))), (25.0 / 4.0)));
    result.push_back(mu);
    mu = ((-4.7751763343393e-10) * std::pow((((2000.0 * bhp) / (bhp_0 * std::pow(rho, 0.84))) - (1000.0 / std::pow(rho, 0.84))), (25.0 / 4.0)));
    result.push_back(mu);
    mu = (4.7751763343393e-10 * std::pow((((2000.0 * bhp) / (bhp_0 * std::pow(rho, 0.84))) - (1000.0 / std::pow(rho, 0.84))), (25.0 / 4.0)));
    result.push_back(mu);
    return result;
}
std::vector<double> LiquidRing_eqn_10_10__rho(double bhp, double bhp_0, double mu) {
    std::vector<double> result;
    double rho = std::pow((((bhp / bhp_0) - 0.5) / (0.0155 * std::pow(mu, 0.16))), (1.0 / 0.84));
    return {rho};
}
std::vector<double> LiquidRing_eqn_10_11__T_c(double T_s) {
    std::vector<double> result;
    double T_c = (T_s + 10.0);
    result.push_back(T_c);
    return result;
}
std::vector<double> LiquidRing_eqn_10_11__T_s(double T_c) {
    std::vector<double> result;
    double T_s = (T_c - 10.0);
    result.push_back(T_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_12__T_c(double T_s) {
    std::vector<double> result;
    double T_c = (T_s + 5.0);
    result.push_back(T_c);
    return result;
}
std::vector<double> LiquidRing_eqn_10_12__T_s(double T_c) {
    std::vector<double> result;
    double T_s = (T_c - 5.0);
    result.push_back(T_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_13__T_c(double T_s) {
    std::vector<double> result;
    double T_c = (T_s + 25.0);
    result.push_back(T_c);
    return result;
}
std::vector<double> LiquidRing_eqn_10_13__T_s(double T_c) {
    std::vector<double> result;
    double T_s = (T_c - 25.0);
    result.push_back(T_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_14__T_c(double T_s) {
    std::vector<double> result;
    double T_c = (T_s + 12.0);
    result.push_back(T_c);
    return result;
}
std::vector<double> LiquidRing_eqn_10_14__T_s(double T_c) {
    std::vector<double> result;
    double T_s = (T_c - 12.0);
    result.push_back(T_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_15__P(double S_Th, double S_p, double p_s) {
    std::vector<double> result;
    double P = ((S_Th * p_s) / (S_Th - S_p));
    result.push_back(P);
    return result;
}
std::vector<double> LiquidRing_eqn_10_15__S_Th(double P, double S_p, double p_s) {
    std::vector<double> result;
    double S_Th = ((P * S_p) / (P - p_s));
    result.push_back(S_Th);
    return result;
}
std::vector<double> LiquidRing_eqn_10_15__S_p(double P, double S_Th, double p_s) {
    std::vector<double> result;
    double S_p = ((S_Th * (P - p_s)) / P);
    result.push_back(S_p);
    return result;
}
std::vector<double> LiquidRing_eqn_10_15__p_s(double P, double S_Th, double S_p) {
    std::vector<double> result;
    double p_s = ((P * (S_Th - S_p)) / S_Th);
    result.push_back(p_s);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_16__P(double S_0, double S_Th, double p_0) {
    std::vector<std::complex<double>> result;
    std::complex<double> P = ((p_0 * std::pow((S_Th / S_0), (5.0 / 3.0))) / (std::pow((S_Th / S_0), 1.66666666666667) - 1.0));
    result.push_back(P);
    P = ((p_0 * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) / (std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0) - 1.0));
    result.push_back(P);
    P = ((p_0 * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) / (std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0) - 1.0));
    result.push_back(P);
    return result;
}
std::vector<double> LiquidRing_eqn_10_16__S_0(double P, double S_Th, double p_0) {
    std::vector<double> result;
    double S_0 = (S_Th / std::pow((P / (P - p_0)), (3.0 / 5.0)));
    result.push_back(S_0);
    return result;
}
std::vector<double> LiquidRing_eqn_10_16__S_Th(double P, double S_0, double p_0) {
    std::vector<double> result;
    double S_Th = (S_0 * std::pow((P / (P - p_0)), (3.0 / 5.0)));
    result.push_back(S_Th);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_16__p_0(double P, double S_0, double S_Th) {
    std::vector<std::complex<double>> result;
    std::complex<double> p_0 = (P - (P / std::pow((S_Th / S_0), (5.0 / 3.0))));
    result.push_back(p_0);
    p_0 = (P - (P / std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)));
    result.push_back(p_0);
    p_0 = (P - (P / std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)));
    result.push_back(p_0);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_17__P(double S_0, double S_Th, double p_0, double p_s) {
    std::vector<std::complex<double>> result;
    std::complex<double> P = (((p_0 * std::pow((S_Th / S_0), (5.0 / 3.0))) - p_s) / (std::pow((S_Th / S_0), 1.66666666666667) - 1.0));
    result.push_back(P);
    P = (((p_0 * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) - p_s) / (std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0) - 1.0));
    result.push_back(P);
    P = (((p_0 * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) - p_s) / (std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0) - 1.0));
    result.push_back(P);
    return result;
}
std::vector<double> LiquidRing_eqn_10_17__S_0(double P, double S_Th, double p_0, double p_s) {
    std::vector<double> result;
    double S_0 = (S_Th / std::pow(((P - p_s) / (P - p_0)), (3.0 / 5.0)));
    result.push_back(S_0);
    return result;
}
std::vector<double> LiquidRing_eqn_10_17__S_Th(double P, double S_0, double p_0, double p_s) {
    std::vector<double> result;
    double S_Th = (S_0 * std::pow(((P - p_s) / (P - p_0)), (3.0 / 5.0)));
    result.push_back(S_Th);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_17__p_0(double P, double S_0, double S_Th, double p_s) {
    std::vector<std::complex<double>> result;
    std::complex<double> p_0 = ((((P * std::pow((S_Th / S_0), (5.0 / 3.0))) - P) + p_s) / std::pow((S_Th / S_0), (5.0 / 3.0)));
    result.push_back(p_0);
    p_0 = ((((P * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) - P) + p_s) / std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0));
    result.push_back(p_0);
    p_0 = ((((P * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) - P) + p_s) / std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0));
    result.push_back(p_0);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_17__p_s(double P, double S_0, double S_Th, double p_0) {
    std::vector<std::complex<double>> result;
    std::complex<double> p_s = ((((-P) * std::pow((S_Th / S_0), (5.0 / 3.0))) + P) + (p_0 * std::pow((S_Th / S_0), (5.0 / 3.0))));
    result.push_back(p_s);
    p_s = ((((-P) * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) + P) + (p_0 * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)));
    result.push_back(p_s);
    p_s = ((((-P) * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)) + P) + (p_0 * std::pow((((-0.5) * std::pow((S_Th / S_0), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_Th / S_0), 0.333333333333333))), 5.0)));
    result.push_back(p_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_18__P(double S_Th, double S_p, double T_e, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double P = ((((((S_Th * T_i) * p_s) + ((460.0 * S_Th) * p_s)) - ((S_p * T_e) * p_c)) - ((460.0 * S_p) * p_c)) / ((((S_Th * T_i) + (460.0 * S_Th)) - (S_p * T_e)) - (460.0 * S_p)));
    result.push_back(P);
    return result;
}
std::vector<double> LiquidRing_eqn_10_18__S_Th(double P, double S_p, double T_e, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double S_Th = ((S_p * ((((P * T_e) + (460.0 * P)) - (T_e * p_c)) - (460.0 * p_c))) / ((((P * T_i) + (460.0 * P)) - (T_i * p_s)) - (460.0 * p_s)));
    result.push_back(S_Th);
    return result;
}
std::vector<double> LiquidRing_eqn_10_18__S_p(double P, double S_Th, double T_e, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double S_p = ((S_Th * ((((P * T_i) + (460.0 * P)) - (T_i * p_s)) - (460.0 * p_s))) / ((((P * T_e) + (460.0 * P)) - (T_e * p_c)) - (460.0 * p_c)));
    result.push_back(S_p);
    return result;
}
std::vector<double> LiquidRing_eqn_10_18__T_e(double P, double S_Th, double S_p, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double T_e = ((((((((P * S_Th) * T_i) + ((460.0 * P) * S_Th)) - ((460.0 * P) * S_p)) - ((S_Th * T_i) * p_s)) - ((460.0 * S_Th) * p_s)) + ((460.0 * S_p) * p_c)) / (S_p * (P - p_c)));
    result.push_back(T_e);
    return result;
}
std::vector<double> LiquidRing_eqn_10_18__T_i(double P, double S_Th, double S_p, double T_e, double p_c, double p_s) {
    std::vector<double> result;
    double T_i = (((((((((-460.0) * P) * S_Th) + ((P * S_p) * T_e)) + ((460.0 * P) * S_p)) + ((460.0 * S_Th) * p_s)) - ((S_p * T_e) * p_c)) - ((460.0 * S_p) * p_c)) / (S_Th * (P - p_s)));
    result.push_back(T_i);
    return result;
}
std::vector<double> LiquidRing_eqn_10_18__p_c(double P, double S_Th, double S_p, double T_e, double T_i, double p_s) {
    std::vector<double> result;
    double p_c = (((((((((-P) * S_Th) * T_i) - ((460.0 * P) * S_Th)) + ((P * S_p) * T_e)) + ((460.0 * P) * S_p)) + ((S_Th * T_i) * p_s)) + ((460.0 * S_Th) * p_s)) / (S_p * (T_e + 460.0)));
    result.push_back(p_c);
    return result;
}
std::vector<double> LiquidRing_eqn_10_18__p_s(double P, double S_Th, double S_p, double T_e, double T_i, double p_c) {
    std::vector<double> result;
    double p_s = ((((((((P * S_Th) * T_i) + ((460.0 * P) * S_Th)) - ((P * S_p) * T_e)) - ((460.0 * P) * S_p)) + ((S_p * T_e) * p_c)) + ((460.0 * S_p) * p_c)) / (S_Th * (T_i + 460.0)));
    result.push_back(p_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_19__P(double S_Th, double S_p, double T_e, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double R = std::pow((S_p / S_Th), 1.666667);
    double P = ((((R * (460.0 + T_e)) * p_c) - ((460.0 + T_i) * p_s)) / ((R * (460.0 + T_e)) - (460.0 + T_i)));
    return {P};
}
std::vector<double> LiquidRing_eqn_10_19__S_Th(double P, double S_p, double T_e, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double S_Th = (S_p / std::pow((((((P * T_i) + (460.0 * P)) - (T_i * p_s)) - (460.0 * p_s)) / ((((P * T_e) + (460.0 * P)) - (T_e * p_c)) - (460.0 * p_c))), (3.0 / 5.0)));
    result.push_back(S_Th);
    return result;
}
std::vector<double> LiquidRing_eqn_10_19__S_p(double P, double S_Th, double T_e, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double S_p = (S_Th * std::pow((((((P * T_i) + (460.0 * P)) - (T_i * p_s)) - (460.0 * p_s)) / ((((P * T_e) + (460.0 * P)) - (T_e * p_c)) - (460.0 * p_c))), (3.0 / 5.0)));
    result.push_back(S_p);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_19__T_e(double P, double S_Th, double S_p, double T_i, double p_c, double p_s) {
    std::vector<std::complex<double>> result;
    std::complex<double> T_e = (((((((P * T_i) - ((460.0 * P) * std::pow((S_p / S_Th), (5.0 / 3.0)))) + (460.0 * P)) - (T_i * p_s)) + ((460.0 * p_c) * std::pow((S_p / S_Th), (5.0 / 3.0)))) - (460.0 * p_s)) / (std::pow((S_p / S_Th), (5.0 / 3.0)) * (P - p_c)));
    result.push_back(T_e);
    T_e = (((((((P * T_i) - ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + (460.0 * P)) - (T_i * p_s)) + ((460.0 * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - (460.0 * p_s)) / ((P - p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)));
    result.push_back(T_e);
    T_e = (((((((P * T_i) - ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + (460.0 * P)) - (T_i * p_s)) + ((460.0 * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - (460.0 * p_s)) / ((P - p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)));
    result.push_back(T_e);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_19__T_i(double P, double S_Th, double S_p, double T_e, double p_c, double p_s) {
    std::vector<std::complex<double>> result;
    std::complex<double> T_i = ((((((((P * T_e) * std::pow((S_p / S_Th), (5.0 / 3.0))) + ((460.0 * P) * std::pow((S_p / S_Th), (5.0 / 3.0)))) - (460.0 * P)) - ((T_e * p_c) * std::pow((S_p / S_Th), (5.0 / 3.0)))) - ((460.0 * p_c) * std::pow((S_p / S_Th), (5.0 / 3.0)))) + (460.0 * p_s)) / (P - p_s));
    result.push_back(T_i);
    T_i = ((((((((P * T_e) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)) + ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - (460.0 * P)) - ((T_e * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - ((460.0 * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + (460.0 * p_s)) / (P - p_s));
    result.push_back(T_i);
    T_i = ((((((((P * T_e) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)) + ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - (460.0 * P)) - ((T_e * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - ((460.0 * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + (460.0 * p_s)) / (P - p_s));
    result.push_back(T_i);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_19__p_c(double P, double S_Th, double S_p, double T_e, double T_i, double p_s) {
    std::vector<std::complex<double>> result;
    std::complex<double> p_c = ((((((((P * T_e) * std::pow((S_p / S_Th), (5.0 / 3.0))) - (P * T_i)) + ((460.0 * P) * std::pow((S_p / S_Th), (5.0 / 3.0)))) - (460.0 * P)) + (T_i * p_s)) + (460.0 * p_s)) / (std::pow((S_p / S_Th), (5.0 / 3.0)) * (T_e + 460.0)));
    result.push_back(p_c);
    p_c = ((((((((P * T_e) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)) - (P * T_i)) + ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - (460.0 * P)) + (T_i * p_s)) + (460.0 * p_s)) / ((T_e + 460.0) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)));
    result.push_back(p_c);
    p_c = ((((((((P * T_e) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)) - (P * T_i)) + ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) - (460.0 * P)) + (T_i * p_s)) + (460.0 * p_s)) / ((T_e + 460.0) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)));
    result.push_back(p_c);
    return result;
}
std::vector<std::complex<double>> LiquidRing_eqn_10_19__p_s(double P, double S_Th, double S_p, double T_e, double T_i, double p_c) {
    std::vector<std::complex<double>> result;
    std::complex<double> p_s = (((((((((-P) * T_e) * std::pow((S_p / S_Th), (5.0 / 3.0))) + (P * T_i)) - ((460.0 * P) * std::pow((S_p / S_Th), (5.0 / 3.0)))) + (460.0 * P)) + ((T_e * p_c) * std::pow((S_p / S_Th), (5.0 / 3.0)))) + ((460.0 * p_c) * std::pow((S_p / S_Th), (5.0 / 3.0)))) / (T_i + 460.0));
    result.push_back(p_s);
    p_s = (((((((((-P) * T_e) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)) + (P * T_i)) - ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + (460.0 * P)) + ((T_e * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + ((460.0 * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) - ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) / (T_i + 460.0));
    result.push_back(p_s);
    p_s = (((((((((-P) * T_e) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0)) + (P * T_i)) - ((460.0 * P) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + (460.0 * P)) + ((T_e * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) + ((460.0 * p_c) * std::pow((((-0.5) * std::pow((S_p / S_Th), 0.333333333333333)) + ((0.866025403784439 * std::complex<double>(0.0, 1.0)) * std::pow((S_p / S_Th), 0.333333333333333))), 5.0))) / (T_i + 460.0));
    result.push_back(p_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_2__PS(double Q_gas, double V, double dP, double dt) {
    std::vector<double> result;
    double PS = (Q_gas - ((V * dP) / dt));
    result.push_back(PS);
    return result;
}
std::vector<double> LiquidRing_eqn_10_2__Q_gas(double PS, double V, double dP, double dt) {
    std::vector<double> result;
    double Q_gas = (PS + ((V * dP) / dt));
    result.push_back(Q_gas);
    return result;
}
std::vector<double> LiquidRing_eqn_10_2__V(double PS, double Q_gas, double dP, double dt) {
    std::vector<double> result;
    double V = ((dt * ((-PS) + Q_gas)) / dP);
    result.push_back(V);
    return result;
}
std::vector<double> LiquidRing_eqn_10_2__dP(double PS, double Q_gas, double V, double dt) {
    std::vector<double> result;
    double dP = ((dt * ((-PS) + Q_gas)) / V);
    result.push_back(dP);
    return result;
}
std::vector<double> LiquidRing_eqn_10_2__dt(double PS, double Q_gas, double V, double dP) {
    std::vector<double> result;
    double dt = (((-V) * dP) / (PS - Q_gas));
    result.push_back(dt);
    return result;
}
std::vector<double> LiquidRing_eqn_10_20__P(double S_0, double S_p, double T_e, double T_i, double p_0, double p_c, double p_s) {
    throw std::runtime_error("LiquidRing_eqn_10_20__P: requires numerical solver (not transpilable)");
}
std::vector<double> LiquidRing_eqn_10_20__S_0(double P, double S_p, double T_e, double T_i, double p_0, double p_c, double p_s) {
    std::vector<double> result;
    double S_0 = (S_p * std::pow((((((((((std::pow(P, 2.0) * T_i) + (460.0 * std::pow(P, 2.0))) - ((P * T_i) * p_0)) - ((P * T_i) * p_c)) - ((460.0 * P) * p_0)) - ((460.0 * P) * p_c)) + ((T_i * p_0) * p_c)) + ((460.0 * p_0) * p_c)) / (P * ((((P * T_e) + (460.0 * P)) - (T_e * p_s)) - (460.0 * p_s)))), (3.0 / 5.0)));
    result.push_back(S_0);
    return result;
}
std::vector<double> LiquidRing_eqn_10_20__S_p(double P, double S_0, double T_e, double T_i, double p_0, double p_c, double p_s) {
    std::vector<double> result;
    double S_p = (S_0 / std::pow((((((((((std::pow(P, 2.0) * T_i) + (460.0 * std::pow(P, 2.0))) - ((P * T_i) * p_0)) - ((P * T_i) * p_c)) - ((460.0 * P) * p_0)) - ((460.0 * P) * p_c)) + ((T_i * p_0) * p_c)) + ((460.0 * p_0) * p_c)) / (P * ((((P * T_e) + (460.0 * P)) - (T_e * p_s)) - (460.0 * p_s)))), (3.0 / 5.0)));
    result.push_back(S_p);
    return result;
}
std::vector<double> LiquidRing_eqn_10_20__T_e(double P, double S_0, double S_p, double T_i, double p_0, double p_c, double p_s) {
    std::vector<double> result;
    double R = std::pow((S_0 / S_p), 1.666667);
    double T_e = (((((P - p_0) * (460.0 + T_i)) * (P - p_c)) / (R * (P * (P - p_s)))) - 460.0);
    return {T_e};
}
std::vector<double> LiquidRing_eqn_10_20__T_i(double P, double S_0, double S_p, double T_e, double p_0, double p_c, double p_s) {
    std::vector<double> result;
    double R = std::pow((S_0 / S_p), 1.666667);
    double T_i = (((R * ((P * (P - p_s)) * (460.0 + T_e))) / ((P - p_0) * (P - p_c))) - 460.0);
    return {T_i};
}
std::vector<double> LiquidRing_eqn_10_20__p_0(double P, double S_0, double S_p, double T_e, double T_i, double p_c, double p_s) {
    std::vector<double> result;
    double R = std::pow((S_0 / S_p), 1.666667);
    double p_0 = (P - ((R * ((P * (P - p_s)) * (460.0 + T_e))) / ((460.0 + T_i) * (P - p_c))));
    return {p_0};
}
std::vector<double> LiquidRing_eqn_10_20__p_c(double P, double S_0, double S_p, double T_e, double T_i, double p_0, double p_s) {
    std::vector<double> result;
    double R = std::pow((S_0 / S_p), 1.666667);
    double p_c = (P - ((R * ((P * (P - p_s)) * (460.0 + T_e))) / ((P - p_0) * (460.0 + T_i))));
    return {p_c};
}
std::vector<double> LiquidRing_eqn_10_20__p_s(double P, double S_0, double S_p, double T_e, double T_i, double p_0, double p_c) {
    std::vector<double> result;
    double R = std::pow((S_0 / S_p), 1.666667);
    double p_s = (P - ((((P - p_0) * (460.0 + T_i)) * (P - p_c)) / (R * (P * (460.0 + T_e)))));
    return {p_s};
}
std::vector<double> LiquidRing_eqn_10_21__P(double P_d, double P_prime) {
    std::vector<double> result;
    double P = ((P_d * P_prime) / 760.0);
    result.push_back(P);
    return result;
}
std::vector<double> LiquidRing_eqn_10_21__P_d(double P, double P_prime) {
    std::vector<double> result;
    double P_d = ((760.0 * P) / P_prime);
    result.push_back(P_d);
    return result;
}
std::vector<double> LiquidRing_eqn_10_21__P_prime(double P, double P_d) {
    std::vector<double> result;
    double P_prime = ((760.0 * P) / P_d);
    result.push_back(P_prime);
    return result;
}
std::vector<double> LiquidRing_eqn_10_3__N_mfw(double Q_gas, double T) {
    std::vector<double> result;
    double N_mfw = ((0.108108108108108 * Q_gas) / T);
    result.push_back(N_mfw);
    return result;
}
std::vector<double> LiquidRing_eqn_10_3__Q_gas(double N_mfw, double T) {
    std::vector<double> result;
    double Q_gas = ((9.25 * N_mfw) * T);
    result.push_back(Q_gas);
    return result;
}
std::vector<double> LiquidRing_eqn_10_3__T(double N_mfw, double Q_gas) {
    std::vector<double> result;
    double T = ((0.108108108108108 * Q_gas) / N_mfw);
    result.push_back(T);
    return result;
}
std::vector<double> LiquidRing_eqn_10_4__Q_gas(double SP_1, double SP_2, double S_p, double V, double t) {
    std::vector<double> result;
    double Q_gas = ((-(SP_1 - (SP_2 * std::exp(((S_p * t) / V))))) / (std::exp(((S_p * t) / V)) - 1.0));
    result.push_back(Q_gas);
    return result;
}
std::vector<double> LiquidRing_eqn_10_4__SP_1(double Q_gas, double SP_2, double S_p, double V, double t) {
    std::vector<double> result;
    double SP_1 = (Q_gas + (((-Q_gas) + SP_2) * std::exp(((S_p * t) / V))));
    result.push_back(SP_1);
    return result;
}
std::vector<double> LiquidRing_eqn_10_4__SP_2(double Q_gas, double SP_1, double S_p, double V, double t) {
    std::vector<double> result;
    double SP_2 = ((((Q_gas * std::exp(((S_p * t) / V))) - Q_gas) + SP_1) * std::exp((((-S_p) * t) / V)));
    result.push_back(SP_2);
    return result;
}
std::vector<double> LiquidRing_eqn_10_4__S_p(double Q_gas, double SP_1, double SP_2, double V, double t) {
    std::vector<double> result;
    double S_p = ((V * std::log(((Q_gas - SP_1) / (Q_gas - SP_2)))) / t);
    result.push_back(S_p);
    return result;
}
std::vector<double> LiquidRing_eqn_10_4__V(double Q_gas, double SP_1, double SP_2, double S_p, double t) {
    std::vector<double> result;
    double V = ((S_p * t) / std::log(((Q_gas - SP_1) / (Q_gas - SP_2))));
    result.push_back(V);
    return result;
}
std::vector<double> LiquidRing_eqn_10_4__t(double Q_gas, double SP_1, double SP_2, double S_p, double V) {
    std::vector<double> result;
    double t = ((V * std::log(((Q_gas - SP_1) / (Q_gas - SP_2)))) / S_p);
    result.push_back(t);
    return result;
}
std::vector<double> LiquidRing_eqn_10_5__P_1(double P_2, double S_p, double V, double t) {
    std::vector<double> result;
    double P_1 = (P_2 * std::exp(((S_p * t) / V)));
    result.push_back(P_1);
    return result;
}
std::vector<double> LiquidRing_eqn_10_5__P_2(double P_1, double S_p, double V, double t) {
    std::vector<double> result;
    double P_2 = (P_1 * std::exp((((-S_p) * t) / V)));
    result.push_back(P_2);
    return result;
}
std::vector<double> LiquidRing_eqn_10_5__S_p(double P_1, double P_2, double V, double t) {
    std::vector<double> result;
    double S_p = ((V * std::log((P_1 / P_2))) / t);
    result.push_back(S_p);
    return result;
}
std::vector<double> LiquidRing_eqn_10_5__V(double P_1, double P_2, double S_p, double t) {
    std::vector<double> result;
    double V = ((S_p * t) / std::log((P_1 / P_2)));
    result.push_back(V);
    return result;
}
std::vector<double> LiquidRing_eqn_10_5__t(double P_1, double P_2, double S_p, double V) {
    std::vector<double> result;
    double t = ((V * std::log((P_1 / P_2))) / S_p);
    result.push_back(t);
    return result;
}
std::vector<double> LiquidRing_eqn_10_6__P_1(double P_2, double S_a, double V, double t) {
    std::vector<double> result;
    double P_1 = (P_2 * std::exp(((S_a * t) / V)));
    result.push_back(P_1);
    return result;
}
std::vector<double> LiquidRing_eqn_10_6__P_2(double P_1, double S_a, double V, double t) {
    std::vector<double> result;
    double P_2 = (P_1 * std::exp((((-S_a) * t) / V)));
    result.push_back(P_2);
    return result;
}
std::vector<double> LiquidRing_eqn_10_6__S_a(double P_1, double P_2, double V, double t) {
    std::vector<double> result;
    double S_a = ((V * std::log((P_1 / P_2))) / t);
    result.push_back(S_a);
    return result;
}
std::vector<double> LiquidRing_eqn_10_6__V(double P_1, double P_2, double S_a, double t) {
    std::vector<double> result;
    double V = ((S_a * t) / std::log((P_1 / P_2)));
    result.push_back(V);
    return result;
}
std::vector<double> LiquidRing_eqn_10_6__t(double P_1, double P_2, double S_a, double V) {
    std::vector<double> result;
    double t = ((V * std::log((P_1 / P_2))) / S_a);
    result.push_back(t);
    return result;
}
std::vector<double> LiquidRing_eqn_10_8__bhp(double c_p, double delta_T, double delta_h_i, double f_a, double rho, double w_i) {
    std::vector<double> result;
    double bhp = (((((0.00315127701375246 * c_p) * delta_T) * f_a) * rho) - ((0.000392927308447937 * delta_h_i) * w_i));
    result.push_back(bhp);
    return result;
}
std::vector<double> LiquidRing_eqn_10_8__c_p(double bhp, double delta_T, double delta_h_i, double f_a, double rho, double w_i) {
    std::vector<double> result;
    double c_p = ((0.124688279301746 * ((2545.0 * bhp) + (delta_h_i * w_i))) / ((delta_T * f_a) * rho));
    result.push_back(c_p);
    return result;
}
std::vector<double> LiquidRing_eqn_10_8__delta_T(double bhp, double c_p, double delta_h_i, double f_a, double rho, double w_i) {
    std::vector<double> result;
    double delta_T = ((0.124688279301746 * ((2545.0 * bhp) + (delta_h_i * w_i))) / ((c_p * f_a) * rho));
    result.push_back(delta_T);
    return result;
}
std::vector<double> LiquidRing_eqn_10_8__delta_h_i(double bhp, double c_p, double delta_T, double f_a, double rho, double w_i) {
    std::vector<double> result;
    double delta_h_i = ((0.02 * (((-127250.0) * bhp) + ((((401.0 * c_p) * delta_T) * f_a) * rho))) / w_i);
    result.push_back(delta_h_i);
    return result;
}
std::vector<double> LiquidRing_eqn_10_8__f_a(double bhp, double c_p, double delta_T, double delta_h_i, double rho, double w_i) {
    std::vector<double> result;
    double f_a = ((0.124688279301746 * ((2545.0 * bhp) + (delta_h_i * w_i))) / ((c_p * delta_T) * rho));
    result.push_back(f_a);
    return result;
}
std::vector<double> LiquidRing_eqn_10_8__rho(double bhp, double c_p, double delta_T, double delta_h_i, double f_a, double w_i) {
    std::vector<double> result;
    double rho = ((0.124688279301746 * ((2545.0 * bhp) + (delta_h_i * w_i))) / ((c_p * delta_T) * f_a));
    result.push_back(rho);
    return result;
}
std::vector<double> LiquidRing_eqn_10_8__w_i(double bhp, double c_p, double delta_T, double delta_h_i, double f_a, double rho) {
    std::vector<double> result;
    double w_i = ((0.02 * (((-127250.0) * bhp) + ((((401.0 * c_p) * delta_T) * f_a) * rho))) / delta_h_i);
    result.push_back(w_i);
    return result;
}
std::vector<double> LiquidRing_eqn_10_9__T_c(double T_s, double delta_T) {
    std::vector<double> result;
    double T_c = (T_s + delta_T);
    result.push_back(T_c);
    return result;
}
std::vector<double> LiquidRing_eqn_10_9__T_s(double T_c, double delta_T) {
    std::vector<double> result;
    double T_s = (T_c - delta_T);
    result.push_back(T_s);
    return result;
}
std::vector<double> LiquidRing_eqn_10_9__delta_T(double T_c, double T_s) {
    std::vector<double> result;
    double delta_T = (T_c - T_s);
    result.push_back(delta_T);
    return result;
}
std::vector<double> Precondensors_eqn_7_1__P(double p_i, double y_i) {
    std::vector<double> result;
    double P = (p_i / y_i);
    result.push_back(P);
    return result;
}
std::vector<double> Precondensors_eqn_7_1__p_i(double P, double y_i) {
    std::vector<double> result;
    double p_i = (P * y_i);
    result.push_back(p_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_1__y_i(double P, double p_i) {
    std::vector<double> result;
    double y_i = (p_i / P);
    result.push_back(y_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_10__L_c_P(double Q_condensor_heat_duty, double del_T) {
    std::vector<double> result;
    double L_c_P = (Q_condensor_heat_duty / (500.0 * del_T));
    result.push_back(L_c_P);
    return result;
}
std::vector<double> Precondensors_eqn_7_10__Q_condensor_heat_duty(double L_c_P, double del_T) {
    std::vector<double> result;
    double Q_condensor_heat_duty = ((500.0 * L_c_P) * del_T);
    result.push_back(Q_condensor_heat_duty);
    return result;
}
std::vector<double> Precondensors_eqn_7_10__del_T(double L_c_P, double Q_condensor_heat_duty) {
    std::vector<double> result;
    double del_T = (Q_condensor_heat_duty / (500.0 * L_c_P));
    result.push_back(del_T);
    return result;
}
std::vector<double> Precondensors_eqn_7_11__Q_condensor_heat_duty(double U_v, double V_c, double del_T_LM) {
    std::vector<double> result;
    double Q_condensor_heat_duty = ((U_v * V_c) * del_T_LM);
    result.push_back(Q_condensor_heat_duty);
    return result;
}
std::vector<double> Precondensors_eqn_7_11__U_v(double Q_condensor_heat_duty, double V_c, double del_T_LM) {
    std::vector<double> result;
    double U_v = (Q_condensor_heat_duty / (V_c * del_T_LM));
    result.push_back(U_v);
    return result;
}
std::vector<double> Precondensors_eqn_7_11__V_c(double Q_condensor_heat_duty, double U_v, double del_T_LM) {
    std::vector<double> result;
    double V_c = (Q_condensor_heat_duty / (U_v * del_T_LM));
    result.push_back(V_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_11__del_T_LM(double Q_condensor_heat_duty, double U_v, double V_c) {
    std::vector<double> result;
    double del_T_LM = (Q_condensor_heat_duty / (U_v * V_c));
    result.push_back(del_T_LM);
    return result;
}
std::vector<double> Precondensors_eqn_7_12__A(double Q_condensor_heat_duty, double U, double del_T) {
    std::vector<double> result;
    double A = (Q_condensor_heat_duty / (U * del_T));
    result.push_back(A);
    return result;
}
std::vector<double> Precondensors_eqn_7_12__Q_condensor_heat_duty(double A, double U, double del_T) {
    std::vector<double> result;
    double Q_condensor_heat_duty = ((A * U) * del_T);
    result.push_back(Q_condensor_heat_duty);
    return result;
}
std::vector<double> Precondensors_eqn_7_12__U(double A, double Q_condensor_heat_duty, double del_T) {
    std::vector<double> result;
    double U = (Q_condensor_heat_duty / (A * del_T));
    result.push_back(U);
    return result;
}
std::vector<double> Precondensors_eqn_7_12__del_T(double A, double Q_condensor_heat_duty, double U) {
    std::vector<double> result;
    double del_T = (Q_condensor_heat_duty / (A * U));
    result.push_back(del_T);
    return result;
}
std::vector<double> Precondensors_eqn_7_14a__A(double Q_condensor_heat_duty, double U, double del_T_LM) {
    std::vector<double> result;
    double A = (Q_condensor_heat_duty / (U * del_T_LM));
    result.push_back(A);
    return result;
}
std::vector<double> Precondensors_eqn_7_14a__Q_condensor_heat_duty(double A, double U, double del_T_LM) {
    std::vector<double> result;
    double Q_condensor_heat_duty = ((A * U) * del_T_LM);
    result.push_back(Q_condensor_heat_duty);
    return result;
}
std::vector<double> Precondensors_eqn_7_14a__U(double A, double Q_condensor_heat_duty, double del_T_LM) {
    std::vector<double> result;
    double U = (Q_condensor_heat_duty / (A * del_T_LM));
    result.push_back(U);
    return result;
}
std::vector<double> Precondensors_eqn_7_14a__del_T_LM(double A, double Q_condensor_heat_duty, double U) {
    std::vector<double> result;
    double del_T_LM = (Q_condensor_heat_duty / (A * U));
    result.push_back(del_T_LM);
    return result;
}
std::vector<double> Precondensors_eqn_7_14b__A(double Q_condensor_heat_duty, double U, double del_T_1, double del_T_2) {
    std::vector<double> result;
    double A = (Q_condensor_heat_duty / ((U * (del_T_1 - del_T_2)) * std::log((del_T_1 - del_T_2))));
    result.push_back(A);
    return result;
}
std::vector<double> Precondensors_eqn_7_14b__Q_condensor_heat_duty(double A, double U, double del_T_1, double del_T_2) {
    std::vector<double> result;
    double Q_condensor_heat_duty = (((A * U) * (del_T_1 - del_T_2)) * std::log((del_T_1 - del_T_2)));
    result.push_back(Q_condensor_heat_duty);
    return result;
}
std::vector<double> Precondensors_eqn_7_14b__U(double A, double Q_condensor_heat_duty, double del_T_1, double del_T_2) {
    std::vector<double> result;
    double U = (Q_condensor_heat_duty / ((A * (del_T_1 - del_T_2)) * std::log((del_T_1 - del_T_2))));
    result.push_back(U);
    return result;
}
std::vector<double> Precondensors_eqn_7_14b__del_T_1(double A, double Q_condensor_heat_duty, double U, double del_T_2) {
    std::vector<double> result;
    double del_T_1 = (del_T_2 + std::exp(_lambertW((Q_condensor_heat_duty / (A * U)))));
    result.push_back(del_T_1);
    return result;
}
std::vector<double> Precondensors_eqn_7_14b__del_T_2(double A, double Q_condensor_heat_duty, double U, double del_T_1) {
    std::vector<double> result;
    double del_T_2 = (del_T_1 - std::exp(_lambertW((Q_condensor_heat_duty / (A * U)))));
    result.push_back(del_T_2);
    return result;
}
std::vector<double> Precondensors_eqn_7_15__U(double sum_R) {
    std::vector<double> result;
    double U = (1.0 / sum_R);
    result.push_back(U);
    return result;
}
std::vector<double> Precondensors_eqn_7_15__sum_R(double U) {
    std::vector<double> result;
    double sum_R = (1.0 / U);
    result.push_back(sum_R);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__D_0(double D_LM, double D_i, double R_f_0, double R_fi, double U_0, double h_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double D_0 = (((((D_LM * D_i) * h_i) * k_w) * (((((-R_f_0) * U_0) * h_0) - U_0) + h_0)) / ((U_0 * h_0) * (((((D_LM * R_fi) * h_i) * k_w) + (D_LM * k_w)) + ((D_i * h_i) * x_w))));
    result.push_back(D_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__D_LM(double D_0, double D_i, double R_f_0, double R_fi, double U_0, double h_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double D_LM = (((((((-D_0) * D_i) * U_0) * h_0) * h_i) * x_w) / (k_w * ((((((((D_0 * R_fi) * U_0) * h_0) * h_i) + ((D_0 * U_0) * h_0)) + ((((D_i * R_f_0) * U_0) * h_0) * h_i)) + ((D_i * U_0) * h_i)) - ((D_i * h_0) * h_i))));
    result.push_back(D_LM);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__D_i(double D_0, double D_LM, double R_f_0, double R_fi, double U_0, double h_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double D_i = (((((((-D_0) * D_LM) * U_0) * h_0) * k_w) * ((R_fi * h_i) + 1.0)) / (h_i * ((((((D_0 * U_0) * h_0) * x_w) + ((((D_LM * R_f_0) * U_0) * h_0) * k_w)) + ((D_LM * U_0) * k_w)) - ((D_LM * h_0) * k_w))));
    result.push_back(D_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__R_f_0(double D_0, double D_LM, double D_i, double R_fi, double U_0, double h_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double R_f_0 = (((((((-D_0) * R_fi) / D_i) - (D_0 / (D_i * h_i))) - ((D_0 * x_w) / (D_LM * k_w))) - (1.0 / h_0)) + (1.0 / U_0));
    result.push_back(R_f_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__R_fi(double D_0, double D_LM, double D_i, double R_f_0, double U_0, double h_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double R_fi = ((((((-1.0) / h_i) - ((D_i * x_w) / (D_LM * k_w))) - ((D_i * R_f_0) / D_0)) - (D_i / (D_0 * h_0))) + (D_i / (D_0 * U_0)));
    result.push_back(R_fi);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__U_0(double D_0, double D_LM, double D_i, double R_f_0, double R_fi, double h_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double U_0 = (((((D_LM * D_i) * h_0) * h_i) * k_w) / (((((((((D_0 * D_LM) * R_fi) * h_0) * h_i) * k_w) + (((D_0 * D_LM) * h_0) * k_w)) + ((((D_0 * D_i) * h_0) * h_i) * x_w)) + (((((D_LM * D_i) * R_f_0) * h_0) * h_i) * k_w)) + (((D_LM * D_i) * h_i) * k_w)));
    result.push_back(U_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__h_0(double D_0, double D_LM, double D_i, double R_f_0, double R_fi, double U_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double h_0 = ((((((-D_LM) * D_i) * U_0) * h_i) * k_w) / (((((((((D_0 * D_LM) * R_fi) * U_0) * h_i) * k_w) + (((D_0 * D_LM) * U_0) * k_w)) + ((((D_0 * D_i) * U_0) * h_i) * x_w)) + (((((D_LM * D_i) * R_f_0) * U_0) * h_i) * k_w)) - (((D_LM * D_i) * h_i) * k_w)));
    result.push_back(h_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__h_i(double D_0, double D_LM, double D_i, double R_f_0, double R_fi, double U_0, double h_0, double k_w, double x_w) {
    std::vector<double> result;
    double h_i = ((((((-D_0) * D_LM) * U_0) * h_0) * k_w) / (((((((((D_0 * D_LM) * R_fi) * U_0) * h_0) * k_w) + ((((D_0 * D_i) * U_0) * h_0) * x_w)) + (((((D_LM * D_i) * R_f_0) * U_0) * h_0) * k_w)) + (((D_LM * D_i) * U_0) * k_w)) - (((D_LM * D_i) * h_0) * k_w)));
    result.push_back(h_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__k_w(double D_0, double D_LM, double D_i, double R_f_0, double R_fi, double U_0, double h_0, double h_i, double x_w) {
    std::vector<double> result;
    double k_w = (((((((-D_0) * D_i) * U_0) * h_0) * h_i) * x_w) / (D_LM * ((((((((D_0 * R_fi) * U_0) * h_0) * h_i) + ((D_0 * U_0) * h_0)) + ((((D_i * R_f_0) * U_0) * h_0) * h_i)) + ((D_i * U_0) * h_i)) - ((D_i * h_0) * h_i))));
    result.push_back(k_w);
    return result;
}
std::vector<double> Precondensors_eqn_7_16__x_w(double D_0, double D_LM, double D_i, double R_f_0, double R_fi, double U_0, double h_0, double h_i, double k_w) {
    std::vector<double> result;
    double x_w = ((((((((-D_LM) * R_fi) * k_w) / D_i) - ((D_LM * k_w) / (D_i * h_i))) - (((D_LM * R_f_0) * k_w) / D_0)) - ((D_LM * k_w) / (D_0 * h_0))) + ((D_LM * k_w) / (D_0 * U_0)));
    result.push_back(x_w);
    return result;
}
std::vector<double> Precondensors_eqn_7_17__R_0(double R_nc, double h_c) {
    std::vector<double> result;
    double R_0 = (R_nc + (1.0 / h_c));
    result.push_back(R_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_17__R_nc(double R_0, double h_c) {
    std::vector<double> result;
    double R_nc = (R_0 - (1.0 / h_c));
    result.push_back(R_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_17__h_c(double R_0, double R_nc) {
    std::vector<double> result;
    double h_c = (1.0 / (R_0 - R_nc));
    result.push_back(h_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__D_0(double D_LM, double D_i, double R_fi, double R_fo, double R_nc, double U_0, double h_c, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double D_0 = (((((D_LM * D_i) * h_i) * k_w) * ((((((-R_fo) * U_0) * h_c) - ((R_nc * U_0) * h_c)) - U_0) + h_c)) / ((U_0 * h_c) * (((((D_LM * R_fi) * h_i) * k_w) + (D_LM * k_w)) + ((D_i * h_i) * x_w))));
    result.push_back(D_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__D_LM(double D_0, double D_i, double R_fi, double R_fo, double R_nc, double U_0, double h_c, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double D_LM = (((((((-D_0) * D_i) * U_0) * h_c) * h_i) * x_w) / (k_w * (((((((((D_0 * R_fi) * U_0) * h_c) * h_i) + ((D_0 * U_0) * h_c)) + ((((D_i * R_fo) * U_0) * h_c) * h_i)) + ((((D_i * R_nc) * U_0) * h_c) * h_i)) + ((D_i * U_0) * h_i)) - ((D_i * h_c) * h_i))));
    result.push_back(D_LM);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__D_i(double D_0, double D_LM, double R_fi, double R_fo, double R_nc, double U_0, double h_c, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double D_i = (((((((-D_0) * D_LM) * U_0) * h_c) * k_w) * ((R_fi * h_i) + 1.0)) / (h_i * (((((((D_0 * U_0) * h_c) * x_w) + ((((D_LM * R_fo) * U_0) * h_c) * k_w)) + ((((D_LM * R_nc) * U_0) * h_c) * k_w)) + ((D_LM * U_0) * k_w)) - ((D_LM * h_c) * k_w))));
    result.push_back(D_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__R_fi(double D_0, double D_LM, double D_i, double R_fo, double R_nc, double U_0, double h_c, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double R_fi = (((((((-1.0) / h_i) - ((D_i * x_w) / (D_LM * k_w))) - ((D_i * R_fo) / D_0)) - ((D_i * R_nc) / D_0)) - (D_i / (D_0 * h_c))) + (D_i / (D_0 * U_0)));
    result.push_back(R_fi);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__R_fo(double D_0, double D_LM, double D_i, double R_fi, double R_nc, double U_0, double h_c, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double R_fo = ((((((((-D_0) * R_fi) / D_i) - (D_0 / (D_i * h_i))) - ((D_0 * x_w) / (D_LM * k_w))) - R_nc) - (1.0 / h_c)) + (1.0 / U_0));
    result.push_back(R_fo);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__R_nc(double D_0, double D_LM, double D_i, double R_fi, double R_fo, double U_0, double h_c, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double R_nc = ((((((((-D_0) * R_fi) / D_i) - (D_0 / (D_i * h_i))) - ((D_0 * x_w) / (D_LM * k_w))) - R_fo) - (1.0 / h_c)) + (1.0 / U_0));
    result.push_back(R_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__U_0(double D_0, double D_LM, double D_i, double R_fi, double R_fo, double R_nc, double h_c, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double U_0 = (((((D_LM * D_i) * h_c) * h_i) * k_w) / ((((((((((D_0 * D_LM) * R_fi) * h_c) * h_i) * k_w) + (((D_0 * D_LM) * h_c) * k_w)) + ((((D_0 * D_i) * h_c) * h_i) * x_w)) + (((((D_LM * D_i) * R_fo) * h_c) * h_i) * k_w)) + (((((D_LM * D_i) * R_nc) * h_c) * h_i) * k_w)) + (((D_LM * D_i) * h_i) * k_w)));
    result.push_back(U_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__h_c(double D_0, double D_LM, double D_i, double R_fi, double R_fo, double R_nc, double U_0, double h_i, double k_w, double x_w) {
    std::vector<double> result;
    double h_c = ((((((-D_LM) * D_i) * U_0) * h_i) * k_w) / ((((((((((D_0 * D_LM) * R_fi) * U_0) * h_i) * k_w) + (((D_0 * D_LM) * U_0) * k_w)) + ((((D_0 * D_i) * U_0) * h_i) * x_w)) + (((((D_LM * D_i) * R_fo) * U_0) * h_i) * k_w)) + (((((D_LM * D_i) * R_nc) * U_0) * h_i) * k_w)) - (((D_LM * D_i) * h_i) * k_w)));
    result.push_back(h_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__h_i(double D_0, double D_LM, double D_i, double R_fi, double R_fo, double R_nc, double U_0, double h_c, double k_w, double x_w) {
    std::vector<double> result;
    double h_i = ((((((-D_0) * D_LM) * U_0) * h_c) * k_w) / ((((((((((D_0 * D_LM) * R_fi) * U_0) * h_c) * k_w) + ((((D_0 * D_i) * U_0) * h_c) * x_w)) + (((((D_LM * D_i) * R_fo) * U_0) * h_c) * k_w)) + (((((D_LM * D_i) * R_nc) * U_0) * h_c) * k_w)) + (((D_LM * D_i) * U_0) * k_w)) - (((D_LM * D_i) * h_c) * k_w)));
    result.push_back(h_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__k_w(double D_0, double D_LM, double D_i, double R_fi, double R_fo, double R_nc, double U_0, double h_c, double h_i, double x_w) {
    std::vector<double> result;
    double k_w = (((((((-D_0) * D_i) * U_0) * h_c) * h_i) * x_w) / (D_LM * (((((((((D_0 * R_fi) * U_0) * h_c) * h_i) + ((D_0 * U_0) * h_c)) + ((((D_i * R_fo) * U_0) * h_c) * h_i)) + ((((D_i * R_nc) * U_0) * h_c) * h_i)) + ((D_i * U_0) * h_i)) - ((D_i * h_c) * h_i))));
    result.push_back(k_w);
    return result;
}
std::vector<double> Precondensors_eqn_7_18__x_w(double D_0, double D_LM, double D_i, double R_fi, double R_fo, double R_nc, double U_0, double h_c, double h_i, double k_w) {
    std::vector<double> result;
    double x_w = (((((((((-D_LM) * R_fi) * k_w) / D_i) - ((D_LM * k_w) / (D_i * h_i))) - (((D_LM * R_fo) * k_w) / D_0)) - (((D_LM * R_nc) * k_w) / D_0)) - ((D_LM * k_w) / (D_0 * h_c))) + ((D_LM * k_w) / (D_0 * U_0)));
    result.push_back(x_w);
    return result;
}
std::vector<double> Precondensors_eqn_7_2__P_i_0(double p_i, double x_i) {
    std::vector<double> result;
    double P_i_0 = (p_i / x_i);
    result.push_back(P_i_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_2__p_i(double P_i_0, double x_i) {
    std::vector<double> result;
    double p_i = (P_i_0 * x_i);
    result.push_back(p_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_2__x_i(double P_i_0, double p_i) {
    std::vector<double> result;
    double x_i = (p_i / P_i_0);
    result.push_back(x_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_3__P_i_0(double epsilon_i, double p_i, double x_i) {
    std::vector<double> result;
    double P_i_0 = (p_i / (epsilon_i * x_i));
    result.push_back(P_i_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_3__epsilon_i(double P_i_0, double p_i, double x_i) {
    std::vector<double> result;
    double epsilon_i = (p_i / (P_i_0 * x_i));
    result.push_back(epsilon_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_3__p_i(double P_i_0, double epsilon_i, double x_i) {
    std::vector<double> result;
    double p_i = ((P_i_0 * epsilon_i) * x_i);
    result.push_back(p_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_3__x_i(double P_i_0, double epsilon_i, double p_i) {
    std::vector<double> result;
    double x_i = (p_i / (P_i_0 * epsilon_i));
    result.push_back(x_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_4a__P(double p_c, double p_nc) {
    std::vector<double> result;
    double P = (p_c + p_nc);
    result.push_back(P);
    return result;
}
std::vector<double> Precondensors_eqn_7_4a__p_c(double P, double p_nc) {
    std::vector<double> result;
    double p_c = (P - p_nc);
    result.push_back(p_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_4a__p_nc(double P, double p_c) {
    std::vector<double> result;
    double p_nc = (P - p_c);
    result.push_back(p_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_4aa__n_i(double n_nc, double p_i, double p_nc) {
    std::vector<double> result;
    double n_i = ((n_nc * p_i) / p_nc);
    result.push_back(n_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_4aa__n_nc(double n_i, double p_i, double p_nc) {
    std::vector<double> result;
    double n_nc = ((n_i * p_nc) / p_i);
    result.push_back(n_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_4aa__p_i(double n_i, double n_nc, double p_nc) {
    std::vector<double> result;
    double p_i = ((n_i * p_nc) / n_nc);
    result.push_back(p_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_4aa__p_nc(double n_i, double n_nc, double p_i) {
    std::vector<double> result;
    double p_nc = ((n_nc * p_i) / n_i);
    result.push_back(p_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ab__P_c(double p, double p_i, double p_nc) {
    std::vector<double> result;
    double P_c = (p - p_nc);
    result.push_back(P_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ab__p(double P_c, double p_i, double p_nc) {
    std::vector<double> result;
    double p = (P_c + p_nc);
    result.push_back(p);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ab__p_i(double P_c, double p, double p_nc) {
    std::vector<double> result;
    double p_i = 0.0;
    result.push_back(p_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ab__p_nc(double P_c, double p, double p_i) {
    std::vector<double> result;
    double p_nc = ((-P_c) + p);
    result.push_back(p_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ac__P_c(double n_i, double n_nc, double p, double p_i) {
    std::vector<double> result;
    double P_c = (p - ((n_nc * p_i) / n_i));
    result.push_back(P_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ac__n_i(double P_c, double n_nc, double p, double p_i) {
    std::vector<double> result;
    double n_i = ((n_nc * p_i) / ((-P_c) + p));
    result.push_back(n_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ac__n_nc(double P_c, double n_i, double p, double p_i) {
    std::vector<double> result;
    double n_nc = ((n_i * ((-P_c) + p)) / p_i);
    result.push_back(n_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ac__p(double P_c, double n_i, double n_nc, double p_i) {
    std::vector<double> result;
    double p = (P_c + ((n_nc * p_i) / n_i));
    result.push_back(p);
    return result;
}
std::vector<double> Precondensors_eqn_7_4ac__p_i(double P_c, double n_i, double n_nc, double p) {
    std::vector<double> result;
    double p_i = ((n_i * ((-P_c) + p)) / n_nc);
    result.push_back(p_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_5__N_i(double N_nc, double P, double P_c, double p_i) {
    std::vector<double> result;
    double N_i = ((N_nc * p_i) / (P - P_c));
    result.push_back(N_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_5__N_nc(double N_i, double P, double P_c, double p_i) {
    std::vector<double> result;
    double N_nc = ((N_i * (P - P_c)) / p_i);
    result.push_back(N_nc);
    return result;
}
std::vector<double> Precondensors_eqn_7_5__P(double N_i, double N_nc, double P_c, double p_i) {
    std::vector<double> result;
    double P = (P_c + ((N_nc * p_i) / N_i));
    result.push_back(P);
    return result;
}
std::vector<double> Precondensors_eqn_7_5__P_c(double N_i, double N_nc, double P, double p_i) {
    std::vector<double> result;
    double P_c = (P - ((N_nc * p_i) / N_i));
    result.push_back(P_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_5__p_i(double N_i, double N_nc, double P, double P_c) {
    std::vector<double> result;
    double p_i = ((N_i * (P - P_c)) / N_nc);
    result.push_back(p_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_6__M(double P, double P_i_0, double W_air, double W_i, double p_c, double x_i) {
    std::vector<double> result;
    double M = (((29.0 * W_i) * (P - p_c)) / ((P_i_0 * W_air) * x_i));
    result.push_back(M);
    return result;
}
std::vector<double> Precondensors_eqn_7_6__P(double M, double P_i_0, double W_air, double W_i, double p_c, double x_i) {
    std::vector<double> result;
    double P = (((((M * P_i_0) * W_air) * x_i) / (29.0 * W_i)) + p_c);
    result.push_back(P);
    return result;
}
std::vector<double> Precondensors_eqn_7_6__P_i_0(double M, double P, double W_air, double W_i, double p_c, double x_i) {
    std::vector<double> result;
    double P_i_0 = (((29.0 * W_i) * (P - p_c)) / ((M * W_air) * x_i));
    result.push_back(P_i_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_6__W_air(double M, double P, double P_i_0, double W_i, double p_c, double x_i) {
    std::vector<double> result;
    double W_air = (((29.0 * W_i) * (P - p_c)) / ((M * P_i_0) * x_i));
    result.push_back(W_air);
    return result;
}
std::vector<double> Precondensors_eqn_7_6__W_i(double M, double P, double P_i_0, double W_air, double p_c, double x_i) {
    std::vector<double> result;
    double W_i = ((((M * P_i_0) * W_air) * x_i) / (29.0 * (P - p_c)));
    result.push_back(W_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_6__p_c(double M, double P, double P_i_0, double W_air, double W_i, double x_i) {
    std::vector<double> result;
    double p_c = ((((((-M) * P_i_0) * W_air) * x_i) / (29.0 * W_i)) + P);
    result.push_back(p_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_6__x_i(double M, double P, double P_i_0, double W_air, double W_i, double p_c) {
    std::vector<double> result;
    double x_i = (((29.0 * W_i) * (P - p_c)) / ((M * P_i_0) * W_air));
    result.push_back(x_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__M(double P, double P_i_0, double W_air, double W_i, double epsilon_i, double p_c, double x_i) {
    std::vector<double> result;
    double M = (((29.0 * W_i) * (P - p_c)) / (((P_i_0 * W_air) * epsilon_i) * x_i));
    result.push_back(M);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__P(double M, double P_i_0, double W_air, double W_i, double epsilon_i, double p_c, double x_i) {
    std::vector<double> result;
    double P = ((((((M * P_i_0) * W_air) * epsilon_i) * x_i) / (29.0 * W_i)) + p_c);
    result.push_back(P);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__P_i_0(double M, double P, double W_air, double W_i, double epsilon_i, double p_c, double x_i) {
    std::vector<double> result;
    double P_i_0 = (((29.0 * W_i) * (P - p_c)) / (((M * W_air) * epsilon_i) * x_i));
    result.push_back(P_i_0);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__W_air(double M, double P, double P_i_0, double W_i, double epsilon_i, double p_c, double x_i) {
    std::vector<double> result;
    double W_air = (((29.0 * W_i) * (P - p_c)) / (((M * P_i_0) * epsilon_i) * x_i));
    result.push_back(W_air);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__W_i(double M, double P, double P_i_0, double W_air, double epsilon_i, double p_c, double x_i) {
    std::vector<double> result;
    double W_i = (((((M * P_i_0) * W_air) * epsilon_i) * x_i) / (29.0 * (P - p_c)));
    result.push_back(W_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__epsilon_i(double M, double P, double P_i_0, double W_air, double W_i, double p_c, double x_i) {
    std::vector<double> result;
    double epsilon_i = (((29.0 * W_i) * (P - p_c)) / (((M * P_i_0) * W_air) * x_i));
    result.push_back(epsilon_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__p_c(double M, double P, double P_i_0, double W_air, double W_i, double epsilon_i, double x_i) {
    std::vector<double> result;
    double p_c = (((((((-M) * P_i_0) * W_air) * epsilon_i) * x_i) / (29.0 * W_i)) + P);
    result.push_back(p_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_7__x_i(double M, double P, double P_i_0, double W_air, double W_i, double epsilon_i, double p_c) {
    std::vector<double> result;
    double x_i = (((29.0 * W_i) * (P - p_c)) / (((M * P_i_0) * W_air) * epsilon_i));
    result.push_back(x_i);
    return result;
}
std::vector<double> Precondensors_eqn_7_8__L_c(double Q_condensor_heat_duty, double c_p, double del_T) {
    std::vector<double> result;
    double L_c = (Q_condensor_heat_duty / (c_p * del_T));
    result.push_back(L_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_8__Q_condensor_heat_duty(double L_c, double c_p, double del_T) {
    std::vector<double> result;
    double Q_condensor_heat_duty = ((L_c * c_p) * del_T);
    result.push_back(Q_condensor_heat_duty);
    return result;
}
std::vector<double> Precondensors_eqn_7_8__c_p(double L_c, double Q_condensor_heat_duty, double del_T) {
    std::vector<double> result;
    double c_p = (Q_condensor_heat_duty / (L_c * del_T));
    result.push_back(c_p);
    return result;
}
std::vector<double> Precondensors_eqn_7_8__del_T(double L_c, double Q_condensor_heat_duty, double c_p) {
    std::vector<double> result;
    double del_T = (Q_condensor_heat_duty / (L_c * c_p));
    result.push_back(del_T);
    return result;
}
std::vector<double> Precondensors_eqn_7_9__L_c(double Q_condensor_heat_duty, double c_p, double del_T, double rho) {
    std::vector<double> result;
    double L_c = ((0.124688279301746 * Q_condensor_heat_duty) / ((c_p * del_T) * rho));
    result.push_back(L_c);
    return result;
}
std::vector<double> Precondensors_eqn_7_9__Q_condensor_heat_duty(double L_c, double c_p, double del_T, double rho) {
    std::vector<double> result;
    double Q_condensor_heat_duty = ((((8.02 * L_c) * c_p) * del_T) * rho);
    result.push_back(Q_condensor_heat_duty);
    return result;
}
std::vector<double> Precondensors_eqn_7_9__c_p(double L_c, double Q_condensor_heat_duty, double del_T, double rho) {
    std::vector<double> result;
    double c_p = ((0.124688279301746 * Q_condensor_heat_duty) / ((L_c * del_T) * rho));
    result.push_back(c_p);
    return result;
}
std::vector<double> Precondensors_eqn_7_9__del_T(double L_c, double Q_condensor_heat_duty, double c_p, double rho) {
    std::vector<double> result;
    double del_T = ((0.124688279301746 * Q_condensor_heat_duty) / ((L_c * c_p) * rho));
    result.push_back(del_T);
    return result;
}
std::vector<double> Precondensors_eqn_7_9__rho(double L_c, double Q_condensor_heat_duty, double c_p, double del_T) {
    std::vector<double> result;
    double rho = ((0.124688279301746 * Q_condensor_heat_duty) / ((L_c * c_p) * del_T));
    result.push_back(rho);
    return result;
}
std::vector<double> PressMgmt_eqn_3_1__Abs_Pressure(double BarometricPressure, double Vacuum) {
    std::vector<double> result;
    double Abs_Pressure = (BarometricPressure - Vacuum);
    result.push_back(Abs_Pressure);
    return result;
}
std::vector<double> PressMgmt_eqn_3_1__BarometricPressure(double Abs_Pressure, double Vacuum) {
    std::vector<double> result;
    double BarometricPressure = (Abs_Pressure + Vacuum);
    result.push_back(BarometricPressure);
    return result;
}
std::vector<double> PressMgmt_eqn_3_1__Vacuum(double Abs_Pressure, double BarometricPressure) {
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
std::vector<double> PressMgmt_eqn_3_13__H_1(double H_2, double KAPPA_2, double P) {
    std::vector<double> result;
    double H_1 = (H_2 - (P / KAPPA_2));
    result.push_back(H_1);
    return result;
}
std::vector<double> PressMgmt_eqn_3_13__H_2(double H_1, double KAPPA_2, double P) {
    std::vector<double> result;
    double H_2 = (H_1 + (P / KAPPA_2));
    result.push_back(H_2);
    return result;
}
std::vector<double> PressMgmt_eqn_3_13__KAPPA_2(double H_1, double H_2, double P) {
    std::vector<double> result;
    double KAPPA_2 = ((-P) / (H_1 - H_2));
    result.push_back(KAPPA_2);
    return result;
}
std::vector<double> PressMgmt_eqn_3_13__P(double H_1, double H_2, double KAPPA_2) {
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
std::vector<double> PressMgmt_eqn_3_2__G(double G_C, double H, double P, double rho) {
    std::vector<double> result;
    double G = (((G_C * H) * P) * rho);
    result.push_back(G);
    return result;
}
std::vector<double> PressMgmt_eqn_3_2__G_C(double G, double H, double P, double rho) {
    std::vector<double> result;
    double G_C = (G / ((H * P) * rho));
    result.push_back(G_C);
    return result;
}
std::vector<double> PressMgmt_eqn_3_2__H(double G, double G_C, double P, double rho) {
    std::vector<double> result;
    double H = (G / ((G_C * P) * rho));
    result.push_back(H);
    return result;
}
std::vector<double> PressMgmt_eqn_3_2__P(double G, double G_C, double H, double rho) {
    std::vector<double> result;
    double P = (G / ((G_C * H) * rho));
    result.push_back(P);
    return result;
}
std::vector<double> PressMgmt_eqn_3_2__rho(double G, double G_C, double H, double P) {
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
std::vector<double> PressMgmt_eqn_3_6__H_1(double H_2, double P, double V, double V_P) {
    std::vector<double> result;
    double H_1 = ((H_2 - ((P * V) / V_P)) + P);
    result.push_back(H_1);
    return result;
}
std::vector<double> PressMgmt_eqn_3_6__H_2(double H_1, double P, double V, double V_P) {
    std::vector<double> result;
    double H_2 = ((H_1 + ((P * V) / V_P)) - P);
    result.push_back(H_2);
    return result;
}
std::vector<double> PressMgmt_eqn_3_6__P(double H_1, double H_2, double V, double V_P) {
    std::vector<double> result;
    double P = ((V_P * ((-H_1) + H_2)) / (V - V_P));
    result.push_back(P);
    return result;
}
std::vector<double> PressMgmt_eqn_3_6__V(double H_1, double H_2, double P, double V_P) {
    std::vector<double> result;
    double V = ((V_P * (((-H_1) + H_2) + P)) / P);
    result.push_back(V);
    return result;
}
std::vector<double> PressMgmt_eqn_3_6__V_P(double H_1, double H_2, double P, double V) {
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
std::vector<double> PressMgmt_eqn_3_9__A_C(double H_1, double H_2, double P, double V) {
    std::vector<double> result;
    double A_C = ((P * V) / (H_2 * (((-H_1) + H_2) + P)));
    result.push_back(A_C);
    return result;
}
std::vector<double> PressMgmt_eqn_3_9__H_1(double A_C, double H_2, double P, double V) {
    std::vector<double> result;
    double H_1 = ((H_2 + P) - ((P * V) / (A_C * H_2)));
    result.push_back(H_1);
    return result;
}
std::vector<double> PressMgmt_eqn_3_9__H_2(double A_C, double H_1, double P, double V) {
    std::vector<double> result;
    double H_2 = (((A_C * (H_1 - P)) - std::sqrt((A_C * ((((A_C * std::pow(H_1, 2.0)) - (((2.0 * A_C) * H_1) * P)) + (A_C * std::pow(P, 2.0))) + ((4.0 * P) * V))))) / (2.0 * A_C));
    result.push_back(H_2);
    H_2 = (((A_C * (H_1 - P)) + std::sqrt((A_C * ((((A_C * std::pow(H_1, 2.0)) - (((2.0 * A_C) * H_1) * P)) + (A_C * std::pow(P, 2.0))) + ((4.0 * P) * V))))) / (2.0 * A_C));
    result.push_back(H_2);
    return result;
}
std::vector<double> PressMgmt_eqn_3_9__P(double A_C, double H_1, double H_2, double V) {
    std::vector<double> result;
    double P = (((A_C * H_2) * (H_1 - H_2)) / ((A_C * H_2) - V));
    result.push_back(P);
    return result;
}
std::vector<double> PressMgmt_eqn_3_9__V(double A_C, double H_1, double H_2, double P) {
    std::vector<double> result;
    double V = (((A_C * H_2) * (((-H_1) + H_2) + P)) / P);
    result.push_back(V);
    return result;
}
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
std::vector<double> ProcessApp1_eqn_5_12__Eff(double N_ES, double N_t, double T) {
    std::vector<double> result;
    double Eff = std::pow((N_ES / N_t), (1.0 / T));
    result.push_back(Eff);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_12__N_ES(double Eff, double N_t, double T) {
    std::vector<double> result;
    double N_ES = (std::pow(Eff, T) * N_t);
    result.push_back(N_ES);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_12__N_t(double Eff, double N_ES, double T) {
    std::vector<double> result;
    double N_t = (N_ES / std::pow(Eff, T));
    result.push_back(N_t);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_12__T(double Eff, double N_ES, double N_t) {
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
    double M = (((294.213699178261 * T) * std::pow(W_E, 2.0)) / std::pow(P_0, 2.0));
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
std::vector<double> ProcessApp1_eqn_5_15__M_1(double M_2, double P_0_1, double P_0_2, double a_M_12) {
    std::vector<double> result;
    double M_1 = ((-M_2) / std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
    result.push_back(M_1);
    M_1 = (M_2 / std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
    result.push_back(M_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_15__M_2(double M_1, double P_0_1, double P_0_2, double a_M_12) {
    std::vector<double> result;
    double M_2 = ((-M_1) * std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
    result.push_back(M_2);
    M_2 = (M_1 * std::pow(((P_0_2 * a_M_12) / P_0_1), (5.0 / 2.0)));
    result.push_back(M_2);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_15__P_0_1(double M_1, double M_2, double P_0_2, double a_M_12) {
    std::vector<double> result;
    double P_0_1 = ((P_0_2 * a_M_12) / std::pow((M_2 / M_1), (2.0 / 5.0)));
    result.push_back(P_0_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_15__P_0_2(double M_1, double M_2, double P_0_1, double a_M_12) {
    std::vector<double> result;
    double P_0_2 = ((P_0_1 * std::pow((M_2 / M_1), (2.0 / 5.0))) / a_M_12);
    result.push_back(P_0_2);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_15__a_M_12(double M_1, double M_2, double P_0_1, double P_0_2) {
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
std::vector<double> ProcessApp1_eqn_5_17__H_2_1(double H_2_3, double H_2_mi, double x_1, double x_3) {
    std::vector<double> result;
    double H_2_1 = std::exp(((((-x_3) * std::log(H_2_3)) + std::log(H_2_mi)) / x_1));
    result.push_back(H_2_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_17__H_2_3(double H_2_1, double H_2_mi, double x_1, double x_3) {
    std::vector<double> result;
    double H_2_3 = std::exp(((((-x_1) * std::log(H_2_1)) + std::log(H_2_mi)) / x_3));
    result.push_back(H_2_3);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_17__H_2_mi(double H_2_1, double H_2_3, double x_1, double x_3) {
    std::vector<double> result;
    double H_2_mi = std::exp(((x_1 * std::log(H_2_1)) + (x_3 * std::log(H_2_3))));
    result.push_back(H_2_mi);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_17__x_1(double H_2_1, double H_2_3, double H_2_mi, double x_3) {
    std::vector<double> result;
    double x_1 = ((((-x_3) * std::log(H_2_3)) + std::log(H_2_mi)) / std::log(H_2_1));
    result.push_back(x_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_17__x_3(double H_2_1, double H_2_3, double H_2_mi, double x_1) {
    std::vector<double> result;
    double x_3 = ((((-x_1) * std::log(H_2_1)) + std::log(H_2_mi)) / std::log(H_2_3));
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
std::vector<double> ProcessApp1_eqn_5_2b__K_1(double K_2, double x_1, double x_2, double y_1, double y_2) {
    std::vector<double> result;
    double K_1 = (((K_2 * x_2) * y_1) / (x_1 * y_2));
    result.push_back(K_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_2b__K_2(double K_1, double x_1, double x_2, double y_1, double y_2) {
    std::vector<double> result;
    double K_2 = (((K_1 * x_1) * y_2) / (x_2 * y_1));
    result.push_back(K_2);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_2b__x_1(double K_1, double K_2, double x_2, double y_1, double y_2) {
    std::vector<double> result;
    double x_1 = (((K_2 * x_2) * y_1) / (K_1 * y_2));
    result.push_back(x_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_2b__x_2(double K_1, double K_2, double x_1, double y_1, double y_2) {
    std::vector<double> result;
    double x_2 = (((K_1 * x_1) * y_2) / (K_2 * y_1));
    result.push_back(x_2);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_2b__y_1(double K_1, double K_2, double x_1, double x_2, double y_2) {
    std::vector<double> result;
    double y_1 = (((K_1 * x_1) * y_2) / (K_2 * x_2));
    result.push_back(y_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_2b__y_2(double K_1, double K_2, double x_1, double x_2, double y_1) {
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
std::vector<double> ProcessApp1_eqn_5_4__P(double P_0_i, double x_i, double y_i) {
    std::vector<double> result;
    double P = ((P_0_i * x_i) / y_i);
    result.push_back(P);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_4__P_0_i(double P, double x_i, double y_i) {
    std::vector<double> result;
    double P_0_i = ((P * y_i) / x_i);
    result.push_back(P_0_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_4__x_i(double P, double P_0_i, double y_i) {
    std::vector<double> result;
    double x_i = ((P * y_i) / P_0_i);
    result.push_back(x_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_4__y_i(double P, double P_0_i, double x_i) {
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
std::vector<double> ProcessApp1_eqn_5_6__P_0_i(double gamma_i, double p_i, double x_i) {
    std::vector<double> result;
    double P_0_i = (p_i / (gamma_i * x_i));
    result.push_back(P_0_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_6__gamma_i(double P_0_i, double p_i, double x_i) {
    std::vector<double> result;
    double gamma_i = (p_i / (P_0_i * x_i));
    result.push_back(gamma_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_6__p_i(double P_0_i, double gamma_i, double x_i) {
    std::vector<double> result;
    double p_i = ((P_0_i * gamma_i) * x_i);
    result.push_back(p_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_6__x_i(double P_0_i, double gamma_i, double p_i) {
    std::vector<double> result;
    double x_i = (p_i / (P_0_i * gamma_i));
    result.push_back(x_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_7__P(double P_0_i, double gamma_i, double x_i, double y_i) {
    std::vector<double> result;
    double P = (((P_0_i * gamma_i) * x_i) / y_i);
    result.push_back(P);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_7__P_0_i(double P, double gamma_i, double x_i, double y_i) {
    std::vector<double> result;
    double P_0_i = ((P * y_i) / (gamma_i * x_i));
    result.push_back(P_0_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_7__gamma_i(double P, double P_0_i, double x_i, double y_i) {
    std::vector<double> result;
    double gamma_i = ((P * y_i) / (P_0_i * x_i));
    result.push_back(gamma_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_7__x_i(double P, double P_0_i, double gamma_i, double y_i) {
    std::vector<double> result;
    double x_i = ((P * y_i) / (P_0_i * gamma_i));
    result.push_back(x_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_7__y_i(double P, double P_0_i, double gamma_i, double x_i) {
    std::vector<double> result;
    double y_i = (((P_0_i * gamma_i) * x_i) / P);
    result.push_back(y_i);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_8__P_0_1(double P_0_2, double alpha_12, double gamma_1, double gamma_2) {
    std::vector<double> result;
    double P_0_1 = (((P_0_2 * alpha_12) * gamma_2) / gamma_1);
    result.push_back(P_0_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_8__P_0_2(double P_0_1, double alpha_12, double gamma_1, double gamma_2) {
    std::vector<double> result;
    double P_0_2 = ((P_0_1 * gamma_1) / (alpha_12 * gamma_2));
    result.push_back(P_0_2);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_8__alpha_12(double P_0_1, double P_0_2, double gamma_1, double gamma_2) {
    std::vector<double> result;
    double alpha_12 = ((P_0_1 * gamma_1) / (P_0_2 * gamma_2));
    result.push_back(alpha_12);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_8__gamma_1(double P_0_1, double P_0_2, double alpha_12, double gamma_2) {
    std::vector<double> result;
    double gamma_1 = (((P_0_2 * alpha_12) * gamma_2) / P_0_1);
    result.push_back(gamma_1);
    return result;
}
std::vector<double> ProcessApp1_eqn_5_8__gamma_2(double P_0_1, double P_0_2, double alpha_12, double gamma_1) {
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
std::vector<double> ProcessApp2_eqn_6_1__T_1(double T_2, double T_R, double c_p, double del_h_v, double w_1, double w_2, double w_v) {
    std::vector<double> result;
    double T_1 = (((((T_R * c_p) * w_1) + ((c_p * w_2) * ((-T_2) + T_R))) + (del_h_v * w_v)) / (c_p * w_1));
    result.push_back(T_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_1__T_2(double T_1, double T_R, double c_p, double del_h_v, double w_1, double w_2, double w_v) {
    std::vector<double> result;
    double T_2 = (((((T_R * c_p) * w_2) + ((c_p * w_1) * ((-T_1) + T_R))) + (del_h_v * w_v)) / (c_p * w_2));
    result.push_back(T_2);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_1__T_R(double T_1, double T_2, double c_p, double del_h_v, double w_1, double w_2, double w_v) {
    std::vector<double> result;
    double T_R = (((((T_1 * c_p) * w_1) + ((T_2 * c_p) * w_2)) - (del_h_v * w_v)) / (c_p * (w_1 + w_2)));
    result.push_back(T_R);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_1__c_p(double T_1, double T_2, double T_R, double del_h_v, double w_1, double w_2, double w_v) {
    std::vector<double> result;
    double c_p = ((del_h_v * w_v) / ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2)));
    result.push_back(c_p);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_1__del_h_v(double T_1, double T_2, double T_R, double c_p, double w_1, double w_2, double w_v) {
    std::vector<double> result;
    double del_h_v = ((c_p * ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2))) / w_v);
    result.push_back(del_h_v);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_1__w_1(double T_1, double T_2, double T_R, double c_p, double del_h_v, double w_2, double w_v) {
    std::vector<double> result;
    double w_1 = ((((((-T_2) * c_p) * w_2) + ((T_R * c_p) * w_2)) + (del_h_v * w_v)) / (c_p * (T_1 - T_R)));
    result.push_back(w_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_1__w_2(double T_1, double T_2, double T_R, double c_p, double del_h_v, double w_1, double w_v) {
    std::vector<double> result;
    double w_2 = ((((((-T_1) * c_p) * w_1) + ((T_R * c_p) * w_1)) + (del_h_v * w_v)) / (c_p * (T_2 - T_R)));
    result.push_back(w_2);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_1__w_v(double T_1, double T_2, double T_R, double c_p, double del_h_v, double w_1, double w_2) {
    std::vector<double> result;
    double w_v = ((c_p * ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2))) / del_h_v);
    result.push_back(w_v);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_10__A(double dV_dt, double delta_P, double mu, double r_c, double s, double tau) {
    std::vector<double> result;
    double A = ((((dV_dt * std::pow(delta_P, (s - 1.0))) * mu) * r_c) * tau);
    result.push_back(A);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_10__dV_dt(double A, double delta_P, double mu, double r_c, double s, double tau) {
    std::vector<double> result;
    double dV_dt = ((A * std::pow(delta_P, (1.0 - s))) / ((mu * r_c) * tau));
    result.push_back(dV_dt);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_10__delta_P(double A, double dV_dt, double mu, double r_c, double s, double tau) {
    std::vector<double> result;
    double delta_P = std::pow(((((dV_dt * mu) * r_c) * tau) / A), ((-1.0) / (s - 1.0)));
    result.push_back(delta_P);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_10__mu(double A, double dV_dt, double delta_P, double r_c, double s, double tau) {
    std::vector<double> result;
    double mu = ((A * std::pow(delta_P, (1.0 - s))) / ((dV_dt * r_c) * tau));
    result.push_back(mu);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_10__r_c(double A, double dV_dt, double delta_P, double mu, double s, double tau) {
    std::vector<double> result;
    double r_c = ((A * std::pow(delta_P, (1.0 - s))) / ((dV_dt * mu) * tau));
    result.push_back(r_c);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_10__s(double A, double dV_dt, double delta_P, double mu, double r_c, double tau) {
    std::vector<double> result;
    double s = (std::log(((A * delta_P) / (((dV_dt * mu) * r_c) * tau))) / std::log(delta_P));
    result.push_back(s);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_10__tau(double A, double dV_dt, double delta_P, double mu, double r_c, double s) {
    std::vector<double> result;
    double tau = ((A * std::pow(delta_P, (1.0 - s))) / ((dV_dt * mu) * r_c));
    result.push_back(tau);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_11a__A_d(double delta_T, double delta_h_i, double delta_m, double h_d, double m_b, double t_R) {
    std::vector<double> result;
    double A_d = (((delta_h_i * delta_m) * m_b) / ((delta_T * h_d) * t_R));
    result.push_back(A_d);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_11a__delta_T(double A_d, double delta_h_i, double delta_m, double h_d, double m_b, double t_R) {
    std::vector<double> result;
    double delta_T = (((delta_h_i * delta_m) * m_b) / ((A_d * h_d) * t_R));
    result.push_back(delta_T);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_11a__delta_h_i(double A_d, double delta_T, double delta_m, double h_d, double m_b, double t_R) {
    std::vector<double> result;
    double delta_h_i = ((((A_d * delta_T) * h_d) * t_R) / (delta_m * m_b));
    result.push_back(delta_h_i);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_11a__delta_m(double A_d, double delta_T, double delta_h_i, double h_d, double m_b, double t_R) {
    std::vector<double> result;
    double delta_m = ((((A_d * delta_T) * h_d) * t_R) / (delta_h_i * m_b));
    result.push_back(delta_m);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_11a__h_d(double A_d, double delta_T, double delta_h_i, double delta_m, double m_b, double t_R) {
    std::vector<double> result;
    double h_d = (((delta_h_i * delta_m) * m_b) / ((A_d * delta_T) * t_R));
    result.push_back(h_d);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_11a__m_b(double A_d, double delta_T, double delta_h_i, double delta_m, double h_d, double t_R) {
    std::vector<double> result;
    double m_b = ((((A_d * delta_T) * h_d) * t_R) / (delta_h_i * delta_m));
    result.push_back(m_b);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_11a__t_R(double A_d, double delta_T, double delta_h_i, double delta_m, double h_d, double m_b) {
    std::vector<double> result;
    double t_R = (((delta_h_i * delta_m) * m_b) / ((A_d * delta_T) * h_d));
    result.push_back(t_R);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_2__Q_v(double T_1, double T_2, double T_R, double c_p, double w_1, double w_2) {
    std::vector<double> result;
    double Q_v = ((c_p * ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2))) / 12000.0);
    result.push_back(Q_v);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_2__T_1(double Q_v, double T_2, double T_R, double c_p, double w_1, double w_2) {
    std::vector<double> result;
    double T_1 = ((((12000.0 * Q_v) + ((T_R * c_p) * w_1)) + ((c_p * w_2) * ((-T_2) + T_R))) / (c_p * w_1));
    result.push_back(T_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_2__T_2(double Q_v, double T_1, double T_R, double c_p, double w_1, double w_2) {
    std::vector<double> result;
    double T_2 = ((((12000.0 * Q_v) + ((T_R * c_p) * w_2)) + ((c_p * w_1) * ((-T_1) + T_R))) / (c_p * w_2));
    result.push_back(T_2);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_2__T_R(double Q_v, double T_1, double T_2, double c_p, double w_1, double w_2) {
    std::vector<double> result;
    double T_R = (((((-12000.0) * Q_v) + ((T_1 * c_p) * w_1)) + ((T_2 * c_p) * w_2)) / (c_p * (w_1 + w_2)));
    result.push_back(T_R);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_2__c_p(double Q_v, double T_1, double T_2, double T_R, double w_1, double w_2) {
    std::vector<double> result;
    double c_p = ((12000.0 * Q_v) / ((((T_1 * w_1) + (T_2 * w_2)) - (T_R * w_1)) - (T_R * w_2)));
    result.push_back(c_p);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_2__w_1(double Q_v, double T_1, double T_2, double T_R, double c_p, double w_2) {
    std::vector<double> result;
    double w_1 = ((((12000.0 * Q_v) - ((T_2 * c_p) * w_2)) + ((T_R * c_p) * w_2)) / (c_p * (T_1 - T_R)));
    result.push_back(w_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_2__w_2(double Q_v, double T_1, double T_2, double T_R, double c_p, double w_1) {
    std::vector<double> result;
    double w_2 = ((((12000.0 * Q_v) - ((T_1 * c_p) * w_1)) + ((T_R * c_p) * w_1)) / (c_p * (T_2 - T_R)));
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
std::vector<double> ProcessApp2_eqn_6_7__C_1(double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double m_b, double m_v) {
    std::vector<double> result;
    double C_1 = (((((C_2 * delta_h_c) * m_b) + ((c_p * m_b) * ((-T_1) + T_2))) + (delta_h_v * m_v)) / (delta_h_c * m_b));
    result.push_back(C_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__C_2(double C_1, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double m_b, double m_v) {
    std::vector<double> result;
    double C_2 = (((((C_1 * delta_h_c) * m_b) + ((c_p * m_b) * (T_1 - T_2))) - (delta_h_v * m_v)) / (delta_h_c * m_b));
    result.push_back(C_2);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__T_1(double C_1, double C_2, double T_2, double c_p, double delta_h_c, double delta_h_v, double m_b, double m_v) {
    std::vector<double> result;
    double T_1 = (((((T_2 * c_p) * m_b) + ((delta_h_c * m_b) * ((-C_1) + C_2))) + (delta_h_v * m_v)) / (c_p * m_b));
    result.push_back(T_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__T_2(double C_1, double C_2, double T_1, double c_p, double delta_h_c, double delta_h_v, double m_b, double m_v) {
    std::vector<double> result;
    double T_2 = (((((T_1 * c_p) * m_b) + ((delta_h_c * m_b) * (C_1 - C_2))) - (delta_h_v * m_v)) / (c_p * m_b));
    result.push_back(T_2);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__c_p(double C_1, double C_2, double T_1, double T_2, double delta_h_c, double delta_h_v, double m_b, double m_v) {
    std::vector<double> result;
    double c_p = ((((((-C_1) * delta_h_c) * m_b) + ((C_2 * delta_h_c) * m_b)) + (delta_h_v * m_v)) / (m_b * (T_1 - T_2)));
    result.push_back(c_p);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__delta_h_c(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_v, double m_b, double m_v) {
    std::vector<double> result;
    double delta_h_c = ((((((-T_1) * c_p) * m_b) + ((T_2 * c_p) * m_b)) + (delta_h_v * m_v)) / (m_b * (C_1 - C_2)));
    result.push_back(delta_h_c);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__delta_h_v(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_c, double m_b, double m_v) {
    std::vector<double> result;
    double delta_h_v = ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p))) / m_v);
    result.push_back(delta_h_v);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__m_b(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double m_v) {
    std::vector<double> result;
    double m_b = ((delta_h_v * m_v) / ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p)));
    result.push_back(m_b);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_7__m_v(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double m_b) {
    std::vector<double> result;
    double m_v = ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p))) / delta_h_v);
    result.push_back(m_v);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__C_1(double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double delta_t, double m_b, double w_v) {
    std::vector<double> result;
    double C_1 = (((((C_2 * delta_h_c) * m_b) + ((c_p * m_b) * ((-T_1) + T_2))) + ((delta_h_v * delta_t) * w_v)) / (delta_h_c * m_b));
    result.push_back(C_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__C_2(double C_1, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double delta_t, double m_b, double w_v) {
    std::vector<double> result;
    double C_2 = (((((C_1 * delta_h_c) * m_b) + ((c_p * m_b) * (T_1 - T_2))) - ((delta_h_v * delta_t) * w_v)) / (delta_h_c * m_b));
    result.push_back(C_2);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__T_1(double C_1, double C_2, double T_2, double c_p, double delta_h_c, double delta_h_v, double delta_t, double m_b, double w_v) {
    std::vector<double> result;
    double T_1 = (((((T_2 * c_p) * m_b) + ((delta_h_c * m_b) * ((-C_1) + C_2))) + ((delta_h_v * delta_t) * w_v)) / (c_p * m_b));
    result.push_back(T_1);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__T_2(double C_1, double C_2, double T_1, double c_p, double delta_h_c, double delta_h_v, double delta_t, double m_b, double w_v) {
    std::vector<double> result;
    double T_2 = (((((T_1 * c_p) * m_b) + ((delta_h_c * m_b) * (C_1 - C_2))) - ((delta_h_v * delta_t) * w_v)) / (c_p * m_b));
    result.push_back(T_2);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__c_p(double C_1, double C_2, double T_1, double T_2, double delta_h_c, double delta_h_v, double delta_t, double m_b, double w_v) {
    std::vector<double> result;
    double c_p = ((((((-C_1) * delta_h_c) * m_b) + ((C_2 * delta_h_c) * m_b)) + ((delta_h_v * delta_t) * w_v)) / (m_b * (T_1 - T_2)));
    result.push_back(c_p);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__delta_h_c(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_v, double delta_t, double m_b, double w_v) {
    std::vector<double> result;
    double delta_h_c = ((((((-T_1) * c_p) * m_b) + ((T_2 * c_p) * m_b)) + ((delta_h_v * delta_t) * w_v)) / (m_b * (C_1 - C_2)));
    result.push_back(delta_h_c);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__delta_h_v(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_t, double m_b, double w_v) {
    std::vector<double> result;
    double delta_h_v = ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p))) / (delta_t * w_v));
    result.push_back(delta_h_v);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__delta_t(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double m_b, double w_v) {
    std::vector<double> result;
    double delta_t = ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p))) / (delta_h_v * w_v));
    result.push_back(delta_t);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__m_b(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double delta_t, double w_v) {
    std::vector<double> result;
    double m_b = (((delta_h_v * delta_t) * w_v) / ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p)));
    result.push_back(m_b);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_8__w_v(double C_1, double C_2, double T_1, double T_2, double c_p, double delta_h_c, double delta_h_v, double delta_t, double m_b) {
    std::vector<double> result;
    double w_v = ((m_b * ((((C_1 * delta_h_c) - (C_2 * delta_h_c)) + (T_1 * c_p)) - (T_2 * c_p))) / (delta_h_v * delta_t));
    result.push_back(w_v);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_9__A(double dV_dt, double delta_P, double m, double mu, double r, double r_M) {
    std::vector<double> result;
    double A = (((dV_dt * r_M) - std::sqrt((dV_dt * ((dV_dt * std::pow(r_M, 2.0)) + ((((4.0 * std::pow(delta_P, 2.0)) * m) * mu) * r))))) / (2.0 * delta_P));
    result.push_back(A);
    A = (((dV_dt * r_M) + std::sqrt((dV_dt * ((dV_dt * std::pow(r_M, 2.0)) + ((((4.0 * std::pow(delta_P, 2.0)) * m) * mu) * r))))) / (2.0 * delta_P));
    result.push_back(A);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_9__dV_dt(double A, double delta_P, double m, double mu, double r, double r_M) {
    std::vector<double> result;
    double dV_dt = ((std::pow(A, 2.0) * delta_P) / ((A * r_M) + (((delta_P * m) * mu) * r)));
    result.push_back(dV_dt);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_9__delta_P(double A, double dV_dt, double m, double mu, double r, double r_M) {
    std::vector<double> result;
    double delta_P = (((A * dV_dt) * r_M) / (std::pow(A, 2.0) - (((dV_dt * m) * mu) * r)));
    result.push_back(delta_P);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_9__m(double A, double dV_dt, double delta_P, double mu, double r, double r_M) {
    std::vector<double> result;
    double m = ((A * ((A * delta_P) - (dV_dt * r_M))) / (((dV_dt * delta_P) * mu) * r));
    result.push_back(m);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_9__mu(double A, double dV_dt, double delta_P, double m, double r, double r_M) {
    std::vector<double> result;
    double mu = ((A * ((A * delta_P) - (dV_dt * r_M))) / (((dV_dt * delta_P) * m) * r));
    result.push_back(mu);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_9__r(double A, double dV_dt, double delta_P, double m, double mu, double r_M) {
    std::vector<double> result;
    double r = ((A * ((A * delta_P) - (dV_dt * r_M))) / (((dV_dt * delta_P) * m) * mu));
    result.push_back(r);
    return result;
}
std::vector<double> ProcessApp2_eqn_6_9__r_M(double A, double dV_dt, double delta_P, double m, double mu, double r) {
    std::vector<double> result;
    double r_M = (((A * delta_P) / dV_dt) - ((((delta_P * m) * mu) * r) / A));
    result.push_back(r_M);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_1__PS(double Q_0, double Q_external_gas_throughput, double V, double dP, double dT) {
    std::vector<double> result;
    double PS = ((Q_0 + Q_external_gas_throughput) - ((V * dP) / dT));
    result.push_back(PS);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_1__Q_0(double PS, double Q_external_gas_throughput, double V, double dP, double dT) {
    std::vector<double> result;
    double Q_0 = ((PS - Q_external_gas_throughput) + ((V * dP) / dT));
    result.push_back(Q_0);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_1__Q_external_gas_throughput(double PS, double Q_0, double V, double dP, double dT) {
    std::vector<double> result;
    double Q_external_gas_throughput = ((PS - Q_0) + ((V * dP) / dT));
    result.push_back(Q_external_gas_throughput);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_1__V(double PS, double Q_0, double Q_external_gas_throughput, double dP, double dT) {
    std::vector<double> result;
    double V = ((dT * (((-PS) + Q_0) + Q_external_gas_throughput)) / dP);
    result.push_back(V);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_1__dP(double PS, double Q_0, double Q_external_gas_throughput, double V, double dT) {
    std::vector<double> result;
    double dP = ((dT * (((-PS) + Q_0) + Q_external_gas_throughput)) / V);
    result.push_back(dP);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_1__dT(double PS, double Q_0, double Q_external_gas_throughput, double V, double dP) {
    std::vector<double> result;
    double dT = ((V * dP) / (((-PS) + Q_0) + Q_external_gas_throughput));
    result.push_back(dT);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__Q(double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double Q = ((((Q_0 + Q_external_gas_throughput) - SP_1) + (((-Q_0) + SP_2) * std::exp(((S_vol_pump_speed * t) / V)))) * std::exp((((-S_vol_pump_speed) * t) / V)));
    result.push_back(Q);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__Q_0(double Q, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double Q_0 = ((-((((Q * std::exp(((S_vol_pump_speed * t) / V))) - Q_external_gas_throughput) + SP_1) - (SP_2 * std::exp(((S_vol_pump_speed * t) / V))))) / (std::exp(((S_vol_pump_speed * t) / V)) - 1.0));
    result.push_back(Q_0);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__Q_external_gas_throughput(double Q, double Q_0, double SP_1, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double Q_external_gas_throughput = (((-Q_0) + SP_1) + (((Q + Q_0) - SP_2) * std::exp(((S_vol_pump_speed * t) / V))));
    result.push_back(Q_external_gas_throughput);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__SP_1(double Q, double Q_0, double Q_external_gas_throughput, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double SP_1 = ((Q_0 + Q_external_gas_throughput) + ((((-Q) - Q_0) + SP_2) * std::exp(((S_vol_pump_speed * t) / V))));
    result.push_back(SP_1);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__SP_2(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double SP_2 = (((((-Q_0) - Q_external_gas_throughput) + SP_1) + ((Q + Q_0) * std::exp(((S_vol_pump_speed * t) / V)))) * std::exp((((-S_vol_pump_speed) * t) / V)));
    result.push_back(SP_2);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__S_vol_pump_speed(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double V, double t) {
    std::vector<double> result;
    double S_vol_pump_speed = ((V * std::log((((Q_0 + Q_external_gas_throughput) - SP_1) / ((Q + Q_0) - SP_2)))) / t);
    result.push_back(S_vol_pump_speed);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__V(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double t) {
    std::vector<double> result;
    double V = ((S_vol_pump_speed * t) / std::log((((Q_0 + Q_external_gas_throughput) - SP_1) / ((Q + Q_0) - SP_2))));
    result.push_back(V);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_2__t(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double V) {
    std::vector<double> result;
    double t = ((V * std::log((((Q_0 + Q_external_gas_throughput) - SP_1) / ((Q + Q_0) - SP_2)))) / S_vol_pump_speed);
    result.push_back(t);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_3__F_s(double t, double t_c) {
    std::vector<double> result;
    double F_s = (t / t_c);
    result.push_back(F_s);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_3__t(double F_s, double t_c) {
    std::vector<double> result;
    double t = (F_s * t_c);
    result.push_back(t);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_3__t_c(double F_s, double t) {
    std::vector<double> result;
    double t_c = (t / F_s);
    result.push_back(t_c);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_4__p_g(double p_s, double p_v) {
    std::vector<double> result;
    double p_g = (p_s - p_v);
    result.push_back(p_g);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_4__p_s(double p_g, double p_v) {
    std::vector<double> result;
    double p_s = (p_g + p_v);
    result.push_back(p_s);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_4__p_v(double p_g, double p_s) {
    std::vector<double> result;
    double p_v = 0.0;
    result.push_back(p_v);
    p_v = ((-p_g) + p_s);
    result.push_back(p_v);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_5__P_0_v(double P_D, double p_g, double p_v_max) {
    std::vector<double> result;
    double P_0_v = ((P_D * p_v_max) / (p_g + p_v_max));
    result.push_back(P_0_v);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_5__P_D(double P_0_v, double p_g, double p_v_max) {
    std::vector<double> result;
    double P_D = ((P_0_v * (p_g + p_v_max)) / p_v_max);
    result.push_back(P_D);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_5__p_g(double P_0_v, double P_D, double p_v_max) {
    std::vector<double> result;
    double p_g = ((p_v_max * ((-P_0_v) + P_D)) / P_0_v);
    result.push_back(p_g);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_5__p_v_max(double P_0_v, double P_D, double p_g) {
    std::vector<double> result;
    double p_v_max = (((-P_0_v) * p_g) / (P_0_v - P_D));
    result.push_back(p_v_max);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__P_0_V(double P_D, double P_v_0, double S_B, double S_D, double p_b, double p_g, double p_v_max) {
    std::vector<double> result;
    double P_0_V = ((((((P_D * S_B) * p_b) + ((P_D * S_D) * p_v_max)) - ((P_v_0 * S_D) * p_g)) - ((P_v_0 * S_D) * p_v_max)) / (P_D * S_B));
    result.push_back(P_0_V);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__P_D(double P_0_V, double P_v_0, double S_B, double S_D, double p_b, double p_g, double p_v_max) {
    std::vector<double> result;
    double P_D = (((P_v_0 * S_D) * (p_g + p_v_max)) / ((((-P_0_V) * S_B) + (S_B * p_b)) + (S_D * p_v_max)));
    result.push_back(P_D);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__P_v_0(double P_0_V, double P_D, double S_B, double S_D, double p_b, double p_g, double p_v_max) {
    std::vector<double> result;
    double P_v_0 = ((P_D * ((((-P_0_V) * S_B) + (S_B * p_b)) + (S_D * p_v_max))) / (S_D * (p_g + p_v_max)));
    result.push_back(P_v_0);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__S_B(double P_0_V, double P_D, double P_v_0, double S_D, double p_b, double p_g, double p_v_max) {
    std::vector<double> result;
    double S_B = ((S_D * (((P_D * p_v_max) - (P_v_0 * p_g)) - (P_v_0 * p_v_max))) / (P_D * (P_0_V - p_b)));
    result.push_back(S_B);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__S_D(double P_0_V, double P_D, double P_v_0, double S_B, double p_b, double p_g, double p_v_max) {
    std::vector<double> result;
    double S_D = (((P_D * S_B) * (P_0_V - p_b)) / (((P_D * p_v_max) - (P_v_0 * p_g)) - (P_v_0 * p_v_max)));
    result.push_back(S_D);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__p_b(double P_0_V, double P_D, double P_v_0, double S_B, double S_D, double p_g, double p_v_max) {
    std::vector<double> result;
    double p_b = ((((((P_0_V * P_D) * S_B) - ((P_D * S_D) * p_v_max)) + ((P_v_0 * S_D) * p_g)) + ((P_v_0 * S_D) * p_v_max)) / (P_D * S_B));
    result.push_back(p_b);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__p_g(double P_0_V, double P_D, double P_v_0, double S_B, double S_D, double p_b, double p_v_max) {
    std::vector<double> result;
    double p_g = (((((((-P_0_V) * P_D) * S_B) + ((P_D * S_B) * p_b)) + ((P_D * S_D) * p_v_max)) - ((P_v_0 * S_D) * p_v_max)) / (P_v_0 * S_D));
    result.push_back(p_g);
    return result;
}
std::vector<double> RotaryPistonVane_eqn_11_6__p_v_max(double P_0_V, double P_D, double P_v_0, double S_B, double S_D, double p_b, double p_g) {
    std::vector<double> result;
    double p_v_max = (((((P_0_V * P_D) * S_B) - ((P_D * S_B) * p_b)) + ((P_v_0 * S_D) * p_g)) / (S_D * (P_D - P_v_0)));
    result.push_back(p_v_max);
    return result;
}
std::vector<double> SelectingPump_eqn_8_1__NC(double NS, double SCON, double installation_cost) {
    std::vector<double> result;
    double NC = (((-0.5) * NS) + ((0.000350630766969363 * installation_cost) / std::pow(SCON, (7.0 / 20.0))));
    result.push_back(NC);
    return result;
}
std::vector<double> SelectingPump_eqn_8_1__NS(double NC, double SCON, double installation_cost) {
    std::vector<double> result;
    double NS = (((-2.0) * NC) + ((0.000701261533938727 * installation_cost) / std::pow(SCON, (7.0 / 20.0))));
    result.push_back(NS);
    return result;
}
std::vector<double> SelectingPump_eqn_8_1__SCON(double NC, double NS, double installation_cost) {
    std::vector<double> result;
    double SCON = (1000.0 * std::pow((installation_cost / (16000.0 * (NS + (2.0 * NC)))), (1.0 / 0.35)));
    return {SCON};
}
std::vector<double> SelectingPump_eqn_8_1__installation_cost(double NC, double NS, double SCON) {
    std::vector<double> result;
    double installation_cost = ((1426.00150101399 * std::pow(SCON, (7.0 / 20.0))) * ((2.0 * NC) + NS));
    result.push_back(installation_cost);
    return result;
}
std::vector<double> SelectingPump_eqn_8_2__hp(double installed_costs) {
    std::vector<double> result;
    double hp = (9.18273645546364e-09 * std::pow(installed_costs, 2.0));
    result.push_back(hp);
    return result;
}
std::vector<double> SelectingPump_eqn_8_2__installed_costs(double hp) {
    std::vector<double> result;
    double installed_costs = (10435.5162785557 * std::sqrt(hp));
    result.push_back(installed_costs);
    return result;
}
std::vector<double> SelectingPump_eqn_8_3__hp(double installed_costs) {
    std::vector<double> result;
    double hp = (10.0 * std::pow((installed_costs / 38000.0), (1.0 / 0.45)));
    return {hp};
}
std::vector<double> SelectingPump_eqn_8_3__installed_costs(double hp) {
    std::vector<double> result;
    double installed_costs = (13482.9087908759 * std::pow(hp, (9.0 / 20.0)));
    result.push_back(installed_costs);
    return result;
}
std::vector<double> SelectingPump_eqn_8_4__hp(double installed_costs) {
    std::vector<double> result;
    double hp = ((-9.1741667595569e-11) * std::pow(installed_costs, (5.0 / 2.0)));
    result.push_back(hp);
    hp = (9.1741667595569e-11 * std::pow(installed_costs, (5.0 / 2.0)));
    result.push_back(hp);
    return result;
}
std::vector<double> SelectingPump_eqn_8_4__installed_costs(double hp) {
    std::vector<double> result;
    double installed_costs = (10350.7864343909 * std::pow(hp, (2.0 / 5.0)));
    result.push_back(installed_costs);
    return result;
}
std::vector<double> SelectingPump_eqn_8_5__Eff(double actual_brake_horsepower, double theoretical_adiabatic_horsepower) {
    std::vector<double> result;
    double Eff = (theoretical_adiabatic_horsepower / actual_brake_horsepower);
    result.push_back(Eff);
    return result;
}
std::vector<double> SelectingPump_eqn_8_5__actual_brake_horsepower(double Eff, double theoretical_adiabatic_horsepower) {
    std::vector<double> result;
    double actual_brake_horsepower = (theoretical_adiabatic_horsepower / Eff);
    result.push_back(actual_brake_horsepower);
    return result;
}
std::vector<double> SelectingPump_eqn_8_5__theoretical_adiabatic_horsepower(double Eff, double actual_brake_horsepower) {
    std::vector<double> result;
    double theoretical_adiabatic_horsepower = (Eff * actual_brake_horsepower);
    result.push_back(theoretical_adiabatic_horsepower);
    return result;
}
std::vector<double> SelectingPump_eqn_8_6__M(double P_1, double P_2, double R, double T, double adiabatic_hp, double k, double w) {
    std::vector<double> result;
    double M = (((((R * T) * k) * w) * (std::pow((P_2 / P_1), ((k - 1.0) / k)) - 1.0)) / ((1980000.0 * adiabatic_hp) * (k - 1.0)));
    result.push_back(M);
    return result;
}
std::vector<double> SelectingPump_eqn_8_6__P_1(double M, double P_2, double R, double T, double adiabatic_hp, double k, double w) {
    std::vector<double> result;
    double P_1 = (P_2 / std::pow((((((1980000.0 * M) * adiabatic_hp) / ((R * T) * w)) - (((1980000.0 * M) * adiabatic_hp) / (((R * T) * k) * w))) + 1.0), (k / (k - 1.0))));
    result.push_back(P_1);
    return result;
}
std::vector<double> SelectingPump_eqn_8_6__P_2(double M, double P_1, double R, double T, double adiabatic_hp, double k, double w) {
    std::vector<double> result;
    double P_2 = (P_1 * std::pow((((((1980000.0 * M) * adiabatic_hp) / ((R * T) * w)) - (((1980000.0 * M) * adiabatic_hp) / (((R * T) * k) * w))) + 1.0), (k / (k - 1.0))));
    result.push_back(P_2);
    return result;
}
std::vector<double> SelectingPump_eqn_8_6__R(double M, double P_1, double P_2, double T, double adiabatic_hp, double k, double w) {
    std::vector<double> result;
    double R = ((((1980000.0 * M) * adiabatic_hp) * (k - 1.0)) / (((T * k) * w) * (std::pow((P_2 / P_1), ((k - 1.0) / k)) - 1.0)));
    result.push_back(R);
    return result;
}
std::vector<double> SelectingPump_eqn_8_6__T(double M, double P_1, double P_2, double R, double adiabatic_hp, double k, double w) {
    std::vector<double> result;
    double T = ((((1980000.0 * M) * adiabatic_hp) * (k - 1.0)) / (((R * k) * w) * (std::pow((P_2 / P_1), ((k - 1.0) / k)) - 1.0)));
    result.push_back(T);
    return result;
}
std::vector<double> SelectingPump_eqn_8_6__adiabatic_hp(double M, double P_1, double P_2, double R, double T, double k, double w) {
    std::vector<double> result;
    double adiabatic_hp = (((((R * T) * k) * w) * (std::pow((P_2 / P_1), ((k - 1.0) / k)) - 1.0)) / ((1980000.0 * M) * (k - 1.0)));
    result.push_back(adiabatic_hp);
    return result;
}
std::vector<double> SelectingPump_eqn_8_6__k(double M, double P_1, double P_2, double R, double T, double adiabatic_hp, double w) {
    throw std::runtime_error("SelectingPump_eqn_8_6__k: requires numerical solver (not transpilable)");
}
std::vector<double> SelectingPump_eqn_8_6__w(double M, double P_1, double P_2, double R, double T, double adiabatic_hp, double k) {
    std::vector<double> result;
    double w = ((((1980000.0 * M) * adiabatic_hp) * (k - 1.0)) / (((R * T) * k) * (std::pow((P_2 / P_1), ((k - 1.0) / k)) - 1.0)));
    result.push_back(w);
    return result;
}
std::vector<double> SelectingPump_eqn_8_7__P_1(double P_2, double adiabatic_hp, double w) {
    std::vector<double> result;
    double P_1 = (P_2 / std::pow(((adiabatic_hp / (w / 20.0)) + 1.0), (1.0 / 0.286)));
    return {P_1};
}
std::vector<double> SelectingPump_eqn_8_7__P_2(double P_1, double adiabatic_hp, double w) {
    std::vector<double> result;
    double P_2 = (P_1 * std::pow(((adiabatic_hp / (w / 20.0)) + 1.0), (1.0 / 0.286)));
    return {P_2};
}
std::vector<double> SelectingPump_eqn_8_7__adiabatic_hp(double P_1, double P_2, double w) {
    std::vector<double> result;
    double adiabatic_hp = ((0.05 * w) * (std::pow((P_2 / P_1), (143.0 / 500.0)) - 1.0));
    result.push_back(adiabatic_hp);
    return result;
}
std::vector<double> SelectingPump_eqn_8_7__w(double P_1, double P_2, double adiabatic_hp) {
    std::vector<double> result;
    double w = ((20.0 * adiabatic_hp) / (std::pow((P_2 / P_1), 0.286) - 1.0));
    result.push_back(w);
    return result;
}
std::vector<double> SelectingPump_eqn_8_8__P_1(double P_2, double adiabatic_power_watts, double f) {
    std::vector<double> result;
    double P_1 = (P_2 / std::pow(((adiabatic_power_watts / (f / 12.0)) + 1.0), (1.0 / 0.286)));
    return {P_1};
}
std::vector<double> SelectingPump_eqn_8_8__P_2(double P_1, double adiabatic_power_watts, double f) {
    std::vector<double> result;
    double P_2 = (P_1 * std::pow(((adiabatic_power_watts / (f / 12.0)) + 1.0), (1.0 / 0.286)));
    return {P_2};
}
std::vector<double> SelectingPump_eqn_8_8__adiabatic_power_watts(double P_1, double P_2, double f) {
    std::vector<double> result;
    double adiabatic_power_watts = ((0.0833333333333333 * f) * (std::pow((P_2 / P_1), (143.0 / 500.0)) - 1.0));
    result.push_back(adiabatic_power_watts);
    return result;
}
std::vector<double> SelectingPump_eqn_8_8__f(double P_1, double P_2, double adiabatic_power_watts) {
    std::vector<double> result;
    double f = ((12.0 * adiabatic_power_watts) / (std::pow((P_2 / P_1), 0.286) - 1.0));
    result.push_back(f);
    return result;
}
std::vector<double> SelectingPump_eqn_8_9__E_j(double E_m, double e, double r, double s) {
    std::vector<double> result;
    double E_j = ((((0.341296928327645 * E_m) * r) * s) / M_E);
    result.push_back(E_j);
    return result;
}
std::vector<double> SelectingPump_eqn_8_9__E_m(double E_j, double e, double r, double s) {
    std::vector<double> result;
    double E_m = (((2.93 * E_j) * M_E) / (r * s));
    result.push_back(E_m);
    return result;
}
std::vector<double> SelectingPump_eqn_8_9__e(double E_j, double E_m, double r, double s) {
    std::vector<double> result;
    double e = ((((0.341296928327645 * E_m) * r) * s) / E_j);
    result.push_back(M_E);
    return result;
}
std::vector<double> SelectingPump_eqn_8_9__r(double E_j, double E_m, double e, double s) {
    std::vector<double> result;
    double r = (((2.93 * E_j) * M_E) / (E_m * s));
    result.push_back(r);
    return result;
}
std::vector<double> SelectingPump_eqn_8_9__s(double E_j, double E_m, double e, double r) {
    std::vector<double> result;
    double s = (((2.93 * E_j) * M_E) / (E_m * r));
    result.push_back(s);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_1__A(double rho_s, double v, double w_s) {
    std::vector<double> result;
    double A = (w_s / (rho_s * v));
    result.push_back(A);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_1__rho_s(double A, double v, double w_s) {
    std::vector<double> result;
    double rho_s = (w_s / (A * v));
    result.push_back(rho_s);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_1__v(double A, double rho_s, double w_s) {
    std::vector<double> result;
    double v = (w_s / (A * rho_s));
    result.push_back(v);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_1__w_s(double A, double rho_s, double v) {
    std::vector<double> result;
    double w_s = ((A * rho_s) * v);
    result.push_back(w_s);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_2__P_m(double d_n, double rho_s, double w_s) {
    std::vector<double> result;
    double P_m = ((1.334027668054e-06 * std::pow(w_s, 2.0)) / (std::pow(d_n, 4.0) * rho_s));
    result.push_back(P_m);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_2__d_n(double P_m, double rho_s, double w_s) {
    std::vector<double> result;
    double d_n = ((-0.0339853079285911) * std::sqrt((w_s / std::pow((P_m * rho_s), 0.5))));
    result.push_back(d_n);
    d_n = (0.0339853079285911 * std::sqrt((w_s / std::pow((P_m * rho_s), 0.5))));
    result.push_back(d_n);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_2__rho_s(double P_m, double d_n, double w_s) {
    std::vector<double> result;
    double rho_s = ((1.334027668054e-06 * std::pow(w_s, 2.0)) / (P_m * std::pow(d_n, 4.0)));
    result.push_back(rho_s);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_2__w_s(double P_m, double d_n, double rho_s) {
    std::vector<double> result;
    double w_s = ((865.8 * std::pow(d_n, 2.0)) * std::sqrt((P_m * rho_s)));
    result.push_back(w_s);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_3__P_s(double V, double t_e, double w_j) {
    std::vector<double> result;
    double P_s = ((33.3333333333333 * ((23.0 * V) - ((10.0 * t_e) * w_j))) / V);
    result.push_back(P_s);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_3__V(double P_s, double t_e, double w_j) {
    std::vector<double> result;
    double V = ((((-1000.0) * t_e) * w_j) / ((3.0 * P_s) - 2300.0));
    result.push_back(V);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_3__t_e(double P_s, double V, double w_j) {
    std::vector<double> result;
    double t_e = (((0.001 * V) * (2300.0 - (3.0 * P_s))) / w_j);
    result.push_back(t_e);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_3__w_j(double P_s, double V, double t_e) {
    std::vector<double> result;
    double w_j = (((0.001 * V) * (2300.0 - (3.0 * P_s))) / t_e);
    result.push_back(w_j);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_4__AEL(double SC, double r, double w_s) {
    std::vector<double> result;
    double AEL = (w_s / (SC * r));
    result.push_back(AEL);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_4__SC(double AEL, double r, double w_s) {
    std::vector<double> result;
    double SC = (w_s / (AEL * r));
    result.push_back(SC);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_4__r(double AEL, double SC, double w_s) {
    std::vector<double> result;
    double r = (w_s / (AEL * SC));
    result.push_back(r);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_4__w_s(double AEL, double SC, double r) {
    std::vector<double> result;
    double w_s = ((AEL * SC) * r);
    result.push_back(w_s);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_5__V(double r_h, double t_h, double w_h) {
    std::vector<double> result;
    double V = ((t_h * w_h) / r_h);
    result.push_back(V);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_5__r_h(double V, double t_h, double w_h) {
    std::vector<double> result;
    double r_h = ((t_h * w_h) / V);
    result.push_back(r_h);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_5__t_h(double V, double r_h, double w_h) {
    std::vector<double> result;
    double t_h = ((V * r_h) / w_h);
    result.push_back(t_h);
    return result;
}
std::vector<double> SteamJetInjectors_eqn_9_5__w_h(double V, double r_h, double t_h) {
    std::vector<double> result;
    double w_h = ((V * r_h) / t_h);
    result.push_back(w_h);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_10__P_1(double P_2, double T_1, double T_2, double V_1, double V_2) {
    std::vector<double> result;
    double P_1 = (((P_2 * T_1) * V_2) / (T_2 * V_1));
    result.push_back(P_1);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_10__P_2(double P_1, double T_1, double T_2, double V_1, double V_2) {
    std::vector<double> result;
    double P_2 = (((P_1 * T_2) * V_1) / (T_1 * V_2));
    result.push_back(P_2);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_10__T_1(double P_1, double P_2, double T_2, double V_1, double V_2) {
    std::vector<double> result;
    double T_1 = (((P_1 * T_2) * V_1) / (P_2 * V_2));
    result.push_back(T_1);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_10__T_2(double P_1, double P_2, double T_1, double V_1, double V_2) {
    std::vector<double> result;
    double T_2 = (((P_2 * T_1) * V_2) / (P_1 * V_1));
    result.push_back(T_2);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_10__V_1(double P_1, double P_2, double T_1, double T_2, double V_2) {
    std::vector<double> result;
    double V_1 = (((P_2 * T_1) * V_2) / (P_1 * T_2));
    result.push_back(V_1);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_10__V_2(double P_1, double P_2, double T_1, double T_2, double V_1) {
    std::vector<double> result;
    double V_2 = (((P_1 * T_2) * V_1) / (P_2 * T_1));
    result.push_back(V_2);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_11__M(double P, double T, double W, double q) {
    std::vector<double> result;
    double M = (((6821.0 * T) * W) / ((738.0 * P) * q));
    result.push_back(M);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_11__P(double M, double T, double W, double q) {
    std::vector<double> result;
    double P = (((6821.0 * T) * W) / ((738.0 * M) * q));
    result.push_back(P);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_11__T(double M, double P, double W, double q) {
    std::vector<double> result;
    double T = ((((738.0 * M) * P) * q) / (6821.0 * W));
    result.push_back(T);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_11__W(double M, double P, double T, double q) {
    std::vector<double> result;
    double W = ((((738.0 * M) * P) * q) / (6821.0 * T));
    result.push_back(W);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_11__q(double M, double P, double T, double W) {
    std::vector<double> result;
    double q = (((6821.0 * T) * W) / ((738.0 * M) * P));
    result.push_back(q);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_12__Total_P(double sum_partial_pressures) {
    std::vector<double> result;
    double Total_P = sum_partial_pressures;
    result.push_back(Total_P);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_12__sum_partial_pressures(double Total_P) {
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
std::vector<double> VacuumTheory_eqn_1_7__R(double T, double V, double n, double p) {
    std::vector<double> result;
    double R = ((V * p) / (T * n));
    result.push_back(R);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_7__T(double R, double V, double n, double p) {
    std::vector<double> result;
    double T = ((V * p) / (R * n));
    result.push_back(T);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_7__V(double R, double T, double n, double p) {
    std::vector<double> result;
    double V = (((R * T) * n) / p);
    result.push_back(V);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_7__n(double R, double T, double V, double p) {
    std::vector<double> result;
    double n = ((V * p) / (R * T));
    result.push_back(n);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_7__p(double R, double T, double V, double n) {
    std::vector<double> result;
    double p = (((R * T) * n) / V);
    result.push_back(p);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_8__M(double P, double R, double T, double V, double m) {
    std::vector<double> result;
    double M = (((R * T) * m) / (P * V));
    result.push_back(M);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_8__P(double M, double R, double T, double V, double m) {
    std::vector<double> result;
    double P = (((R * T) * m) / (M * V));
    result.push_back(P);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_8__R(double M, double P, double T, double V, double m) {
    std::vector<double> result;
    double R = (((M * P) * V) / (T * m));
    result.push_back(R);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_8__T(double M, double P, double R, double V, double m) {
    std::vector<double> result;
    double T = (((M * P) * V) / (R * m));
    result.push_back(T);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_8__V(double M, double P, double R, double T, double m) {
    std::vector<double> result;
    double V = (((R * T) * m) / (M * P));
    result.push_back(V);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_8__m(double M, double P, double R, double T, double V) {
    std::vector<double> result;
    double m = (((M * P) * V) / (R * T));
    result.push_back(m);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_9__M(double P, double R, double T, double rho) {
    std::vector<double> result;
    double M = (((R * T) * rho) / P);
    result.push_back(M);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_9__P(double M, double R, double T, double rho) {
    std::vector<double> result;
    double P = (((R * T) * rho) / M);
    result.push_back(P);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_9__R(double M, double P, double T, double rho) {
    std::vector<double> result;
    double R = ((M * P) / (T * rho));
    result.push_back(R);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_9__T(double M, double P, double R, double rho) {
    std::vector<double> result;
    double T = ((M * P) / (R * rho));
    result.push_back(T);
    return result;
}
std::vector<double> VacuumTheory_eqn_1_9__rho(double M, double P, double R, double T) {
    std::vector<double> result;
    double rho = ((M * P) / (R * T));
    result.push_back(rho);
    return result;
}

int main() {
    std::cout << "Running test suite..." << std::endl;
    int pass = 0, fail = 0;
    std::cout << "Testing AirLeak_eqn_4_10__T... ";
    try { auto r = AirLeak_eqn_4_10__T(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing AirLeak_eqn_4_10__V... ";
    try { auto r = AirLeak_eqn_4_10__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing AirLeak_eqn_4_10__del_P... ";
    try { auto r = AirLeak_eqn_4_10__del_P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing AirLeak_eqn_4_10__leakage... ";
    try { auto r = AirLeak_eqn_4_10__leakage(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing AirLeak_eqn_4_10__t... ";
    try { auto r = AirLeak_eqn_4_10__t(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing AirLeak_eqn_4_7__W... ";
    try { auto r = AirLeak_eqn_4_7__W(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing AirLeak_eqn_4_7__W_T... ";
    try { auto r = AirLeak_eqn_4_7__W_T(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing AirLeak_eqn_4_7__sum_individual_leak_rates... ";
    try { auto r = AirLeak_eqn_4_7__sum_individual_leak_rates(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__D... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_1__D(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__Re... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_1__Re(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_1__mu(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__rho... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_1__rho(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__v... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_1__v(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_10__Suc_Pres... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_10__Suc_Pres(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_10__delta_P... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_10__delta_P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_10__oper_press... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_10__oper_press(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__D... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_11__D(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_11__L(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__f... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_11__f(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__g_c... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_11__g_c(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__h_r... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_11__h_r(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__v... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_11__v(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_12__L(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__d... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_12__d(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__delta_P... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_12__delta_P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__f... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_12__f(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__g... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_12__g(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__rho... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_12__rho(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__v... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_12__v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_13__L(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__d... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_13__d(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__delta_P... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_13__delta_P(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__f... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_13__f(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__q... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_13__q(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__rho... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_13__rho(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__M... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_14__M(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__R... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_14__R(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__T... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_14__T(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__g_c... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_14__g_c(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__k... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_14__k(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__v_s... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_14__v_s(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_15__Re... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_15__Re(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_15__f... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_15__f(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_16__Re... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_16__Re(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_16__f... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_16__f(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_17__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_17__L(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_17__d... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_17__d(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_17__delta_P... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_17__delta_P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_17__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_17__mu(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_17__q... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_17__q(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_17__v... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_17__v(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_18a__D_eq... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_18a__D_eq(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_18a__R_ll... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_18a__R_ll(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_18b__R_ll... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_18b__R_ll(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_18b__h... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_18b__h(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_18b__w... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_18b__w(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__R_ll... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19a__R_ll(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__Re... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19a__Re(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19a__mu(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__rho... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19a__rho(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__v... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19a__v(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__Re... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19b__Re(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__h... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19b__h(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19b__mu(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__rho... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19b__rho(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__v... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19b__v(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__w... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_19b__w(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_2__delta... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_2__delta(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_2__lambd... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_2__lambd(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_2__psi... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_2__psi(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_20__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_20__L(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_20__sum_equivalent_length... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_20__sum_equivalent_length(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_20__sum_pipe... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_20__sum_pipe(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_22__P_s... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_22__P_s(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_22__Q_throughput... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_22__Q_throughput(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_22__S_p... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_22__S_p(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_25__C(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__P_1... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_25__P_1(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__P_2... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_25__P_2(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__Q_throughput... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_25__Q_throughput(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__D... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_26__D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_26__L(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__P_downstream... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_26__P_downstream(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__P_p... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_26__P_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__P_upstream... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_26__P_upstream(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_26__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__q... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_26__q(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_28__C(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__D... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_28__D(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_28__L(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__P_p... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_28__P_p(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_28__mu(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_29__C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_29__C(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_29__S_1... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_29__S_1(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_29__S_2... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_29__S_2(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_3__D... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_3__D(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_3__kn... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_3__kn(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_3__lambd... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_3__lambd(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_31__C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_31__C(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_31__S_p... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_31__S_p(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_31__S_pump_speed... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_31__S_pump_speed(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_32__C_series... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_32__C_series(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_32__geometric_sum_C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_32__geometric_sum_C(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_33__C_paralell... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_33__C_paralell(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_33__arithmetic_sum_C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_33__arithmetic_sum_C(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_34__C(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__C_1... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_34__C_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__C_2... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_34__C_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__D... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_34__D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_34__L(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__P_p... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_34__P_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_34__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_35__C_L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_35__C_L(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_35__C_T... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_35__C_T(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_35__F_p... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_35__F_p(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_36__C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_36__C(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_36__C_0... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_36__C_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_36__F_t... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_36__F_t(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__A... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_37__A(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__C... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_37__C(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__F_t... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_37__F_t(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__M... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_37__M(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__T... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_37__T(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_4___beta... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_4___beta(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_4__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_4__mu(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_4__vel_grad... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_4__vel_grad(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__D... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_5__D(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__L... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_5__L(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__delta_P... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_5__delta_P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_5__mu(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__q... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_5__q(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__lambd... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_6__lambd(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__mu... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_6__mu(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__rho... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_6__rho(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__v_a... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_6__v_a(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__T... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_7__T(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__k... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_7__k(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__m... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_7__m(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__v_a... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_7__v_a(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__M... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_8__M(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__P_c... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_8__P_c(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__T_c... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_8__T_c(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__mu_c... ";
    try { auto r = FluidFlowVacuumLines_eqn_2_8__mu_c(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_1__D_r... ";
    try { auto r = LiquidRing_eqn_10_1__D_r(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_1__sig_R... ";
    try { auto r = LiquidRing_eqn_10_1__sig_R(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_1__w... ";
    try { auto r = LiquidRing_eqn_10_1__w(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_10__bhp... ";
    try { auto r = LiquidRing_eqn_10_10__bhp(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_10__bhp_0... ";
    try { auto r = LiquidRing_eqn_10_10__bhp_0(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_10__mu... ";
    try { auto r = LiquidRing_eqn_10_10__mu(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_10__rho... ";
    try { auto r = LiquidRing_eqn_10_10__rho(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_11__T_c... ";
    try { auto r = LiquidRing_eqn_10_11__T_c(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_11__T_s... ";
    try { auto r = LiquidRing_eqn_10_11__T_s(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_12__T_c... ";
    try { auto r = LiquidRing_eqn_10_12__T_c(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_12__T_s... ";
    try { auto r = LiquidRing_eqn_10_12__T_s(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_13__T_c... ";
    try { auto r = LiquidRing_eqn_10_13__T_c(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_13__T_s... ";
    try { auto r = LiquidRing_eqn_10_13__T_s(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_14__T_c... ";
    try { auto r = LiquidRing_eqn_10_14__T_c(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_14__T_s... ";
    try { auto r = LiquidRing_eqn_10_14__T_s(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_15__P... ";
    try { auto r = LiquidRing_eqn_10_15__P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_15__S_Th... ";
    try { auto r = LiquidRing_eqn_10_15__S_Th(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_15__S_p... ";
    try { auto r = LiquidRing_eqn_10_15__S_p(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_15__p_s... ";
    try { auto r = LiquidRing_eqn_10_15__p_s(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_16__P... ";
    try { auto r = LiquidRing_eqn_10_16__P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_16__S_0... ";
    try { auto r = LiquidRing_eqn_10_16__S_0(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_16__S_Th... ";
    try { auto r = LiquidRing_eqn_10_16__S_Th(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_16__p_0... ";
    try { auto r = LiquidRing_eqn_10_16__p_0(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_17__P... ";
    try { auto r = LiquidRing_eqn_10_17__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_17__S_0... ";
    try { auto r = LiquidRing_eqn_10_17__S_0(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_17__S_Th... ";
    try { auto r = LiquidRing_eqn_10_17__S_Th(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_17__p_0... ";
    try { auto r = LiquidRing_eqn_10_17__p_0(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_17__p_s... ";
    try { auto r = LiquidRing_eqn_10_17__p_s(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_18__P... ";
    try { auto r = LiquidRing_eqn_10_18__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_18__S_Th... ";
    try { auto r = LiquidRing_eqn_10_18__S_Th(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_18__S_p... ";
    try { auto r = LiquidRing_eqn_10_18__S_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_18__T_e... ";
    try { auto r = LiquidRing_eqn_10_18__T_e(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_18__T_i... ";
    try { auto r = LiquidRing_eqn_10_18__T_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_18__p_c... ";
    try { auto r = LiquidRing_eqn_10_18__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_18__p_s... ";
    try { auto r = LiquidRing_eqn_10_18__p_s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_19__P... ";
    try { auto r = LiquidRing_eqn_10_19__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_19__S_Th... ";
    try { auto r = LiquidRing_eqn_10_19__S_Th(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_19__S_p... ";
    try { auto r = LiquidRing_eqn_10_19__S_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_19__T_e... ";
    try { auto r = LiquidRing_eqn_10_19__T_e(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_19__T_i... ";
    try { auto r = LiquidRing_eqn_10_19__T_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_19__p_c... ";
    try { auto r = LiquidRing_eqn_10_19__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_19__p_s... ";
    try { auto r = LiquidRing_eqn_10_19__p_s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_2__PS... ";
    try { auto r = LiquidRing_eqn_10_2__PS(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_2__Q_gas... ";
    try { auto r = LiquidRing_eqn_10_2__Q_gas(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_2__V... ";
    try { auto r = LiquidRing_eqn_10_2__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_2__dP... ";
    try { auto r = LiquidRing_eqn_10_2__dP(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_2__dt... ";
    try { auto r = LiquidRing_eqn_10_2__dt(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__P... ";
    try { auto r = LiquidRing_eqn_10_20__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__S_0... ";
    try { auto r = LiquidRing_eqn_10_20__S_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__S_p... ";
    try { auto r = LiquidRing_eqn_10_20__S_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__T_e... ";
    try { auto r = LiquidRing_eqn_10_20__T_e(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__T_i... ";
    try { auto r = LiquidRing_eqn_10_20__T_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__p_0... ";
    try { auto r = LiquidRing_eqn_10_20__p_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__p_c... ";
    try { auto r = LiquidRing_eqn_10_20__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_20__p_s... ";
    try { auto r = LiquidRing_eqn_10_20__p_s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_21__P... ";
    try { auto r = LiquidRing_eqn_10_21__P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_21__P_d... ";
    try { auto r = LiquidRing_eqn_10_21__P_d(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_21__P_prime... ";
    try { auto r = LiquidRing_eqn_10_21__P_prime(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_3__N_mfw... ";
    try { auto r = LiquidRing_eqn_10_3__N_mfw(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_3__Q_gas... ";
    try { auto r = LiquidRing_eqn_10_3__Q_gas(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_3__T... ";
    try { auto r = LiquidRing_eqn_10_3__T(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_4__Q_gas... ";
    try { auto r = LiquidRing_eqn_10_4__Q_gas(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_4__SP_1... ";
    try { auto r = LiquidRing_eqn_10_4__SP_1(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_4__SP_2... ";
    try { auto r = LiquidRing_eqn_10_4__SP_2(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_4__S_p... ";
    try { auto r = LiquidRing_eqn_10_4__S_p(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_4__V... ";
    try { auto r = LiquidRing_eqn_10_4__V(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_4__t... ";
    try { auto r = LiquidRing_eqn_10_4__t(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_5__P_1... ";
    try { auto r = LiquidRing_eqn_10_5__P_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_5__P_2... ";
    try { auto r = LiquidRing_eqn_10_5__P_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_5__S_p... ";
    try { auto r = LiquidRing_eqn_10_5__S_p(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_5__V... ";
    try { auto r = LiquidRing_eqn_10_5__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_5__t... ";
    try { auto r = LiquidRing_eqn_10_5__t(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_6__P_1... ";
    try { auto r = LiquidRing_eqn_10_6__P_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_6__P_2... ";
    try { auto r = LiquidRing_eqn_10_6__P_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_6__S_a... ";
    try { auto r = LiquidRing_eqn_10_6__S_a(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_6__V... ";
    try { auto r = LiquidRing_eqn_10_6__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_6__t... ";
    try { auto r = LiquidRing_eqn_10_6__t(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_8__bhp... ";
    try { auto r = LiquidRing_eqn_10_8__bhp(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_8__c_p... ";
    try { auto r = LiquidRing_eqn_10_8__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_8__delta_T... ";
    try { auto r = LiquidRing_eqn_10_8__delta_T(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_8__delta_h_i... ";
    try { auto r = LiquidRing_eqn_10_8__delta_h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_8__f_a... ";
    try { auto r = LiquidRing_eqn_10_8__f_a(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_8__rho... ";
    try { auto r = LiquidRing_eqn_10_8__rho(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_8__w_i... ";
    try { auto r = LiquidRing_eqn_10_8__w_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_9__T_c... ";
    try { auto r = LiquidRing_eqn_10_9__T_c(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_9__T_s... ";
    try { auto r = LiquidRing_eqn_10_9__T_s(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LiquidRing_eqn_10_9__delta_T... ";
    try { auto r = LiquidRing_eqn_10_9__delta_T(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_1__P... ";
    try { auto r = Precondensors_eqn_7_1__P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_1__p_i... ";
    try { auto r = Precondensors_eqn_7_1__p_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_1__y_i... ";
    try { auto r = Precondensors_eqn_7_1__y_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_10__L_c_P... ";
    try { auto r = Precondensors_eqn_7_10__L_c_P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_10__Q_condensor_heat_duty... ";
    try { auto r = Precondensors_eqn_7_10__Q_condensor_heat_duty(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_10__del_T... ";
    try { auto r = Precondensors_eqn_7_10__del_T(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_11__Q_condensor_heat_duty... ";
    try { auto r = Precondensors_eqn_7_11__Q_condensor_heat_duty(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_11__U_v... ";
    try { auto r = Precondensors_eqn_7_11__U_v(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_11__V_c... ";
    try { auto r = Precondensors_eqn_7_11__V_c(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_11__del_T_LM... ";
    try { auto r = Precondensors_eqn_7_11__del_T_LM(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_12__A... ";
    try { auto r = Precondensors_eqn_7_12__A(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_12__Q_condensor_heat_duty... ";
    try { auto r = Precondensors_eqn_7_12__Q_condensor_heat_duty(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_12__U... ";
    try { auto r = Precondensors_eqn_7_12__U(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_12__del_T... ";
    try { auto r = Precondensors_eqn_7_12__del_T(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14a__A... ";
    try { auto r = Precondensors_eqn_7_14a__A(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14a__Q_condensor_heat_duty... ";
    try { auto r = Precondensors_eqn_7_14a__Q_condensor_heat_duty(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14a__U... ";
    try { auto r = Precondensors_eqn_7_14a__U(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14a__del_T_LM... ";
    try { auto r = Precondensors_eqn_7_14a__del_T_LM(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14b__A... ";
    try { auto r = Precondensors_eqn_7_14b__A(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14b__Q_condensor_heat_duty... ";
    try { auto r = Precondensors_eqn_7_14b__Q_condensor_heat_duty(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14b__U... ";
    try { auto r = Precondensors_eqn_7_14b__U(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14b__del_T_1... ";
    try { auto r = Precondensors_eqn_7_14b__del_T_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_14b__del_T_2... ";
    try { auto r = Precondensors_eqn_7_14b__del_T_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_15__U... ";
    try { auto r = Precondensors_eqn_7_15__U(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_15__sum_R... ";
    try { auto r = Precondensors_eqn_7_15__sum_R(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__D_0... ";
    try { auto r = Precondensors_eqn_7_16__D_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__D_LM... ";
    try { auto r = Precondensors_eqn_7_16__D_LM(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__D_i... ";
    try { auto r = Precondensors_eqn_7_16__D_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__R_f_0... ";
    try { auto r = Precondensors_eqn_7_16__R_f_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__R_fi... ";
    try { auto r = Precondensors_eqn_7_16__R_fi(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__U_0... ";
    try { auto r = Precondensors_eqn_7_16__U_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__h_0... ";
    try { auto r = Precondensors_eqn_7_16__h_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__h_i... ";
    try { auto r = Precondensors_eqn_7_16__h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__k_w... ";
    try { auto r = Precondensors_eqn_7_16__k_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_16__x_w... ";
    try { auto r = Precondensors_eqn_7_16__x_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_17__R_0... ";
    try { auto r = Precondensors_eqn_7_17__R_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_17__R_nc... ";
    try { auto r = Precondensors_eqn_7_17__R_nc(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_17__h_c... ";
    try { auto r = Precondensors_eqn_7_17__h_c(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__D_0... ";
    try { auto r = Precondensors_eqn_7_18__D_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__D_LM... ";
    try { auto r = Precondensors_eqn_7_18__D_LM(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__D_i... ";
    try { auto r = Precondensors_eqn_7_18__D_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__R_fi... ";
    try { auto r = Precondensors_eqn_7_18__R_fi(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__R_fo... ";
    try { auto r = Precondensors_eqn_7_18__R_fo(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__R_nc... ";
    try { auto r = Precondensors_eqn_7_18__R_nc(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__U_0... ";
    try { auto r = Precondensors_eqn_7_18__U_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__h_c... ";
    try { auto r = Precondensors_eqn_7_18__h_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__h_i... ";
    try { auto r = Precondensors_eqn_7_18__h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__k_w... ";
    try { auto r = Precondensors_eqn_7_18__k_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_18__x_w... ";
    try { auto r = Precondensors_eqn_7_18__x_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_2__P_i_0... ";
    try { auto r = Precondensors_eqn_7_2__P_i_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_2__p_i... ";
    try { auto r = Precondensors_eqn_7_2__p_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_2__x_i... ";
    try { auto r = Precondensors_eqn_7_2__x_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_3__P_i_0... ";
    try { auto r = Precondensors_eqn_7_3__P_i_0(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_3__epsilon_i... ";
    try { auto r = Precondensors_eqn_7_3__epsilon_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_3__p_i... ";
    try { auto r = Precondensors_eqn_7_3__p_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_3__x_i... ";
    try { auto r = Precondensors_eqn_7_3__x_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4a__P... ";
    try { auto r = Precondensors_eqn_7_4a__P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4a__p_c... ";
    try { auto r = Precondensors_eqn_7_4a__p_c(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4a__p_nc... ";
    try { auto r = Precondensors_eqn_7_4a__p_nc(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4aa__n_i... ";
    try { auto r = Precondensors_eqn_7_4aa__n_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4aa__n_nc... ";
    try { auto r = Precondensors_eqn_7_4aa__n_nc(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4aa__p_i... ";
    try { auto r = Precondensors_eqn_7_4aa__p_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4aa__p_nc... ";
    try { auto r = Precondensors_eqn_7_4aa__p_nc(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ab__P_c... ";
    try { auto r = Precondensors_eqn_7_4ab__P_c(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ab__p... ";
    try { auto r = Precondensors_eqn_7_4ab__p(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ab__p_i... ";
    try { auto r = Precondensors_eqn_7_4ab__p_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ab__p_nc... ";
    try { auto r = Precondensors_eqn_7_4ab__p_nc(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ac__P_c... ";
    try { auto r = Precondensors_eqn_7_4ac__P_c(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ac__n_i... ";
    try { auto r = Precondensors_eqn_7_4ac__n_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ac__n_nc... ";
    try { auto r = Precondensors_eqn_7_4ac__n_nc(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ac__p... ";
    try { auto r = Precondensors_eqn_7_4ac__p(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_4ac__p_i... ";
    try { auto r = Precondensors_eqn_7_4ac__p_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_5__N_i... ";
    try { auto r = Precondensors_eqn_7_5__N_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_5__N_nc... ";
    try { auto r = Precondensors_eqn_7_5__N_nc(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_5__P... ";
    try { auto r = Precondensors_eqn_7_5__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_5__P_c... ";
    try { auto r = Precondensors_eqn_7_5__P_c(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_5__p_i... ";
    try { auto r = Precondensors_eqn_7_5__p_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_6__M... ";
    try { auto r = Precondensors_eqn_7_6__M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_6__P... ";
    try { auto r = Precondensors_eqn_7_6__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_6__P_i_0... ";
    try { auto r = Precondensors_eqn_7_6__P_i_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_6__W_air... ";
    try { auto r = Precondensors_eqn_7_6__W_air(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_6__W_i... ";
    try { auto r = Precondensors_eqn_7_6__W_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_6__p_c... ";
    try { auto r = Precondensors_eqn_7_6__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_6__x_i... ";
    try { auto r = Precondensors_eqn_7_6__x_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__M... ";
    try { auto r = Precondensors_eqn_7_7__M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__P... ";
    try { auto r = Precondensors_eqn_7_7__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__P_i_0... ";
    try { auto r = Precondensors_eqn_7_7__P_i_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__W_air... ";
    try { auto r = Precondensors_eqn_7_7__W_air(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__W_i... ";
    try { auto r = Precondensors_eqn_7_7__W_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__epsilon_i... ";
    try { auto r = Precondensors_eqn_7_7__epsilon_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__p_c... ";
    try { auto r = Precondensors_eqn_7_7__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_7__x_i... ";
    try { auto r = Precondensors_eqn_7_7__x_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_8__L_c... ";
    try { auto r = Precondensors_eqn_7_8__L_c(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_8__Q_condensor_heat_duty... ";
    try { auto r = Precondensors_eqn_7_8__Q_condensor_heat_duty(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_8__c_p... ";
    try { auto r = Precondensors_eqn_7_8__c_p(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_8__del_T... ";
    try { auto r = Precondensors_eqn_7_8__del_T(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_9__L_c... ";
    try { auto r = Precondensors_eqn_7_9__L_c(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_9__Q_condensor_heat_duty... ";
    try { auto r = Precondensors_eqn_7_9__Q_condensor_heat_duty(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_9__c_p... ";
    try { auto r = Precondensors_eqn_7_9__c_p(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_9__del_T... ";
    try { auto r = Precondensors_eqn_7_9__del_T(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Precondensors_eqn_7_9__rho... ";
    try { auto r = Precondensors_eqn_7_9__rho(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_1__Abs_Pressure... ";
    try { auto r = PressMgmt_eqn_3_1__Abs_Pressure(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_1__BarometricPressure... ";
    try { auto r = PressMgmt_eqn_3_1__BarometricPressure(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_1__Vacuum... ";
    try { auto r = PressMgmt_eqn_3_1__Vacuum(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_11__A_C... ";
    try { auto r = PressMgmt_eqn_3_11__A_C(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_11__H_2... ";
    try { auto r = PressMgmt_eqn_3_11__H_2(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_11__P... ";
    try { auto r = PressMgmt_eqn_3_11__P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_11__V... ";
    try { auto r = PressMgmt_eqn_3_11__V(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_12__H_2... ";
    try { auto r = PressMgmt_eqn_3_12__H_2(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_12__KAPPA_1... ";
    try { auto r = PressMgmt_eqn_3_12__KAPPA_1(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_12__P... ";
    try { auto r = PressMgmt_eqn_3_12__P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_13__H_1... ";
    try { auto r = PressMgmt_eqn_3_13__H_1(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_13__H_2... ";
    try { auto r = PressMgmt_eqn_3_13__H_2(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_13__KAPPA_2... ";
    try { auto r = PressMgmt_eqn_3_13__KAPPA_2(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_13__P... ";
    try { auto r = PressMgmt_eqn_3_13__P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_15__V_PMIN... ";
    try { auto r = PressMgmt_eqn_3_15__V_PMIN(); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_16__V_div_V_P_MAX... ";
    try { auto r = PressMgmt_eqn_3_16__V_div_V_P_MAX(); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_17__P_MIN... ";
    try { auto r = PressMgmt_eqn_3_17__P_MIN(); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_2__G... ";
    try { auto r = PressMgmt_eqn_3_2__G(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_2__G_C... ";
    try { auto r = PressMgmt_eqn_3_2__G_C(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_2__H... ";
    try { auto r = PressMgmt_eqn_3_2__H(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_2__P... ";
    try { auto r = PressMgmt_eqn_3_2__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_2__rho... ";
    try { auto r = PressMgmt_eqn_3_2__rho(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_3__H_1... ";
    try { auto r = PressMgmt_eqn_3_3__H_1(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_3__H_2... ";
    try { auto r = PressMgmt_eqn_3_3__H_2(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_3__P... ";
    try { auto r = PressMgmt_eqn_3_3__P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_3__P_P... ";
    try { auto r = PressMgmt_eqn_3_3__P_P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_4__KAPPA... ";
    try { auto r = PressMgmt_eqn_3_4__KAPPA(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_4__P... ";
    try { auto r = PressMgmt_eqn_3_4__P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_4__V... ";
    try { auto r = PressMgmt_eqn_3_4__V(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_5__P... ";
    try { auto r = PressMgmt_eqn_3_5__P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_5__P_P... ";
    try { auto r = PressMgmt_eqn_3_5__P_P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_5__V... ";
    try { auto r = PressMgmt_eqn_3_5__V(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_5__V_P... ";
    try { auto r = PressMgmt_eqn_3_5__V_P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_6__H_1... ";
    try { auto r = PressMgmt_eqn_3_6__H_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_6__H_2... ";
    try { auto r = PressMgmt_eqn_3_6__H_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_6__P... ";
    try { auto r = PressMgmt_eqn_3_6__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_6__V... ";
    try { auto r = PressMgmt_eqn_3_6__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_6__V_P... ";
    try { auto r = PressMgmt_eqn_3_6__V_P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_8__A_C... ";
    try { auto r = PressMgmt_eqn_3_8__A_C(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_8__H_2... ";
    try { auto r = PressMgmt_eqn_3_8__H_2(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_8__V_P... ";
    try { auto r = PressMgmt_eqn_3_8__V_P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_9__A_C... ";
    try { auto r = PressMgmt_eqn_3_9__A_C(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_9__H_1... ";
    try { auto r = PressMgmt_eqn_3_9__H_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_9__H_2... ";
    try { auto r = PressMgmt_eqn_3_9__H_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_9__P... ";
    try { auto r = PressMgmt_eqn_3_9__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing PressMgmt_eqn_3_9__V... ";
    try { auto r = PressMgmt_eqn_3_9__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_1__K_i... ";
    try { auto r = ProcessApp1_eqn_5_1__K_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_1__x_i... ";
    try { auto r = ProcessApp1_eqn_5_1__x_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_1__y_i... ";
    try { auto r = ProcessApp1_eqn_5_1__y_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10a__D... ";
    try { auto r = ProcessApp1_eqn_5_10a__D(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10a__L_0... ";
    try { auto r = ProcessApp1_eqn_5_10a__L_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10a__V_1... ";
    try { auto r = ProcessApp1_eqn_5_10a__V_1(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10b__L_0... ";
    try { auto r = ProcessApp1_eqn_5_10b__L_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10b__R... ";
    try { auto r = ProcessApp1_eqn_5_10b__R(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10b__V_1... ";
    try { auto r = ProcessApp1_eqn_5_10b__V_1(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10c__D... ";
    try { auto r = ProcessApp1_eqn_5_10c__D(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10c__L_0... ";
    try { auto r = ProcessApp1_eqn_5_10c__L_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_10c__R... ";
    try { auto r = ProcessApp1_eqn_5_10c__R(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_11__B... ";
    try { auto r = ProcessApp1_eqn_5_11__B(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_11__L_N... ";
    try { auto r = ProcessApp1_eqn_5_11__L_N(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_11__V_0... ";
    try { auto r = ProcessApp1_eqn_5_11__V_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_12__Eff... ";
    try { auto r = ProcessApp1_eqn_5_12__Eff(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_12__N_ES... ";
    try { auto r = ProcessApp1_eqn_5_12__N_ES(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_12__N_t... ";
    try { auto r = ProcessApp1_eqn_5_12__N_t(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_12__T... ";
    try { auto r = ProcessApp1_eqn_5_12__T(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_13__HETP... ";
    try { auto r = ProcessApp1_eqn_5_13__HETP(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_13__H_p... ";
    try { auto r = ProcessApp1_eqn_5_13__H_p(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_13__N_ES... ";
    try { auto r = ProcessApp1_eqn_5_13__N_ES(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_14__M... ";
    try { auto r = ProcessApp1_eqn_5_14__M(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_14__P_0... ";
    try { auto r = ProcessApp1_eqn_5_14__P_0(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_14__T... ";
    try { auto r = ProcessApp1_eqn_5_14__T(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_14__W_E... ";
    try { auto r = ProcessApp1_eqn_5_14__W_E(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_15__M_1... ";
    try { auto r = ProcessApp1_eqn_5_15__M_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_15__M_2... ";
    try { auto r = ProcessApp1_eqn_5_15__M_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_15__P_0_1... ";
    try { auto r = ProcessApp1_eqn_5_15__P_0_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_15__P_0_2... ";
    try { auto r = ProcessApp1_eqn_5_15__P_0_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_15__a_M_12... ";
    try { auto r = ProcessApp1_eqn_5_15__a_M_12(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_16__H_i... ";
    try { auto r = ProcessApp1_eqn_5_16__H_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_16__p_i... ";
    try { auto r = ProcessApp1_eqn_5_16__p_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_16__x_i... ";
    try { auto r = ProcessApp1_eqn_5_16__x_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_17__H_2_1... ";
    try { auto r = ProcessApp1_eqn_5_17__H_2_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_17__H_2_3... ";
    try { auto r = ProcessApp1_eqn_5_17__H_2_3(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_17__H_2_mi... ";
    try { auto r = ProcessApp1_eqn_5_17__H_2_mi(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_17__x_1... ";
    try { auto r = ProcessApp1_eqn_5_17__x_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_17__x_3... ";
    try { auto r = ProcessApp1_eqn_5_17__x_3(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2a__K_1... ";
    try { auto r = ProcessApp1_eqn_5_2a__K_1(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2a__K_2... ";
    try { auto r = ProcessApp1_eqn_5_2a__K_2(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2a__alpha_1_2... ";
    try { auto r = ProcessApp1_eqn_5_2a__alpha_1_2(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2b__K_1... ";
    try { auto r = ProcessApp1_eqn_5_2b__K_1(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2b__K_2... ";
    try { auto r = ProcessApp1_eqn_5_2b__K_2(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2b__x_1... ";
    try { auto r = ProcessApp1_eqn_5_2b__x_1(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2b__x_2... ";
    try { auto r = ProcessApp1_eqn_5_2b__x_2(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2b__y_1... ";
    try { auto r = ProcessApp1_eqn_5_2b__y_1(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_2b__y_2... ";
    try { auto r = ProcessApp1_eqn_5_2b__y_2(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_3__P_0_i... ";
    try { auto r = ProcessApp1_eqn_5_3__P_0_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_3__p_i... ";
    try { auto r = ProcessApp1_eqn_5_3__p_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_3__x_i... ";
    try { auto r = ProcessApp1_eqn_5_3__x_i(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_4__P... ";
    try { auto r = ProcessApp1_eqn_5_4__P(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_4__P_0_i... ";
    try { auto r = ProcessApp1_eqn_5_4__P_0_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_4__x_i... ";
    try { auto r = ProcessApp1_eqn_5_4__x_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_4__y_i... ";
    try { auto r = ProcessApp1_eqn_5_4__y_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_5__P_0_1... ";
    try { auto r = ProcessApp1_eqn_5_5__P_0_1(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_5__P_0_2... ";
    try { auto r = ProcessApp1_eqn_5_5__P_0_2(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_5__alpha_12... ";
    try { auto r = ProcessApp1_eqn_5_5__alpha_12(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_6__P_0_i... ";
    try { auto r = ProcessApp1_eqn_5_6__P_0_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_6__gamma_i... ";
    try { auto r = ProcessApp1_eqn_5_6__gamma_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_6__p_i... ";
    try { auto r = ProcessApp1_eqn_5_6__p_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_6__x_i... ";
    try { auto r = ProcessApp1_eqn_5_6__x_i(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_7__P... ";
    try { auto r = ProcessApp1_eqn_5_7__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_7__P_0_i... ";
    try { auto r = ProcessApp1_eqn_5_7__P_0_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_7__gamma_i... ";
    try { auto r = ProcessApp1_eqn_5_7__gamma_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_7__x_i... ";
    try { auto r = ProcessApp1_eqn_5_7__x_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_7__y_i... ";
    try { auto r = ProcessApp1_eqn_5_7__y_i(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_8__P_0_1... ";
    try { auto r = ProcessApp1_eqn_5_8__P_0_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_8__P_0_2... ";
    try { auto r = ProcessApp1_eqn_5_8__P_0_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_8__alpha_12... ";
    try { auto r = ProcessApp1_eqn_5_8__alpha_12(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_8__gamma_1... ";
    try { auto r = ProcessApp1_eqn_5_8__gamma_1(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_8__gamma_2... ";
    try { auto r = ProcessApp1_eqn_5_8__gamma_2(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_9__D... ";
    try { auto r = ProcessApp1_eqn_5_9__D(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_9__L_0... ";
    try { auto r = ProcessApp1_eqn_5_9__L_0(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp1_eqn_5_9__V_1... ";
    try { auto r = ProcessApp1_eqn_5_9__V_1(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__T_1... ";
    try { auto r = ProcessApp2_eqn_6_1__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__T_2... ";
    try { auto r = ProcessApp2_eqn_6_1__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__T_R... ";
    try { auto r = ProcessApp2_eqn_6_1__T_R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__c_p... ";
    try { auto r = ProcessApp2_eqn_6_1__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__del_h_v... ";
    try { auto r = ProcessApp2_eqn_6_1__del_h_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__w_1... ";
    try { auto r = ProcessApp2_eqn_6_1__w_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__w_2... ";
    try { auto r = ProcessApp2_eqn_6_1__w_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_1__w_v... ";
    try { auto r = ProcessApp2_eqn_6_1__w_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_10__A... ";
    try { auto r = ProcessApp2_eqn_6_10__A(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_10__dV_dt... ";
    try { auto r = ProcessApp2_eqn_6_10__dV_dt(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_10__delta_P... ";
    try { auto r = ProcessApp2_eqn_6_10__delta_P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_10__mu... ";
    try { auto r = ProcessApp2_eqn_6_10__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_10__r_c... ";
    try { auto r = ProcessApp2_eqn_6_10__r_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_10__s... ";
    try { auto r = ProcessApp2_eqn_6_10__s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_10__tau... ";
    try { auto r = ProcessApp2_eqn_6_10__tau(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_11a__A_d... ";
    try { auto r = ProcessApp2_eqn_6_11a__A_d(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_11a__delta_T... ";
    try { auto r = ProcessApp2_eqn_6_11a__delta_T(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_11a__delta_h_i... ";
    try { auto r = ProcessApp2_eqn_6_11a__delta_h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_11a__delta_m... ";
    try { auto r = ProcessApp2_eqn_6_11a__delta_m(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_11a__h_d... ";
    try { auto r = ProcessApp2_eqn_6_11a__h_d(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_11a__m_b... ";
    try { auto r = ProcessApp2_eqn_6_11a__m_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_11a__t_R... ";
    try { auto r = ProcessApp2_eqn_6_11a__t_R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_2__Q_v... ";
    try { auto r = ProcessApp2_eqn_6_2__Q_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_2__T_1... ";
    try { auto r = ProcessApp2_eqn_6_2__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_2__T_2... ";
    try { auto r = ProcessApp2_eqn_6_2__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_2__T_R... ";
    try { auto r = ProcessApp2_eqn_6_2__T_R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_2__c_p... ";
    try { auto r = ProcessApp2_eqn_6_2__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_2__w_1... ";
    try { auto r = ProcessApp2_eqn_6_2__w_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_2__w_2... ";
    try { auto r = ProcessApp2_eqn_6_2__w_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_4__Q_v... ";
    try { auto r = ProcessApp2_eqn_6_4__Q_v(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_4__delta_h_v... ";
    try { auto r = ProcessApp2_eqn_6_4__delta_h_v(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_4__w_v... ";
    try { auto r = ProcessApp2_eqn_6_4__w_v(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__C_1... ";
    try { auto r = ProcessApp2_eqn_6_7__C_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__C_2... ";
    try { auto r = ProcessApp2_eqn_6_7__C_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__T_1... ";
    try { auto r = ProcessApp2_eqn_6_7__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__T_2... ";
    try { auto r = ProcessApp2_eqn_6_7__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__c_p... ";
    try { auto r = ProcessApp2_eqn_6_7__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__delta_h_c... ";
    try { auto r = ProcessApp2_eqn_6_7__delta_h_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__delta_h_v... ";
    try { auto r = ProcessApp2_eqn_6_7__delta_h_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__m_b... ";
    try { auto r = ProcessApp2_eqn_6_7__m_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_7__m_v... ";
    try { auto r = ProcessApp2_eqn_6_7__m_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__C_1... ";
    try { auto r = ProcessApp2_eqn_6_8__C_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__C_2... ";
    try { auto r = ProcessApp2_eqn_6_8__C_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__T_1... ";
    try { auto r = ProcessApp2_eqn_6_8__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__T_2... ";
    try { auto r = ProcessApp2_eqn_6_8__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__c_p... ";
    try { auto r = ProcessApp2_eqn_6_8__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__delta_h_c... ";
    try { auto r = ProcessApp2_eqn_6_8__delta_h_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__delta_h_v... ";
    try { auto r = ProcessApp2_eqn_6_8__delta_h_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__delta_t... ";
    try { auto r = ProcessApp2_eqn_6_8__delta_t(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__m_b... ";
    try { auto r = ProcessApp2_eqn_6_8__m_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_8__w_v... ";
    try { auto r = ProcessApp2_eqn_6_8__w_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_9__A... ";
    try { auto r = ProcessApp2_eqn_6_9__A(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_9__dV_dt... ";
    try { auto r = ProcessApp2_eqn_6_9__dV_dt(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_9__delta_P... ";
    try { auto r = ProcessApp2_eqn_6_9__delta_P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_9__m... ";
    try { auto r = ProcessApp2_eqn_6_9__m(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_9__mu... ";
    try { auto r = ProcessApp2_eqn_6_9__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_9__r... ";
    try { auto r = ProcessApp2_eqn_6_9__r(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing ProcessApp2_eqn_6_9__r_M... ";
    try { auto r = ProcessApp2_eqn_6_9__r_M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_1__PS... ";
    try { auto r = RotaryPistonVane_eqn_11_1__PS(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_1__Q_0... ";
    try { auto r = RotaryPistonVane_eqn_11_1__Q_0(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_1__Q_external_gas_throughput... ";
    try { auto r = RotaryPistonVane_eqn_11_1__Q_external_gas_throughput(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_1__V... ";
    try { auto r = RotaryPistonVane_eqn_11_1__V(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_1__dP... ";
    try { auto r = RotaryPistonVane_eqn_11_1__dP(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_1__dT... ";
    try { auto r = RotaryPistonVane_eqn_11_1__dT(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__Q... ";
    try { auto r = RotaryPistonVane_eqn_11_2__Q(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__Q_0... ";
    try { auto r = RotaryPistonVane_eqn_11_2__Q_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__Q_external_gas_throughput... ";
    try { auto r = RotaryPistonVane_eqn_11_2__Q_external_gas_throughput(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__SP_1... ";
    try { auto r = RotaryPistonVane_eqn_11_2__SP_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__SP_2... ";
    try { auto r = RotaryPistonVane_eqn_11_2__SP_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__S_vol_pump_speed... ";
    try { auto r = RotaryPistonVane_eqn_11_2__S_vol_pump_speed(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__V... ";
    try { auto r = RotaryPistonVane_eqn_11_2__V(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_2__t... ";
    try { auto r = RotaryPistonVane_eqn_11_2__t(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_3__F_s... ";
    try { auto r = RotaryPistonVane_eqn_11_3__F_s(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_3__t... ";
    try { auto r = RotaryPistonVane_eqn_11_3__t(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_3__t_c... ";
    try { auto r = RotaryPistonVane_eqn_11_3__t_c(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_4__p_g... ";
    try { auto r = RotaryPistonVane_eqn_11_4__p_g(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_4__p_s... ";
    try { auto r = RotaryPistonVane_eqn_11_4__p_s(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_4__p_v... ";
    try { auto r = RotaryPistonVane_eqn_11_4__p_v(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_5__P_0_v... ";
    try { auto r = RotaryPistonVane_eqn_11_5__P_0_v(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_5__P_D... ";
    try { auto r = RotaryPistonVane_eqn_11_5__P_D(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_5__p_g... ";
    try { auto r = RotaryPistonVane_eqn_11_5__p_g(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_5__p_v_max... ";
    try { auto r = RotaryPistonVane_eqn_11_5__p_v_max(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__P_0_V... ";
    try { auto r = RotaryPistonVane_eqn_11_6__P_0_V(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__P_D... ";
    try { auto r = RotaryPistonVane_eqn_11_6__P_D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__P_v_0... ";
    try { auto r = RotaryPistonVane_eqn_11_6__P_v_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__S_B... ";
    try { auto r = RotaryPistonVane_eqn_11_6__S_B(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__S_D... ";
    try { auto r = RotaryPistonVane_eqn_11_6__S_D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__p_b... ";
    try { auto r = RotaryPistonVane_eqn_11_6__p_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__p_g... ";
    try { auto r = RotaryPistonVane_eqn_11_6__p_g(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing RotaryPistonVane_eqn_11_6__p_v_max... ";
    try { auto r = RotaryPistonVane_eqn_11_6__p_v_max(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_1__NC... ";
    try { auto r = SelectingPump_eqn_8_1__NC(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_1__NS... ";
    try { auto r = SelectingPump_eqn_8_1__NS(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_1__SCON... ";
    try { auto r = SelectingPump_eqn_8_1__SCON(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_1__installation_cost... ";
    try { auto r = SelectingPump_eqn_8_1__installation_cost(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_2__hp... ";
    try { auto r = SelectingPump_eqn_8_2__hp(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_2__installed_costs... ";
    try { auto r = SelectingPump_eqn_8_2__installed_costs(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_3__hp... ";
    try { auto r = SelectingPump_eqn_8_3__hp(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_3__installed_costs... ";
    try { auto r = SelectingPump_eqn_8_3__installed_costs(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_4__hp... ";
    try { auto r = SelectingPump_eqn_8_4__hp(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_4__installed_costs... ";
    try { auto r = SelectingPump_eqn_8_4__installed_costs(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_5__Eff... ";
    try { auto r = SelectingPump_eqn_8_5__Eff(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_5__actual_brake_horsepower... ";
    try { auto r = SelectingPump_eqn_8_5__actual_brake_horsepower(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_5__theoretical_adiabatic_horsepower... ";
    try { auto r = SelectingPump_eqn_8_5__theoretical_adiabatic_horsepower(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__M... ";
    try { auto r = SelectingPump_eqn_8_6__M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__P_1... ";
    try { auto r = SelectingPump_eqn_8_6__P_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__P_2... ";
    try { auto r = SelectingPump_eqn_8_6__P_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__R... ";
    try { auto r = SelectingPump_eqn_8_6__R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__T... ";
    try { auto r = SelectingPump_eqn_8_6__T(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__adiabatic_hp... ";
    try { auto r = SelectingPump_eqn_8_6__adiabatic_hp(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__k... ";
    try { auto r = SelectingPump_eqn_8_6__k(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_6__w... ";
    try { auto r = SelectingPump_eqn_8_6__w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_7__P_1... ";
    try { auto r = SelectingPump_eqn_8_7__P_1(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_7__P_2... ";
    try { auto r = SelectingPump_eqn_8_7__P_2(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_7__adiabatic_hp... ";
    try { auto r = SelectingPump_eqn_8_7__adiabatic_hp(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_7__w... ";
    try { auto r = SelectingPump_eqn_8_7__w(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_8__P_1... ";
    try { auto r = SelectingPump_eqn_8_8__P_1(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_8__P_2... ";
    try { auto r = SelectingPump_eqn_8_8__P_2(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_8__adiabatic_power_watts... ";
    try { auto r = SelectingPump_eqn_8_8__adiabatic_power_watts(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_8__f... ";
    try { auto r = SelectingPump_eqn_8_8__f(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_9__E_j... ";
    try { auto r = SelectingPump_eqn_8_9__E_j(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_9__E_m... ";
    try { auto r = SelectingPump_eqn_8_9__E_m(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_9__e... ";
    try { auto r = SelectingPump_eqn_8_9__e(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_9__r... ";
    try { auto r = SelectingPump_eqn_8_9__r(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SelectingPump_eqn_8_9__s... ";
    try { auto r = SelectingPump_eqn_8_9__s(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_1__A... ";
    try { auto r = SteamJetInjectors_eqn_9_1__A(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_1__rho_s... ";
    try { auto r = SteamJetInjectors_eqn_9_1__rho_s(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_1__v... ";
    try { auto r = SteamJetInjectors_eqn_9_1__v(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_1__w_s... ";
    try { auto r = SteamJetInjectors_eqn_9_1__w_s(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_2__P_m... ";
    try { auto r = SteamJetInjectors_eqn_9_2__P_m(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_2__d_n... ";
    try { auto r = SteamJetInjectors_eqn_9_2__d_n(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_2__rho_s... ";
    try { auto r = SteamJetInjectors_eqn_9_2__rho_s(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_2__w_s... ";
    try { auto r = SteamJetInjectors_eqn_9_2__w_s(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_3__P_s... ";
    try { auto r = SteamJetInjectors_eqn_9_3__P_s(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_3__V... ";
    try { auto r = SteamJetInjectors_eqn_9_3__V(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_3__t_e... ";
    try { auto r = SteamJetInjectors_eqn_9_3__t_e(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_3__w_j... ";
    try { auto r = SteamJetInjectors_eqn_9_3__w_j(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_4__AEL... ";
    try { auto r = SteamJetInjectors_eqn_9_4__AEL(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_4__SC... ";
    try { auto r = SteamJetInjectors_eqn_9_4__SC(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_4__r... ";
    try { auto r = SteamJetInjectors_eqn_9_4__r(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_4__w_s... ";
    try { auto r = SteamJetInjectors_eqn_9_4__w_s(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_5__V... ";
    try { auto r = SteamJetInjectors_eqn_9_5__V(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_5__r_h... ";
    try { auto r = SteamJetInjectors_eqn_9_5__r_h(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_5__t_h... ";
    try { auto r = SteamJetInjectors_eqn_9_5__t_h(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing SteamJetInjectors_eqn_9_5__w_h... ";
    try { auto r = SteamJetInjectors_eqn_9_5__w_h(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_10__P_1... ";
    try { auto r = VacuumTheory_eqn_1_10__P_1(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_10__P_2... ";
    try { auto r = VacuumTheory_eqn_1_10__P_2(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_10__T_1... ";
    try { auto r = VacuumTheory_eqn_1_10__T_1(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_10__T_2... ";
    try { auto r = VacuumTheory_eqn_1_10__T_2(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_10__V_1... ";
    try { auto r = VacuumTheory_eqn_1_10__V_1(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_10__V_2... ";
    try { auto r = VacuumTheory_eqn_1_10__V_2(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_11__M... ";
    try { auto r = VacuumTheory_eqn_1_11__M(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_11__P... ";
    try { auto r = VacuumTheory_eqn_1_11__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_11__T... ";
    try { auto r = VacuumTheory_eqn_1_11__T(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_11__W... ";
    try { auto r = VacuumTheory_eqn_1_11__W(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_11__q... ";
    try { auto r = VacuumTheory_eqn_1_11__q(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_12__Total_P... ";
    try { auto r = VacuumTheory_eqn_1_12__Total_P(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_12__sum_partial_pressures... ";
    try { auto r = VacuumTheory_eqn_1_12__sum_partial_pressures(1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_13a__n... ";
    try { auto r = VacuumTheory_eqn_1_13a__n(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_13a__n_a... ";
    try { auto r = VacuumTheory_eqn_1_13a__n_a(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_13a__y_a... ";
    try { auto r = VacuumTheory_eqn_1_13a__y_a(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_13b__P... ";
    try { auto r = VacuumTheory_eqn_1_13b__P(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_13b__p_a... ";
    try { auto r = VacuumTheory_eqn_1_13b__p_a(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_13b__y_a... ";
    try { auto r = VacuumTheory_eqn_1_13b__y_a(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_3__T... ";
    try { auto r = VacuumTheory_eqn_1_3__T(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_3__k... ";
    try { auto r = VacuumTheory_eqn_1_3__k(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_3__m... ";
    try { auto r = VacuumTheory_eqn_1_3__m(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_3__v... ";
    try { auto r = VacuumTheory_eqn_1_3__v(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_7__R... ";
    try { auto r = VacuumTheory_eqn_1_7__R(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_7__T... ";
    try { auto r = VacuumTheory_eqn_1_7__T(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_7__V... ";
    try { auto r = VacuumTheory_eqn_1_7__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_7__n... ";
    try { auto r = VacuumTheory_eqn_1_7__n(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_7__p... ";
    try { auto r = VacuumTheory_eqn_1_7__p(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_8__M... ";
    try { auto r = VacuumTheory_eqn_1_8__M(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_8__P... ";
    try { auto r = VacuumTheory_eqn_1_8__P(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_8__R... ";
    try { auto r = VacuumTheory_eqn_1_8__R(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_8__T... ";
    try { auto r = VacuumTheory_eqn_1_8__T(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_8__V... ";
    try { auto r = VacuumTheory_eqn_1_8__V(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_8__m... ";
    try { auto r = VacuumTheory_eqn_1_8__m(1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_9__M... ";
    try { auto r = VacuumTheory_eqn_1_9__M(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_9__P... ";
    try { auto r = VacuumTheory_eqn_1_9__P(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_9__R... ";
    try { auto r = VacuumTheory_eqn_1_9__R(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_9__T... ";
    try { auto r = VacuumTheory_eqn_1_9__T(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing VacuumTheory_eqn_1_9__rho... ";
    try { auto r = VacuumTheory_eqn_1_9__rho(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }

    std::cout << "\n" << pass << " passed, " << fail << " failed." << std::endl;
    return fail > 0 ? 1 : 0;
}