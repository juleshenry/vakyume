#ifndef VAKYUME_STEAMJETINJECTORS_HPP
#define VAKYUME_STEAMJETINJECTORS_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

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

#endif // VAKYUME_STEAMJETINJECTORS_HPP
