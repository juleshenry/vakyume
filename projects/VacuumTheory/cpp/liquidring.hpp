#ifndef VAKYUME_LIQUIDRING_HPP
#define VAKYUME_LIQUIDRING_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

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

#endif // VAKYUME_LIQUIDRING_HPP
