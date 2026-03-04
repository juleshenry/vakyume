#ifndef VAKYUME_SELECTINGPUMP_HPP
#define VAKYUME_SELECTINGPUMP_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

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

#endif // VAKYUME_SELECTINGPUMP_HPP
