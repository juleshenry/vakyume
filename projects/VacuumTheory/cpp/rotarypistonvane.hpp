#ifndef VAKYUME_ROTARYPISTONVANE_HPP
#define VAKYUME_ROTARYPISTONVANE_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

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

#endif // VAKYUME_ROTARYPISTONVANE_HPP
