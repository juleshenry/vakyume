#ifndef VAKYUME_ROTARY_HPP
#define VAKYUME_ROTARY_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

std::vector<double> Rotary_eqn_11_2__Q(double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double Q = ((((Q_0 + Q_external_gas_throughput) - SP_1) + (((-Q_0) + SP_2) * std::exp(((S_vol_pump_speed * t) / V)))) * std::exp((((-S_vol_pump_speed) * t) / V)));
    result.push_back(Q);
    return result;
}

std::vector<double> Rotary_eqn_11_2__Q_0(double Q, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double Q_0 = ((-((((Q * std::exp(((S_vol_pump_speed * t) / V))) - Q_external_gas_throughput) + SP_1) - (SP_2 * std::exp(((S_vol_pump_speed * t) / V))))) / (std::exp(((S_vol_pump_speed * t) / V)) - 1.0));
    result.push_back(Q_0);
    return result;
}

std::vector<double> Rotary_eqn_11_2__Q_external_gas_throughput(double Q, double Q_0, double SP_1, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double Q_external_gas_throughput = (((-Q_0) + SP_1) + (((Q + Q_0) - SP_2) * std::exp(((S_vol_pump_speed * t) / V))));
    result.push_back(Q_external_gas_throughput);
    return result;
}

std::vector<double> Rotary_eqn_11_2__SP_1(double Q, double Q_0, double Q_external_gas_throughput, double SP_2, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double SP_1 = ((Q_0 + Q_external_gas_throughput) + ((((-Q) - Q_0) + SP_2) * std::exp(((S_vol_pump_speed * t) / V))));
    result.push_back(SP_1);
    return result;
}

std::vector<double> Rotary_eqn_11_2__SP_2(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double S_vol_pump_speed, double V, double t) {
    std::vector<double> result;
    double SP_2 = (((((-Q_0) - Q_external_gas_throughput) + SP_1) + ((Q + Q_0) * std::exp(((S_vol_pump_speed * t) / V)))) * std::exp((((-S_vol_pump_speed) * t) / V)));
    result.push_back(SP_2);
    return result;
}

std::vector<double> Rotary_eqn_11_2__S_vol_pump_speed(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double V, double t) {
    std::vector<double> result;
    double S_vol_pump_speed = ((V * std::log((((Q_0 + Q_external_gas_throughput) - SP_1) / ((Q + Q_0) - SP_2)))) / t);
    result.push_back(S_vol_pump_speed);
    return result;
}

std::vector<double> Rotary_eqn_11_2__V(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double t) {
    std::vector<double> result;
    double V = ((S_vol_pump_speed * t) / std::log((((Q_0 + Q_external_gas_throughput) - SP_1) / ((Q + Q_0) - SP_2))));
    result.push_back(V);
    return result;
}

std::vector<double> Rotary_eqn_11_2__t(double Q, double Q_0, double Q_external_gas_throughput, double SP_1, double SP_2, double S_vol_pump_speed, double V) {
    std::vector<double> result;
    double t = ((V * std::log((((Q_0 + Q_external_gas_throughput) - SP_1) / ((Q + Q_0) - SP_2)))) / S_vol_pump_speed);
    result.push_back(t);
    return result;
}

#endif // VAKYUME_ROTARY_HPP
