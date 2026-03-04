#ifndef VAKYUME_PRECONDENSORS_HPP
#define VAKYUME_PRECONDENSORS_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

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
    throw std::runtime_error("Precondensors_eqn_7_14b__del_T_1: requires numerical solver (not transpilable)");
}

std::vector<double> Precondensors_eqn_7_14b__del_T_2(double A, double Q_condensor_heat_duty, double U, double del_T_1) {
    throw std::runtime_error("Precondensors_eqn_7_14b__del_T_2: requires numerical solver (not transpilable)");
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

#endif // VAKYUME_PRECONDENSORS_HPP
