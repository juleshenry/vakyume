#ifndef VAKYUME_KINEMATICS_HPP
#define VAKYUME_KINEMATICS_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

std::vector<double> Kinematics_eqn_1_1__a(double t, double v, double v0) {
    std::vector<double> result;
    double a = ((v - v0) / t);
    result.push_back(a);
    return result;
}

std::vector<double> Kinematics_eqn_1_1__t(double a, double v, double v0) {
    std::vector<double> result;
    double t = ((v - v0) / a);
    result.push_back(t);
    return result;
}

std::vector<double> Kinematics_eqn_1_1__v(double a, double t, double v0) {
    std::vector<double> result;
    double v = ((a * t) + v0);
    result.push_back(v);
    return result;
}

std::vector<double> Kinematics_eqn_1_1__v0(double a, double t, double v) {
    std::vector<double> result;
    double v0 = (((-a) * t) + v);
    result.push_back(v0);
    return result;
}

std::vector<double> Kinematics_eqn_1_2__a(double t, double v0, double x, double x0) {
    std::vector<double> result;
    double a = ((2.0 * ((((-t) * v0) + x) - x0)) / std::pow(t, 2.0));
    result.push_back(a);
    return result;
}

std::vector<double> Kinematics_eqn_1_2__t(double a, double v0, double x, double x0) {
    std::vector<double> result;
    double t = (((-v0) - std::sqrt(((((2.0 * a) * x) - ((2.0 * a) * x0)) + std::pow(v0, 2.0)))) / a);
    result.push_back(t);
    t = (((-v0) + std::sqrt(((((2.0 * a) * x) - ((2.0 * a) * x0)) + std::pow(v0, 2.0)))) / a);
    result.push_back(t);
    return result;
}

std::vector<double> Kinematics_eqn_1_2__v0(double a, double t, double x, double x0) {
    std::vector<double> result;
    double v0 = ((((((-0.5) * a) * std::pow(t, 2.0)) + x) - x0) / t);
    result.push_back(v0);
    return result;
}

std::vector<double> Kinematics_eqn_1_2__x(double a, double t, double v0, double x0) {
    std::vector<double> result;
    double x = ((((0.5 * a) * std::pow(t, 2.0)) + (t * v0)) + x0);
    result.push_back(x);
    return result;
}

std::vector<double> Kinematics_eqn_1_2__x0(double a, double t, double v0, double x) {
    std::vector<double> result;
    double x0 = (((((-0.5) * a) * std::pow(t, 2.0)) - (t * v0)) + x);
    result.push_back(x0);
    return result;
}

std::vector<double> Kinematics_eqn_1_3__a(double dx, double v, double v0) {
    std::vector<double> result;
    double a = ((std::pow(v, 2.0) - std::pow(v0, 2.0)) / (2.0 * dx));
    result.push_back(a);
    return result;
}

std::vector<double> Kinematics_eqn_1_3__dx(double a, double v, double v0) {
    std::vector<double> result;
    double dx = ((std::pow(v, 2.0) - std::pow(v0, 2.0)) / (2.0 * a));
    result.push_back(dx);
    return result;
}

std::vector<double> Kinematics_eqn_1_3__v(double a, double dx, double v0) {
    std::vector<double> result;
    double v = (-std::sqrt((((2.0 * a) * dx) + std::pow(v0, 2.0))));
    result.push_back(v);
    v = std::sqrt((((2.0 * a) * dx) + std::pow(v0, 2.0)));
    result.push_back(v);
    return result;
}

std::vector<double> Kinematics_eqn_1_3__v0(double a, double dx, double v) {
    std::vector<double> result;
    double v0 = (-std::sqrt(((((-2.0) * a) * dx) + std::pow(v, 2.0))));
    result.push_back(v0);
    v0 = std::sqrt(((((-2.0) * a) * dx) + std::pow(v, 2.0)));
    result.push_back(v0);
    return result;
}

#endif // VAKYUME_KINEMATICS_HPP
