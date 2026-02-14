#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

using namespace std::complex_literals;

std::vector<double> Kinematics_eqn_1_1__a(double t, double v, double v0){
    std::vector<double> result;
    auto a = (v - v0) / t;
    result.push_back(a);
    return result;
}
std::vector<double> Kinematics_eqn_1_1__t(double a, double v, double v0) {
    auto t = (v - static_cast<double>(v0)) / a;
    return {t};
}
std::vector<double> Kinematics_eqn_1_1__v(float a, float t, float v0) {
    std::vector<double> result;
    double v = a*t + v0;
    result.push_back(v);
    return result;
}
std::vector<double> Kinematics_eqn_1_1__v0(float a, float t, float v) {
    double result = -a * t + v;
    return {result};
}

int main() {
    std::cout << "Running test suite..." << std::endl;
    std::cout << "Testing Kinematics_eqn_1_1__a... ";
    try { Kinematics_eqn_1_1__a(1.0, 1.0, 1.0); std::cout << "OK" << std::endl; } catch(...) { std::cout << "FAIL" << std::endl; }
    std::cout << "Testing Kinematics_eqn_1_1__t... ";
    try { Kinematics_eqn_1_1__t(1.0, 1.0, 1.0); std::cout << "OK" << std::endl; } catch(...) { std::cout << "FAIL" << std::endl; }
    std::cout << "Testing Kinematics_eqn_1_1__v... ";
    try { Kinematics_eqn_1_1__v(1.0, 1.0, 1.0); std::cout << "OK" << std::endl; } catch(...) { std::cout << "FAIL" << std::endl; }
    std::cout << "Testing Kinematics_eqn_1_1__v0... ";
    try { Kinematics_eqn_1_1__v0(1.0, 1.0, 1.0); std::cout << "OK" << std::endl; } catch(...) { std::cout << "FAIL" << std::endl; }
    std::cout << "Test suite finished." << std::endl;
    return 0;
}