#include <iostream>
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

int main() {
    std::cout << "Running test suite..." << std::endl;
    int pass = 0, fail = 0;
    std::cout << "Testing Kinematics_eqn_1_1__a... ";
    try { auto r = Kinematics_eqn_1_1__a(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_1__t... ";
    try { auto r = Kinematics_eqn_1_1__t(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_1__v... ";
    try { auto r = Kinematics_eqn_1_1__v(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_1__v0... ";
    try { auto r = Kinematics_eqn_1_1__v0(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }

    std::cout << "\n" << pass << " passed, " << fail << " failed." << std::endl;
    return fail > 0 ? 1 : 0;
}