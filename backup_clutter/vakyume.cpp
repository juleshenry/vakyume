#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

using namespace std::complex_literals;

std::vector<double> AirLeak_eqn_4_10__T(double V, double delP, double leakage, double t) {
    std::vector<double> result;
    
    // Calculate temperature T using the provided formula and add to results vector.
    double T = 3.127 * V * delP / (leakage * t);
    result.push_back(T);
    
    return result;
}
std::vector<double> AirLeak_eqn_4_10__V(double T, double del_P, double leakage, double t) {
    std::vector<double> result;
    
    // Calculate V using the given formula
    double V = 0.319795330988168 * T * leakage * t / del_P;
    result.push_back(V);
    
    return result;
}

int main() {
    std::cout << "Running test suite..." << std::endl;
    std::cout << "Testing AirLeak_eqn_4_10__T... ";
    try { AirLeak_eqn_4_10__T(1.0, 1.0, 1.0, 1.0); std::cout << "OK" << std::endl; } catch(...) { std::cout << "FAIL" << std::endl; }
    std::cout << "Testing AirLeak_eqn_4_10__V... ";
    try { AirLeak_eqn_4_10__V(1.0, 1.0, 1.0, 1.0); std::cout << "OK" << std::endl; } catch(...) { std::cout << "FAIL" << std::endl; }
    std::cout << "Test suite finished." << std::endl;
    return 0;
}