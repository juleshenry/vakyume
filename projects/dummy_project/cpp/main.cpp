#include <iostream>

#include "basic.hpp"
#include "rotary.hpp"

int main() {
    std::cout << "Running test suite..." << std::endl;
    int pass = 0, fail = 0;
    std::cout << "Testing Basic_eqn_1_1__x... ";
    try { auto r = Basic_eqn_1_1__x(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Basic_eqn_1_1__y... ";
    try { auto r = Basic_eqn_1_1__y(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Basic_eqn_1_1__z... ";
    try { auto r = Basic_eqn_1_1__z(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__Q... ";
    try { auto r = Rotary_eqn_11_2__Q(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__Q_0... ";
    try { auto r = Rotary_eqn_11_2__Q_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__Q_external_gas_throughput... ";
    try { auto r = Rotary_eqn_11_2__Q_external_gas_throughput(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__SP_1... ";
    try { auto r = Rotary_eqn_11_2__SP_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__SP_2... ";
    try { auto r = Rotary_eqn_11_2__SP_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__S_vol_pump_speed... ";
    try { auto r = Rotary_eqn_11_2__S_vol_pump_speed(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__V... ";
    try { auto r = Rotary_eqn_11_2__V(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Rotary_eqn_11_2__t... ";
    try { auto r = Rotary_eqn_11_2__t(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }

    std::cout << "\n" << pass << " passed, " << fail << " failed." << std::endl;
    return fail > 0 ? 1 : 0;
}