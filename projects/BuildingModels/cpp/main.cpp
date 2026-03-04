#include <iostream>

#include "kinematics.hpp"
#include "linearmomentumandthecentreofmass.hpp"

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
    std::cout << "Testing Kinematics_eqn_1_2__a... ";
    try { auto r = Kinematics_eqn_1_2__a(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_2__t... ";
    try { auto r = Kinematics_eqn_1_2__t(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_2__v0... ";
    try { auto r = Kinematics_eqn_1_2__v0(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_2__x... ";
    try { auto r = Kinematics_eqn_1_2__x(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_2__x0... ";
    try { auto r = Kinematics_eqn_1_2__x0(1.0, 1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_3__a... ";
    try { auto r = Kinematics_eqn_1_3__a(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_3__dx... ";
    try { auto r = Kinematics_eqn_1_3__dx(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_3__v... ";
    try { auto r = Kinematics_eqn_1_3__v(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing Kinematics_eqn_1_3__v0... ";
    try { auto r = Kinematics_eqn_1_3__v0(1.0, 1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_1__m... ";
    try { auto r = LinearMomentumAndTheCentreOfMass_eqn_10_1__m(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_1__p... ";
    try { auto r = LinearMomentumAndTheCentreOfMass_eqn_10_1__p(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_1__v... ";
    try { auto r = LinearMomentumAndTheCentreOfMass_eqn_10_1__v(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_11__F_ext... ";
    try { auto r = LinearMomentumAndTheCentreOfMass_eqn_10_11__F_ext(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_11__M... ";
    try { auto r = LinearMomentumAndTheCentreOfMass_eqn_10_11__M(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }
    std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_11__a_CM... ";
    try { auto r = LinearMomentumAndTheCentreOfMass_eqn_10_11__a_CM(1.0, 1.0); std::cout << "OK (" << r.size() << " results)" << std::endl; pass++; } catch(std::exception& e) { std::cout << "FAIL: " << e.what() << std::endl; fail++; } catch(...) { std::cout << "FAIL" << std::endl; fail++; }

    std::cout << "\n" << pass << " passed, " << fail << " failed." << std::endl;
    return fail > 0 ? 1 : 0;
}