#include <iostream>

#include "airleak.hpp"
#include "fluidflowvacuumlines.hpp"
#include "liquidring.hpp"
#include "precondensors.hpp"
#include "pressmgmt.hpp"
#include "processapp1.hpp"
#include "processapp2.hpp"
#include "rotarypistonvane.hpp"
#include "selectingpump.hpp"
#include "steamjetinjectors.hpp"
#include "vacuumtheory.hpp"

int main() {
  std::cout << "Running test suite..." << std::endl;
  int pass = 0, fail = 0;
  std::cout << "Testing AirLeak_eqn_4_10__T... ";
  try {
    auto r = AirLeak_eqn_4_10__T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing AirLeak_eqn_4_10__V... ";
  try {
    auto r = AirLeak_eqn_4_10__V(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing AirLeak_eqn_4_10__del_P... ";
  try {
    auto r = AirLeak_eqn_4_10__del_P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing AirLeak_eqn_4_10__leakage... ";
  try {
    auto r = AirLeak_eqn_4_10__leakage(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing AirLeak_eqn_4_7__W... ";
  try {
    auto r = AirLeak_eqn_4_7__W(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing AirLeak_eqn_4_7__W_T... ";
  try {
    auto r = AirLeak_eqn_4_7__W_T(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing AirLeak_eqn_4_7__sum_individual_leak_rates... ";
  try {
    auto r = AirLeak_eqn_4_7__sum_individual_leak_rates(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__D... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_1__D(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__Re... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_1__Re(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_1__mu(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__rho... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_1__rho(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_1__v... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_1__v(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_10__Suc_Pres... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_10__Suc_Pres(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_10__delta_P... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_10__delta_P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_10__oper_press... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_10__oper_press(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__D... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_11__D(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_11__L(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__f... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_11__f(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__g_c... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_11__g_c(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__h_r... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_11__h_r(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_11__v... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_11__v(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_12__L(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__d... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_12__d(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__delta_P... ";
  try {
    auto r =
        FluidFlowVacuumLines_eqn_2_12__delta_P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__f... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_12__f(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__g... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_12__g(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__rho... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_12__rho(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_12__v... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_12__v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_13__L(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__d... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_13__d(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__delta_P... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_13__delta_P(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__f... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_13__f(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__q... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_13__q(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_13__rho... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_13__rho(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__M... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_14__M(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__R... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_14__R(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__T... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_14__T(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__g_c... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_14__g_c(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__k... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_14__k(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_14__v_s... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_14__v_s(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_15__Re... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_15__Re(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_15__f... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_15__f(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_16__Re... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_16__Re(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_16__f... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_16__f(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17a__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17a__L(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17a__d... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17a__d(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17a__delta_P... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17a__delta_P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17a__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17a__mu(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17a__v... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17a__v(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17b__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17b__L(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17b__d... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17b__d(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17b__delta_P... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17b__delta_P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17b__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17b__mu(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_17b__q... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_17b__q(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_18a__D_eq... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_18a__D_eq(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_18a__R_ll... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_18a__R_ll(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_18b__R_ll... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_18b__R_ll(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_18b__h... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_18b__h(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_18b__w... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_18b__w(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__R_ll... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19a__R_ll(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__Re... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19a__Re(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19a__mu(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__rho... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19a__rho(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19a__v... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19a__v(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__Re... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19b__Re(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__h... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19b__h(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19b__mu(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__rho... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19b__rho(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__v... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19b__v(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_19b__w... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_19b__w(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_2__delta... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_2__delta(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_2__lambd... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_2__lambd(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_2__psi... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_2__psi(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_20__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_20__L(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout
      << "Testing FluidFlowVacuumLines_eqn_2_20__sum_equivalent_length... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_20__sum_equivalent_length(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_20__sum_pipe... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_20__sum_pipe(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_22__P_s... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_22__P_s(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_22__Q_throughput... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_22__Q_throughput(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_22__S_p... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_22__S_p(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_25__C(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__P_1... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_25__P_1(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__P_2... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_25__P_2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_25__Q_throughput... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_25__Q_throughput(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__D... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_26__D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_26__L(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__P_downstream... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_26__P_downstream(1.0, 1.0, 1.0, 1.0,
                                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__P_p... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_26__P_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__P_upstream... ";
  try {
    auto r =
        FluidFlowVacuumLines_eqn_2_26__P_upstream(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_26__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_26__q... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_26__q(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_28__C(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__D... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_28__D(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_28__L(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__P_p... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_28__P_p(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_28__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_28__mu(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_29__C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_29__C(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_29__S_1... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_29__S_1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_29__S_2... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_29__S_2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_3__D... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_3__D(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_3__kn... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_3__kn(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_3__lambd... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_3__lambd(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_31__C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_31__C(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_31__S_p... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_31__S_p(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_31__S_pump_speed... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_31__S_pump_speed(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_32__C_series... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_32__C_series(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_32__geometric_sum_C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_32__geometric_sum_C(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_33__C_paralell... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_33__C_paralell(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_33__arithmetic_sum_C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_33__arithmetic_sum_C(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_34__C(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__C_1... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_34__C_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__C_2... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_34__C_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__D... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_34__D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_34__L(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__P_p... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_34__P_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_34__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_34__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_35__C_L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_35__C_L(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_35__C_T... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_35__C_T(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_35__F_p... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_35__F_p(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_36__C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_36__C(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_36__C_0... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_36__C_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_36__F_t... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_36__F_t(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__A... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_37__A(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__C... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_37__C(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__F_t... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_37__F_t(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__M... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_37__M(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_37__T... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_37__T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_4___beta... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_4___beta(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_4__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_4__mu(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_4__vel_grad... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_4__vel_grad(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__D... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_5__D(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__L... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_5__L(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__delta_P... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_5__delta_P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_5__mu(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_5__q... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_5__q(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__lambd... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_6__lambd(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__mu... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_6__mu(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__rho... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_6__rho(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_6__v_a... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_6__v_a(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__T... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_7__T(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__k... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_7__k(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__m... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_7__m(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_7__v_a... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_7__v_a(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__M... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_8__M(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__P_c... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_8__P_c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__T_c... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_8__T_c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidFlowVacuumLines_eqn_2_8__mu_c... ";
  try {
    auto r = FluidFlowVacuumLines_eqn_2_8__mu_c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_1__D_r... ";
  try {
    auto r = LiquidRing_eqn_10_1__D_r(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_1__sig_R... ";
  try {
    auto r = LiquidRing_eqn_10_1__sig_R(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_1__w... ";
  try {
    auto r = LiquidRing_eqn_10_1__w(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_10__bhp... ";
  try {
    auto r = LiquidRing_eqn_10_10__bhp(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_10__bhp_0... ";
  try {
    auto r = LiquidRing_eqn_10_10__bhp_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_10__mu... ";
  try {
    auto r = LiquidRing_eqn_10_10__mu(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_10__rho... ";
  try {
    auto r = LiquidRing_eqn_10_10__rho(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_11__T_c... ";
  try {
    auto r = LiquidRing_eqn_10_11__T_c(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_11__T_s... ";
  try {
    auto r = LiquidRing_eqn_10_11__T_s(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_12__T_c... ";
  try {
    auto r = LiquidRing_eqn_10_12__T_c(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_12__T_s... ";
  try {
    auto r = LiquidRing_eqn_10_12__T_s(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_13__T_c... ";
  try {
    auto r = LiquidRing_eqn_10_13__T_c(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_13__T_s... ";
  try {
    auto r = LiquidRing_eqn_10_13__T_s(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_14__T_c... ";
  try {
    auto r = LiquidRing_eqn_10_14__T_c(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_14__T_s... ";
  try {
    auto r = LiquidRing_eqn_10_14__T_s(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_15__P... ";
  try {
    auto r = LiquidRing_eqn_10_15__P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_15__S_Th... ";
  try {
    auto r = LiquidRing_eqn_10_15__S_Th(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_15__S_p... ";
  try {
    auto r = LiquidRing_eqn_10_15__S_p(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_15__p_s... ";
  try {
    auto r = LiquidRing_eqn_10_15__p_s(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_16__P... ";
  try {
    auto r = LiquidRing_eqn_10_16__P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_16__S_0... ";
  try {
    auto r = LiquidRing_eqn_10_16__S_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_16__S_Th... ";
  try {
    auto r = LiquidRing_eqn_10_16__S_Th(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_16__p_0... ";
  try {
    auto r = LiquidRing_eqn_10_16__p_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_17__P... ";
  try {
    auto r = LiquidRing_eqn_10_17__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_17__S_0... ";
  try {
    auto r = LiquidRing_eqn_10_17__S_0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_17__S_Th... ";
  try {
    auto r = LiquidRing_eqn_10_17__S_Th(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_17__p_0... ";
  try {
    auto r = LiquidRing_eqn_10_17__p_0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_17__p_s... ";
  try {
    auto r = LiquidRing_eqn_10_17__p_s(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_18__P... ";
  try {
    auto r = LiquidRing_eqn_10_18__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_18__S_Th... ";
  try {
    auto r = LiquidRing_eqn_10_18__S_Th(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_18__S_p... ";
  try {
    auto r = LiquidRing_eqn_10_18__S_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_18__T_e... ";
  try {
    auto r = LiquidRing_eqn_10_18__T_e(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_18__T_i... ";
  try {
    auto r = LiquidRing_eqn_10_18__T_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_18__p_c... ";
  try {
    auto r = LiquidRing_eqn_10_18__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_18__p_s... ";
  try {
    auto r = LiquidRing_eqn_10_18__p_s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_19__P... ";
  try {
    auto r = LiquidRing_eqn_10_19__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_19__S_Th... ";
  try {
    auto r = LiquidRing_eqn_10_19__S_Th(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_19__S_p... ";
  try {
    auto r = LiquidRing_eqn_10_19__S_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_19__T_e... ";
  try {
    auto r = LiquidRing_eqn_10_19__T_e(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_19__T_i... ";
  try {
    auto r = LiquidRing_eqn_10_19__T_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_19__p_c... ";
  try {
    auto r = LiquidRing_eqn_10_19__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_19__p_s... ";
  try {
    auto r = LiquidRing_eqn_10_19__p_s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_2__PS... ";
  try {
    auto r = LiquidRing_eqn_10_2__PS(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_2__Q_gas... ";
  try {
    auto r = LiquidRing_eqn_10_2__Q_gas(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_2__V... ";
  try {
    auto r = LiquidRing_eqn_10_2__V(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_2__dP... ";
  try {
    auto r = LiquidRing_eqn_10_2__dP(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_2__dt... ";
  try {
    auto r = LiquidRing_eqn_10_2__dt(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__P... ";
  try {
    auto r = LiquidRing_eqn_10_20__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__S_0... ";
  try {
    auto r = LiquidRing_eqn_10_20__S_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__S_p... ";
  try {
    auto r = LiquidRing_eqn_10_20__S_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__T_e... ";
  try {
    auto r = LiquidRing_eqn_10_20__T_e(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__T_i... ";
  try {
    auto r = LiquidRing_eqn_10_20__T_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__p_0... ";
  try {
    auto r = LiquidRing_eqn_10_20__p_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__p_c... ";
  try {
    auto r = LiquidRing_eqn_10_20__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_20__p_s... ";
  try {
    auto r = LiquidRing_eqn_10_20__p_s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_21__P... ";
  try {
    auto r = LiquidRing_eqn_10_21__P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_21__P_d... ";
  try {
    auto r = LiquidRing_eqn_10_21__P_d(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_21__P_prime... ";
  try {
    auto r = LiquidRing_eqn_10_21__P_prime(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_3__N_mfw... ";
  try {
    auto r = LiquidRing_eqn_10_3__N_mfw(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_3__Q_gas... ";
  try {
    auto r = LiquidRing_eqn_10_3__Q_gas(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_3__T... ";
  try {
    auto r = LiquidRing_eqn_10_3__T(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_4__Q_gas... ";
  try {
    auto r = LiquidRing_eqn_10_4__Q_gas(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_4__SP_1... ";
  try {
    auto r = LiquidRing_eqn_10_4__SP_1(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_4__SP_2... ";
  try {
    auto r = LiquidRing_eqn_10_4__SP_2(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_4__S_p... ";
  try {
    auto r = LiquidRing_eqn_10_4__S_p(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_4__V... ";
  try {
    auto r = LiquidRing_eqn_10_4__V(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_4__t... ";
  try {
    auto r = LiquidRing_eqn_10_4__t(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_5__P_1... ";
  try {
    auto r = LiquidRing_eqn_10_5__P_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_5__P_2... ";
  try {
    auto r = LiquidRing_eqn_10_5__P_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_5__S_p... ";
  try {
    auto r = LiquidRing_eqn_10_5__S_p(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_5__V... ";
  try {
    auto r = LiquidRing_eqn_10_5__V(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_5__t... ";
  try {
    auto r = LiquidRing_eqn_10_5__t(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_6__P_1... ";
  try {
    auto r = LiquidRing_eqn_10_6__P_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_6__P_2... ";
  try {
    auto r = LiquidRing_eqn_10_6__P_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_6__S_a... ";
  try {
    auto r = LiquidRing_eqn_10_6__S_a(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_6__V... ";
  try {
    auto r = LiquidRing_eqn_10_6__V(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_6__t... ";
  try {
    auto r = LiquidRing_eqn_10_6__t(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_8__bhp... ";
  try {
    auto r = LiquidRing_eqn_10_8__bhp(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_8__c_p... ";
  try {
    auto r = LiquidRing_eqn_10_8__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_8__delta_T... ";
  try {
    auto r = LiquidRing_eqn_10_8__delta_T(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_8__delta_h_i... ";
  try {
    auto r = LiquidRing_eqn_10_8__delta_h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_8__f_a... ";
  try {
    auto r = LiquidRing_eqn_10_8__f_a(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_8__rho... ";
  try {
    auto r = LiquidRing_eqn_10_8__rho(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_8__w_i... ";
  try {
    auto r = LiquidRing_eqn_10_8__w_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_9__T_c... ";
  try {
    auto r = LiquidRing_eqn_10_9__T_c(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_9__T_s... ";
  try {
    auto r = LiquidRing_eqn_10_9__T_s(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LiquidRing_eqn_10_9__delta_T... ";
  try {
    auto r = LiquidRing_eqn_10_9__delta_T(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_1__P... ";
  try {
    auto r = Precondensors_eqn_7_1__P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_1__p_i... ";
  try {
    auto r = Precondensors_eqn_7_1__p_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_1__y_i... ";
  try {
    auto r = Precondensors_eqn_7_1__y_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_10__L_c_P... ";
  try {
    auto r = Precondensors_eqn_7_10__L_c_P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_10__Q_condensor_heat_duty... ";
  try {
    auto r = Precondensors_eqn_7_10__Q_condensor_heat_duty(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_10__del_T... ";
  try {
    auto r = Precondensors_eqn_7_10__del_T(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_11__Q_condensor_heat_duty... ";
  try {
    auto r = Precondensors_eqn_7_11__Q_condensor_heat_duty(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_11__U_v... ";
  try {
    auto r = Precondensors_eqn_7_11__U_v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_11__V_c... ";
  try {
    auto r = Precondensors_eqn_7_11__V_c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_11__del_T_LM... ";
  try {
    auto r = Precondensors_eqn_7_11__del_T_LM(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_12__A... ";
  try {
    auto r = Precondensors_eqn_7_12__A(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_12__Q_condensor_heat_duty... ";
  try {
    auto r = Precondensors_eqn_7_12__Q_condensor_heat_duty(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_12__U... ";
  try {
    auto r = Precondensors_eqn_7_12__U(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_12__del_T... ";
  try {
    auto r = Precondensors_eqn_7_12__del_T(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14a__A... ";
  try {
    auto r = Precondensors_eqn_7_14a__A(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14a__Q_condensor_heat_duty... ";
  try {
    auto r = Precondensors_eqn_7_14a__Q_condensor_heat_duty(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14a__U... ";
  try {
    auto r = Precondensors_eqn_7_14a__U(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14a__del_T_LM... ";
  try {
    auto r = Precondensors_eqn_7_14a__del_T_LM(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14b__A... ";
  try {
    auto r = Precondensors_eqn_7_14b__A(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14b__Q_condensor_heat_duty... ";
  try {
    auto r = Precondensors_eqn_7_14b__Q_condensor_heat_duty(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14b__U... ";
  try {
    auto r = Precondensors_eqn_7_14b__U(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14b__del_T_1... ";
  try {
    auto r = Precondensors_eqn_7_14b__del_T_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_14b__del_T_2... ";
  try {
    auto r = Precondensors_eqn_7_14b__del_T_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_15__U... ";
  try {
    auto r = Precondensors_eqn_7_15__U(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_15__sum_R... ";
  try {
    auto r = Precondensors_eqn_7_15__sum_R(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__D_0... ";
  try {
    auto r = Precondensors_eqn_7_16__D_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__D_LM... ";
  try {
    auto r = Precondensors_eqn_7_16__D_LM(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__D_i... ";
  try {
    auto r = Precondensors_eqn_7_16__D_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__R_f_0... ";
  try {
    auto r = Precondensors_eqn_7_16__R_f_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                           1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__R_fi... ";
  try {
    auto r = Precondensors_eqn_7_16__R_fi(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__U_0... ";
  try {
    auto r = Precondensors_eqn_7_16__U_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__h_0... ";
  try {
    auto r = Precondensors_eqn_7_16__h_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__h_i... ";
  try {
    auto r = Precondensors_eqn_7_16__h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__k_w... ";
  try {
    auto r = Precondensors_eqn_7_16__k_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_16__x_w... ";
  try {
    auto r = Precondensors_eqn_7_16__x_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_17__R_0... ";
  try {
    auto r = Precondensors_eqn_7_17__R_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_17__R_nc... ";
  try {
    auto r = Precondensors_eqn_7_17__R_nc(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_17__h_c... ";
  try {
    auto r = Precondensors_eqn_7_17__h_c(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__D_0... ";
  try {
    auto r = Precondensors_eqn_7_18__D_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__D_LM... ";
  try {
    auto r = Precondensors_eqn_7_18__D_LM(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__D_i... ";
  try {
    auto r = Precondensors_eqn_7_18__D_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__R_fi... ";
  try {
    auto r = Precondensors_eqn_7_18__R_fi(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__R_fo... ";
  try {
    auto r = Precondensors_eqn_7_18__R_fo(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__R_nc... ";
  try {
    auto r = Precondensors_eqn_7_18__R_nc(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__U_0... ";
  try {
    auto r = Precondensors_eqn_7_18__U_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__h_c... ";
  try {
    auto r = Precondensors_eqn_7_18__h_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__h_i... ";
  try {
    auto r = Precondensors_eqn_7_18__h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__k_w... ";
  try {
    auto r = Precondensors_eqn_7_18__k_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_18__x_w... ";
  try {
    auto r = Precondensors_eqn_7_18__x_w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                         1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_2__P_i_0... ";
  try {
    auto r = Precondensors_eqn_7_2__P_i_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_2__p_i... ";
  try {
    auto r = Precondensors_eqn_7_2__p_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_2__x_i... ";
  try {
    auto r = Precondensors_eqn_7_2__x_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_3__P_i_0... ";
  try {
    auto r = Precondensors_eqn_7_3__P_i_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_3__epsilon_i... ";
  try {
    auto r = Precondensors_eqn_7_3__epsilon_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_3__p_i... ";
  try {
    auto r = Precondensors_eqn_7_3__p_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_3__x_i... ";
  try {
    auto r = Precondensors_eqn_7_3__x_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4a__P... ";
  try {
    auto r = Precondensors_eqn_7_4a__P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4a__p_c... ";
  try {
    auto r = Precondensors_eqn_7_4a__p_c(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4a__p_nc... ";
  try {
    auto r = Precondensors_eqn_7_4a__p_nc(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4aa__n_i... ";
  try {
    auto r = Precondensors_eqn_7_4aa__n_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4aa__n_nc... ";
  try {
    auto r = Precondensors_eqn_7_4aa__n_nc(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4aa__p_i... ";
  try {
    auto r = Precondensors_eqn_7_4aa__p_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4aa__p_nc... ";
  try {
    auto r = Precondensors_eqn_7_4aa__p_nc(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ab__P_c... ";
  try {
    auto r = Precondensors_eqn_7_4ab__P_c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ab__p... ";
  try {
    auto r = Precondensors_eqn_7_4ab__p(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ab__p_i... ";
  try {
    auto r = Precondensors_eqn_7_4ab__p_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ab__p_nc... ";
  try {
    auto r = Precondensors_eqn_7_4ab__p_nc(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ac__P_c... ";
  try {
    auto r = Precondensors_eqn_7_4ac__P_c(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ac__n_i... ";
  try {
    auto r = Precondensors_eqn_7_4ac__n_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ac__n_nc... ";
  try {
    auto r = Precondensors_eqn_7_4ac__n_nc(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ac__p... ";
  try {
    auto r = Precondensors_eqn_7_4ac__p(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_4ac__p_i... ";
  try {
    auto r = Precondensors_eqn_7_4ac__p_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_5__N_i... ";
  try {
    auto r = Precondensors_eqn_7_5__N_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_5__N_nc... ";
  try {
    auto r = Precondensors_eqn_7_5__N_nc(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_5__P... ";
  try {
    auto r = Precondensors_eqn_7_5__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_5__P_c... ";
  try {
    auto r = Precondensors_eqn_7_5__P_c(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_5__p_i... ";
  try {
    auto r = Precondensors_eqn_7_5__p_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_6__M... ";
  try {
    auto r = Precondensors_eqn_7_6__M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_6__P... ";
  try {
    auto r = Precondensors_eqn_7_6__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_6__P_i_0... ";
  try {
    auto r = Precondensors_eqn_7_6__P_i_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_6__W_air... ";
  try {
    auto r = Precondensors_eqn_7_6__W_air(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_6__W_i... ";
  try {
    auto r = Precondensors_eqn_7_6__W_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_6__p_c... ";
  try {
    auto r = Precondensors_eqn_7_6__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_6__x_i... ";
  try {
    auto r = Precondensors_eqn_7_6__x_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__M... ";
  try {
    auto r = Precondensors_eqn_7_7__M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__P... ";
  try {
    auto r = Precondensors_eqn_7_7__P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__P_i_0... ";
  try {
    auto r = Precondensors_eqn_7_7__P_i_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__W_air... ";
  try {
    auto r = Precondensors_eqn_7_7__W_air(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__W_i... ";
  try {
    auto r = Precondensors_eqn_7_7__W_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__epsilon_i... ";
  try {
    auto r =
        Precondensors_eqn_7_7__epsilon_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__p_c... ";
  try {
    auto r = Precondensors_eqn_7_7__p_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_7__x_i... ";
  try {
    auto r = Precondensors_eqn_7_7__x_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_8__L_c... ";
  try {
    auto r = Precondensors_eqn_7_8__L_c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_8__Q_condensor_heat_duty... ";
  try {
    auto r = Precondensors_eqn_7_8__Q_condensor_heat_duty(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_8__c_p... ";
  try {
    auto r = Precondensors_eqn_7_8__c_p(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_8__del_T... ";
  try {
    auto r = Precondensors_eqn_7_8__del_T(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_9__L_c... ";
  try {
    auto r = Precondensors_eqn_7_9__L_c(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_9__Q_condensor_heat_duty... ";
  try {
    auto r = Precondensors_eqn_7_9__Q_condensor_heat_duty(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_9__c_p... ";
  try {
    auto r = Precondensors_eqn_7_9__c_p(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_9__del_T... ";
  try {
    auto r = Precondensors_eqn_7_9__del_T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Precondensors_eqn_7_9__rho... ";
  try {
    auto r = Precondensors_eqn_7_9__rho(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_1__Abs_Pressure... ";
  try {
    auto r = PressMgmt_eqn_3_1__Abs_Pressure(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_1__BarometricPressure... ";
  try {
    auto r = PressMgmt_eqn_3_1__BarometricPressure(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_1__Vacuum... ";
  try {
    auto r = PressMgmt_eqn_3_1__Vacuum(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_11__A_C... ";
  try {
    auto r = PressMgmt_eqn_3_11__A_C(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_11__H_2... ";
  try {
    auto r = PressMgmt_eqn_3_11__H_2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_11__P... ";
  try {
    auto r = PressMgmt_eqn_3_11__P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_11__V... ";
  try {
    auto r = PressMgmt_eqn_3_11__V(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_12__H_2... ";
  try {
    auto r = PressMgmt_eqn_3_12__H_2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_12__KAPPA_1... ";
  try {
    auto r = PressMgmt_eqn_3_12__KAPPA_1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_12__P... ";
  try {
    auto r = PressMgmt_eqn_3_12__P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_13__H_1... ";
  try {
    auto r = PressMgmt_eqn_3_13__H_1(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_13__H_2... ";
  try {
    auto r = PressMgmt_eqn_3_13__H_2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_13__KAPPA_2... ";
  try {
    auto r = PressMgmt_eqn_3_13__KAPPA_2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_13__P... ";
  try {
    auto r = PressMgmt_eqn_3_13__P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_15__V_PMIN... ";
  try {
    auto r = PressMgmt_eqn_3_15__V_PMIN();
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_16__V_div_V_P_MAX... ";
  try {
    auto r = PressMgmt_eqn_3_16__V_div_V_P_MAX();
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_17__P_MIN... ";
  try {
    auto r = PressMgmt_eqn_3_17__P_MIN();
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_2__G... ";
  try {
    auto r = PressMgmt_eqn_3_2__G(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_2__G_C... ";
  try {
    auto r = PressMgmt_eqn_3_2__G_C(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_2__H... ";
  try {
    auto r = PressMgmt_eqn_3_2__H(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_2__P... ";
  try {
    auto r = PressMgmt_eqn_3_2__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_2__rho... ";
  try {
    auto r = PressMgmt_eqn_3_2__rho(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_3__H_1... ";
  try {
    auto r = PressMgmt_eqn_3_3__H_1(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_3__H_2... ";
  try {
    auto r = PressMgmt_eqn_3_3__H_2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_3__P... ";
  try {
    auto r = PressMgmt_eqn_3_3__P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_3__P_P... ";
  try {
    auto r = PressMgmt_eqn_3_3__P_P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_4__KAPPA... ";
  try {
    auto r = PressMgmt_eqn_3_4__KAPPA(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_4__P... ";
  try {
    auto r = PressMgmt_eqn_3_4__P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_4__V... ";
  try {
    auto r = PressMgmt_eqn_3_4__V(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_5__P... ";
  try {
    auto r = PressMgmt_eqn_3_5__P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_5__P_P... ";
  try {
    auto r = PressMgmt_eqn_3_5__P_P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_5__V... ";
  try {
    auto r = PressMgmt_eqn_3_5__V(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_5__V_P... ";
  try {
    auto r = PressMgmt_eqn_3_5__V_P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_6__H_1... ";
  try {
    auto r = PressMgmt_eqn_3_6__H_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_6__H_2... ";
  try {
    auto r = PressMgmt_eqn_3_6__H_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_6__P... ";
  try {
    auto r = PressMgmt_eqn_3_6__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_6__V... ";
  try {
    auto r = PressMgmt_eqn_3_6__V(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_6__V_P... ";
  try {
    auto r = PressMgmt_eqn_3_6__V_P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_8__A_C... ";
  try {
    auto r = PressMgmt_eqn_3_8__A_C(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_8__H_2... ";
  try {
    auto r = PressMgmt_eqn_3_8__H_2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_8__V_P... ";
  try {
    auto r = PressMgmt_eqn_3_8__V_P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_9__A_C... ";
  try {
    auto r = PressMgmt_eqn_3_9__A_C(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_9__H_1... ";
  try {
    auto r = PressMgmt_eqn_3_9__H_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_9__H_2... ";
  try {
    auto r = PressMgmt_eqn_3_9__H_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_9__P... ";
  try {
    auto r = PressMgmt_eqn_3_9__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PressMgmt_eqn_3_9__V... ";
  try {
    auto r = PressMgmt_eqn_3_9__V(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_1__K_i... ";
  try {
    auto r = ProcessApp1_eqn_5_1__K_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_1__x_i... ";
  try {
    auto r = ProcessApp1_eqn_5_1__x_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_1__y_i... ";
  try {
    auto r = ProcessApp1_eqn_5_1__y_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10a__D... ";
  try {
    auto r = ProcessApp1_eqn_5_10a__D(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10a__L_0... ";
  try {
    auto r = ProcessApp1_eqn_5_10a__L_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10a__V_1... ";
  try {
    auto r = ProcessApp1_eqn_5_10a__V_1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10b__L_0... ";
  try {
    auto r = ProcessApp1_eqn_5_10b__L_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10b__R... ";
  try {
    auto r = ProcessApp1_eqn_5_10b__R(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10b__V_1... ";
  try {
    auto r = ProcessApp1_eqn_5_10b__V_1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10c__D... ";
  try {
    auto r = ProcessApp1_eqn_5_10c__D(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10c__L_0... ";
  try {
    auto r = ProcessApp1_eqn_5_10c__L_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_10c__R... ";
  try {
    auto r = ProcessApp1_eqn_5_10c__R(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_11__B... ";
  try {
    auto r = ProcessApp1_eqn_5_11__B(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_11__L_N... ";
  try {
    auto r = ProcessApp1_eqn_5_11__L_N(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_11__V_0... ";
  try {
    auto r = ProcessApp1_eqn_5_11__V_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_12__Eff... ";
  try {
    auto r = ProcessApp1_eqn_5_12__Eff(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_12__N_ES... ";
  try {
    auto r = ProcessApp1_eqn_5_12__N_ES(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_12__N_t... ";
  try {
    auto r = ProcessApp1_eqn_5_12__N_t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_12__T... ";
  try {
    auto r = ProcessApp1_eqn_5_12__T(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_13__HETP... ";
  try {
    auto r = ProcessApp1_eqn_5_13__HETP(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_13__H_p... ";
  try {
    auto r = ProcessApp1_eqn_5_13__H_p(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_13__N_ES... ";
  try {
    auto r = ProcessApp1_eqn_5_13__N_ES(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_14__M... ";
  try {
    auto r = ProcessApp1_eqn_5_14__M(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_14__P_0... ";
  try {
    auto r = ProcessApp1_eqn_5_14__P_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_14__T... ";
  try {
    auto r = ProcessApp1_eqn_5_14__T(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_14__W_E... ";
  try {
    auto r = ProcessApp1_eqn_5_14__W_E(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_15__M_1... ";
  try {
    auto r = ProcessApp1_eqn_5_15__M_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_15__M_2... ";
  try {
    auto r = ProcessApp1_eqn_5_15__M_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_15__P_0_1... ";
  try {
    auto r = ProcessApp1_eqn_5_15__P_0_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_15__P_0_2... ";
  try {
    auto r = ProcessApp1_eqn_5_15__P_0_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_15__a_M_12... ";
  try {
    auto r = ProcessApp1_eqn_5_15__a_M_12(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_16__H_i... ";
  try {
    auto r = ProcessApp1_eqn_5_16__H_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_16__p_i... ";
  try {
    auto r = ProcessApp1_eqn_5_16__p_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_16__x_i... ";
  try {
    auto r = ProcessApp1_eqn_5_16__x_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_17__H_2_1... ";
  try {
    auto r = ProcessApp1_eqn_5_17__H_2_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_17__H_2_3... ";
  try {
    auto r = ProcessApp1_eqn_5_17__H_2_3(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_17__H_2_mi... ";
  try {
    auto r = ProcessApp1_eqn_5_17__H_2_mi(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_17__x_1... ";
  try {
    auto r = ProcessApp1_eqn_5_17__x_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_17__x_3... ";
  try {
    auto r = ProcessApp1_eqn_5_17__x_3(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2a__K_1... ";
  try {
    auto r = ProcessApp1_eqn_5_2a__K_1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2a__K_2... ";
  try {
    auto r = ProcessApp1_eqn_5_2a__K_2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2a__alpha_1_2... ";
  try {
    auto r = ProcessApp1_eqn_5_2a__alpha_1_2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2b__K_1... ";
  try {
    auto r = ProcessApp1_eqn_5_2b__K_1(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2b__K_2... ";
  try {
    auto r = ProcessApp1_eqn_5_2b__K_2(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2b__x_1... ";
  try {
    auto r = ProcessApp1_eqn_5_2b__x_1(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2b__x_2... ";
  try {
    auto r = ProcessApp1_eqn_5_2b__x_2(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2b__y_1... ";
  try {
    auto r = ProcessApp1_eqn_5_2b__y_1(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_2b__y_2... ";
  try {
    auto r = ProcessApp1_eqn_5_2b__y_2(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_3__P_0_i... ";
  try {
    auto r = ProcessApp1_eqn_5_3__P_0_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_3__p_i... ";
  try {
    auto r = ProcessApp1_eqn_5_3__p_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_3__x_i... ";
  try {
    auto r = ProcessApp1_eqn_5_3__x_i(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_4__P... ";
  try {
    auto r = ProcessApp1_eqn_5_4__P(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_4__P_0_i... ";
  try {
    auto r = ProcessApp1_eqn_5_4__P_0_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_4__x_i... ";
  try {
    auto r = ProcessApp1_eqn_5_4__x_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_4__y_i... ";
  try {
    auto r = ProcessApp1_eqn_5_4__y_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_5__P_0_1... ";
  try {
    auto r = ProcessApp1_eqn_5_5__P_0_1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_5__P_0_2... ";
  try {
    auto r = ProcessApp1_eqn_5_5__P_0_2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_5__alpha_12... ";
  try {
    auto r = ProcessApp1_eqn_5_5__alpha_12(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_6__P_0_i... ";
  try {
    auto r = ProcessApp1_eqn_5_6__P_0_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_6__gamma_i... ";
  try {
    auto r = ProcessApp1_eqn_5_6__gamma_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_6__p_i... ";
  try {
    auto r = ProcessApp1_eqn_5_6__p_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_6__x_i... ";
  try {
    auto r = ProcessApp1_eqn_5_6__x_i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_7__P... ";
  try {
    auto r = ProcessApp1_eqn_5_7__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_7__P_0_i... ";
  try {
    auto r = ProcessApp1_eqn_5_7__P_0_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_7__gamma_i... ";
  try {
    auto r = ProcessApp1_eqn_5_7__gamma_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_7__x_i... ";
  try {
    auto r = ProcessApp1_eqn_5_7__x_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_7__y_i... ";
  try {
    auto r = ProcessApp1_eqn_5_7__y_i(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_8__P_0_1... ";
  try {
    auto r = ProcessApp1_eqn_5_8__P_0_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_8__P_0_2... ";
  try {
    auto r = ProcessApp1_eqn_5_8__P_0_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_8__alpha_12... ";
  try {
    auto r = ProcessApp1_eqn_5_8__alpha_12(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_8__gamma_1... ";
  try {
    auto r = ProcessApp1_eqn_5_8__gamma_1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_8__gamma_2... ";
  try {
    auto r = ProcessApp1_eqn_5_8__gamma_2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_9__D... ";
  try {
    auto r = ProcessApp1_eqn_5_9__D(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_9__L_0... ";
  try {
    auto r = ProcessApp1_eqn_5_9__L_0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp1_eqn_5_9__V_1... ";
  try {
    auto r = ProcessApp1_eqn_5_9__V_1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__T_1... ";
  try {
    auto r = ProcessApp2_eqn_6_1__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__T_2... ";
  try {
    auto r = ProcessApp2_eqn_6_1__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__T_R... ";
  try {
    auto r = ProcessApp2_eqn_6_1__T_R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__c_p... ";
  try {
    auto r = ProcessApp2_eqn_6_1__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__del_h_v... ";
  try {
    auto r = ProcessApp2_eqn_6_1__del_h_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__w_1... ";
  try {
    auto r = ProcessApp2_eqn_6_1__w_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__w_2... ";
  try {
    auto r = ProcessApp2_eqn_6_1__w_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_1__w_v... ";
  try {
    auto r = ProcessApp2_eqn_6_1__w_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_10__A... ";
  try {
    auto r = ProcessApp2_eqn_6_10__A(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_10__dV_dt... ";
  try {
    auto r = ProcessApp2_eqn_6_10__dV_dt(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_10__delta_P... ";
  try {
    auto r = ProcessApp2_eqn_6_10__delta_P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_10__mu... ";
  try {
    auto r = ProcessApp2_eqn_6_10__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_10__r_c... ";
  try {
    auto r = ProcessApp2_eqn_6_10__r_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_10__s... ";
  try {
    auto r = ProcessApp2_eqn_6_10__s(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_10__tau... ";
  try {
    auto r = ProcessApp2_eqn_6_10__tau(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_11a__A_d... ";
  try {
    auto r = ProcessApp2_eqn_6_11a__A_d(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_11a__delta_T... ";
  try {
    auto r = ProcessApp2_eqn_6_11a__delta_T(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_11a__delta_h_i... ";
  try {
    auto r = ProcessApp2_eqn_6_11a__delta_h_i(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_11a__delta_m... ";
  try {
    auto r = ProcessApp2_eqn_6_11a__delta_m(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_11a__h_d... ";
  try {
    auto r = ProcessApp2_eqn_6_11a__h_d(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_11a__m_b... ";
  try {
    auto r = ProcessApp2_eqn_6_11a__m_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_11a__t_R... ";
  try {
    auto r = ProcessApp2_eqn_6_11a__t_R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_2__Q_v... ";
  try {
    auto r = ProcessApp2_eqn_6_2__Q_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_2__T_1... ";
  try {
    auto r = ProcessApp2_eqn_6_2__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_2__T_2... ";
  try {
    auto r = ProcessApp2_eqn_6_2__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_2__T_R... ";
  try {
    auto r = ProcessApp2_eqn_6_2__T_R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_2__c_p... ";
  try {
    auto r = ProcessApp2_eqn_6_2__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_2__w_1... ";
  try {
    auto r = ProcessApp2_eqn_6_2__w_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_2__w_2... ";
  try {
    auto r = ProcessApp2_eqn_6_2__w_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_4__Q_v... ";
  try {
    auto r = ProcessApp2_eqn_6_4__Q_v(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_4__delta_h_v... ";
  try {
    auto r = ProcessApp2_eqn_6_4__delta_h_v(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_4__w_v... ";
  try {
    auto r = ProcessApp2_eqn_6_4__w_v(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__C_1... ";
  try {
    auto r = ProcessApp2_eqn_6_7__C_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__C_2... ";
  try {
    auto r = ProcessApp2_eqn_6_7__C_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__T_1... ";
  try {
    auto r = ProcessApp2_eqn_6_7__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__T_2... ";
  try {
    auto r = ProcessApp2_eqn_6_7__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__c_p... ";
  try {
    auto r = ProcessApp2_eqn_6_7__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__delta_h_c... ";
  try {
    auto r =
        ProcessApp2_eqn_6_7__delta_h_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__delta_h_v... ";
  try {
    auto r =
        ProcessApp2_eqn_6_7__delta_h_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__m_b... ";
  try {
    auto r = ProcessApp2_eqn_6_7__m_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_7__m_v... ";
  try {
    auto r = ProcessApp2_eqn_6_7__m_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__C_1... ";
  try {
    auto r =
        ProcessApp2_eqn_6_8__C_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__C_2... ";
  try {
    auto r =
        ProcessApp2_eqn_6_8__C_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__T_1... ";
  try {
    auto r =
        ProcessApp2_eqn_6_8__T_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__T_2... ";
  try {
    auto r =
        ProcessApp2_eqn_6_8__T_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__c_p... ";
  try {
    auto r =
        ProcessApp2_eqn_6_8__c_p(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__delta_h_c... ";
  try {
    auto r = ProcessApp2_eqn_6_8__delta_h_c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__delta_h_v... ";
  try {
    auto r = ProcessApp2_eqn_6_8__delta_h_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                            1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__delta_t... ";
  try {
    auto r = ProcessApp2_eqn_6_8__delta_t(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__m_b... ";
  try {
    auto r =
        ProcessApp2_eqn_6_8__m_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_8__w_v... ";
  try {
    auto r =
        ProcessApp2_eqn_6_8__w_v(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_9__A... ";
  try {
    auto r = ProcessApp2_eqn_6_9__A(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_9__dV_dt... ";
  try {
    auto r = ProcessApp2_eqn_6_9__dV_dt(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_9__delta_P... ";
  try {
    auto r = ProcessApp2_eqn_6_9__delta_P(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_9__m... ";
  try {
    auto r = ProcessApp2_eqn_6_9__m(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_9__mu... ";
  try {
    auto r = ProcessApp2_eqn_6_9__mu(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_9__r... ";
  try {
    auto r = ProcessApp2_eqn_6_9__r(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ProcessApp2_eqn_6_9__r_M... ";
  try {
    auto r = ProcessApp2_eqn_6_9__r_M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_1__PS... ";
  try {
    auto r = RotaryPistonVane_eqn_11_1__PS(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_1__Q_0... ";
  try {
    auto r = RotaryPistonVane_eqn_11_1__Q_0(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout
      << "Testing RotaryPistonVane_eqn_11_1__Q_external_gas_throughput... ";
  try {
    auto r = RotaryPistonVane_eqn_11_1__Q_external_gas_throughput(1.0, 1.0, 1.0,
                                                                  1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_1__V... ";
  try {
    auto r = RotaryPistonVane_eqn_11_1__V(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_1__dP... ";
  try {
    auto r = RotaryPistonVane_eqn_11_1__dP(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_1__dT... ";
  try {
    auto r = RotaryPistonVane_eqn_11_1__dT(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_2__Q... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__Q(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_2__Q_0... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__Q_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout
      << "Testing RotaryPistonVane_eqn_11_2__Q_external_gas_throughput... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__Q_external_gas_throughput(
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_2__SP_1... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__SP_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_2__SP_2... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__SP_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_2__S_vol_pump_speed... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__S_vol_pump_speed(1.0, 1.0, 1.0, 1.0,
                                                         1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_2__V... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__V(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_2__t... ";
  try {
    auto r = RotaryPistonVane_eqn_11_2__t(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_3__F_s... ";
  try {
    auto r = RotaryPistonVane_eqn_11_3__F_s(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_3__t... ";
  try {
    auto r = RotaryPistonVane_eqn_11_3__t(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_3__t_c... ";
  try {
    auto r = RotaryPistonVane_eqn_11_3__t_c(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_4__p_g... ";
  try {
    auto r = RotaryPistonVane_eqn_11_4__p_g(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_4__p_s... ";
  try {
    auto r = RotaryPistonVane_eqn_11_4__p_s(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_4__p_v... ";
  try {
    auto r = RotaryPistonVane_eqn_11_4__p_v(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_5__P_0_v... ";
  try {
    auto r = RotaryPistonVane_eqn_11_5__P_0_v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_5__P_D... ";
  try {
    auto r = RotaryPistonVane_eqn_11_5__P_D(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_5__p_g... ";
  try {
    auto r = RotaryPistonVane_eqn_11_5__p_g(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_5__p_v_max... ";
  try {
    auto r = RotaryPistonVane_eqn_11_5__p_v_max(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__P_0_V... ";
  try {
    auto r =
        RotaryPistonVane_eqn_11_6__P_0_V(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__P_D... ";
  try {
    auto r = RotaryPistonVane_eqn_11_6__P_D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__P_v_0... ";
  try {
    auto r =
        RotaryPistonVane_eqn_11_6__P_v_0(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__S_B... ";
  try {
    auto r = RotaryPistonVane_eqn_11_6__S_B(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__S_D... ";
  try {
    auto r = RotaryPistonVane_eqn_11_6__S_D(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__p_b... ";
  try {
    auto r = RotaryPistonVane_eqn_11_6__p_b(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__p_g... ";
  try {
    auto r = RotaryPistonVane_eqn_11_6__p_g(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotaryPistonVane_eqn_11_6__p_v_max... ";
  try {
    auto r =
        RotaryPistonVane_eqn_11_6__p_v_max(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_1__NC... ";
  try {
    auto r = SelectingPump_eqn_8_1__NC(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_1__NS... ";
  try {
    auto r = SelectingPump_eqn_8_1__NS(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_1__SCON... ";
  try {
    auto r = SelectingPump_eqn_8_1__SCON(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_1__installation_cost... ";
  try {
    auto r = SelectingPump_eqn_8_1__installation_cost(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_2__hp... ";
  try {
    auto r = SelectingPump_eqn_8_2__hp(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_2__installed_costs... ";
  try {
    auto r = SelectingPump_eqn_8_2__installed_costs(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_3__hp... ";
  try {
    auto r = SelectingPump_eqn_8_3__hp(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_3__installed_costs... ";
  try {
    auto r = SelectingPump_eqn_8_3__installed_costs(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_4__hp... ";
  try {
    auto r = SelectingPump_eqn_8_4__hp(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_4__installed_costs... ";
  try {
    auto r = SelectingPump_eqn_8_4__installed_costs(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_5__Eff... ";
  try {
    auto r = SelectingPump_eqn_8_5__Eff(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_5__actual_brake_horsepower... ";
  try {
    auto r = SelectingPump_eqn_8_5__actual_brake_horsepower(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout
      << "Testing SelectingPump_eqn_8_5__theoretical_adiabatic_horsepower... ";
  try {
    auto r = SelectingPump_eqn_8_5__theoretical_adiabatic_horsepower(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__M... ";
  try {
    auto r = SelectingPump_eqn_8_6__M(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__P_1... ";
  try {
    auto r = SelectingPump_eqn_8_6__P_1(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__P_2... ";
  try {
    auto r = SelectingPump_eqn_8_6__P_2(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__R... ";
  try {
    auto r = SelectingPump_eqn_8_6__R(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__T... ";
  try {
    auto r = SelectingPump_eqn_8_6__T(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__adiabatic_hp... ";
  try {
    auto r =
        SelectingPump_eqn_8_6__adiabatic_hp(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__k... ";
  try {
    auto r = SelectingPump_eqn_8_6__k(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_6__w... ";
  try {
    auto r = SelectingPump_eqn_8_6__w(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_7__P_1... ";
  try {
    auto r = SelectingPump_eqn_8_7__P_1(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_7__P_2... ";
  try {
    auto r = SelectingPump_eqn_8_7__P_2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_7__adiabatic_hp... ";
  try {
    auto r = SelectingPump_eqn_8_7__adiabatic_hp(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_7__w... ";
  try {
    auto r = SelectingPump_eqn_8_7__w(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_8__P_1... ";
  try {
    auto r = SelectingPump_eqn_8_8__P_1(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_8__P_2... ";
  try {
    auto r = SelectingPump_eqn_8_8__P_2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_8__adiabatic_power_watts... ";
  try {
    auto r = SelectingPump_eqn_8_8__adiabatic_power_watts(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_8__f... ";
  try {
    auto r = SelectingPump_eqn_8_8__f(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_9__E_j... ";
  try {
    auto r = SelectingPump_eqn_8_9__E_j(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_9__E_m... ";
  try {
    auto r = SelectingPump_eqn_8_9__E_m(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_9__e... ";
  try {
    auto r = SelectingPump_eqn_8_9__e(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_9__r... ";
  try {
    auto r = SelectingPump_eqn_8_9__r(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SelectingPump_eqn_8_9__s... ";
  try {
    auto r = SelectingPump_eqn_8_9__s(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_1__A... ";
  try {
    auto r = SteamJetInjectors_eqn_9_1__A(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_1__rho_s... ";
  try {
    auto r = SteamJetInjectors_eqn_9_1__rho_s(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_1__v... ";
  try {
    auto r = SteamJetInjectors_eqn_9_1__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_1__w_s... ";
  try {
    auto r = SteamJetInjectors_eqn_9_1__w_s(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_2__P_m... ";
  try {
    auto r = SteamJetInjectors_eqn_9_2__P_m(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_2__d_n... ";
  try {
    auto r = SteamJetInjectors_eqn_9_2__d_n(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_2__rho_s... ";
  try {
    auto r = SteamJetInjectors_eqn_9_2__rho_s(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_2__w_s... ";
  try {
    auto r = SteamJetInjectors_eqn_9_2__w_s(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_3__P_s... ";
  try {
    auto r = SteamJetInjectors_eqn_9_3__P_s(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_3__V... ";
  try {
    auto r = SteamJetInjectors_eqn_9_3__V(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_3__t_e... ";
  try {
    auto r = SteamJetInjectors_eqn_9_3__t_e(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_3__w_j... ";
  try {
    auto r = SteamJetInjectors_eqn_9_3__w_j(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_4__AEL... ";
  try {
    auto r = SteamJetInjectors_eqn_9_4__AEL(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_4__SC... ";
  try {
    auto r = SteamJetInjectors_eqn_9_4__SC(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_4__r... ";
  try {
    auto r = SteamJetInjectors_eqn_9_4__r(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_4__w_s... ";
  try {
    auto r = SteamJetInjectors_eqn_9_4__w_s(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_5__V... ";
  try {
    auto r = SteamJetInjectors_eqn_9_5__V(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_5__r_h... ";
  try {
    auto r = SteamJetInjectors_eqn_9_5__r_h(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_5__t_h... ";
  try {
    auto r = SteamJetInjectors_eqn_9_5__t_h(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SteamJetInjectors_eqn_9_5__w_h... ";
  try {
    auto r = SteamJetInjectors_eqn_9_5__w_h(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_10__P_1... ";
  try {
    auto r = VacuumTheory_eqn_1_10__P_1(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_10__P_2... ";
  try {
    auto r = VacuumTheory_eqn_1_10__P_2(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_10__T_1... ";
  try {
    auto r = VacuumTheory_eqn_1_10__T_1(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_10__T_2... ";
  try {
    auto r = VacuumTheory_eqn_1_10__T_2(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_10__V_1... ";
  try {
    auto r = VacuumTheory_eqn_1_10__V_1(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_10__V_2... ";
  try {
    auto r = VacuumTheory_eqn_1_10__V_2(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_11__M... ";
  try {
    auto r = VacuumTheory_eqn_1_11__M(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_11__P... ";
  try {
    auto r = VacuumTheory_eqn_1_11__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_11__T... ";
  try {
    auto r = VacuumTheory_eqn_1_11__T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_11__W... ";
  try {
    auto r = VacuumTheory_eqn_1_11__W(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_11__q... ";
  try {
    auto r = VacuumTheory_eqn_1_11__q(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_12__Total_P... ";
  try {
    auto r = VacuumTheory_eqn_1_12__Total_P(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_12__sum_partial_pressures... ";
  try {
    auto r = VacuumTheory_eqn_1_12__sum_partial_pressures(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_13a__n... ";
  try {
    auto r = VacuumTheory_eqn_1_13a__n(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_13a__n_a... ";
  try {
    auto r = VacuumTheory_eqn_1_13a__n_a(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_13a__y_a... ";
  try {
    auto r = VacuumTheory_eqn_1_13a__y_a(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_13b__P... ";
  try {
    auto r = VacuumTheory_eqn_1_13b__P(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_13b__p_a... ";
  try {
    auto r = VacuumTheory_eqn_1_13b__p_a(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_13b__y_a... ";
  try {
    auto r = VacuumTheory_eqn_1_13b__y_a(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_3__T... ";
  try {
    auto r = VacuumTheory_eqn_1_3__T(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_3__k... ";
  try {
    auto r = VacuumTheory_eqn_1_3__k(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_3__m... ";
  try {
    auto r = VacuumTheory_eqn_1_3__m(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_3__v... ";
  try {
    auto r = VacuumTheory_eqn_1_3__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_7__R... ";
  try {
    auto r = VacuumTheory_eqn_1_7__R(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_7__T... ";
  try {
    auto r = VacuumTheory_eqn_1_7__T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_7__V... ";
  try {
    auto r = VacuumTheory_eqn_1_7__V(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_7__n... ";
  try {
    auto r = VacuumTheory_eqn_1_7__n(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_7__p... ";
  try {
    auto r = VacuumTheory_eqn_1_7__p(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_8__M... ";
  try {
    auto r = VacuumTheory_eqn_1_8__M(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_8__P... ";
  try {
    auto r = VacuumTheory_eqn_1_8__P(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_8__R... ";
  try {
    auto r = VacuumTheory_eqn_1_8__R(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_8__T... ";
  try {
    auto r = VacuumTheory_eqn_1_8__T(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_8__V... ";
  try {
    auto r = VacuumTheory_eqn_1_8__V(1.0, 1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_9__M... ";
  try {
    auto r = VacuumTheory_eqn_1_9__M(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_9__P... ";
  try {
    auto r = VacuumTheory_eqn_1_9__P(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_9__R... ";
  try {
    auto r = VacuumTheory_eqn_1_9__R(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_9__T... ";
  try {
    auto r = VacuumTheory_eqn_1_9__T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing VacuumTheory_eqn_1_9__rho... ";
  try {
    auto r = VacuumTheory_eqn_1_9__rho(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }

  std::cout << "\n" << pass << " passed, " << fail << " failed." << std::endl;
  return fail > 0 ? 1 : 0;
}