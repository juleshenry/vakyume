#include <iostream>

#include "applyingnewtonslaws.hpp"
#include "describingmotioninmultipledimensions.hpp"
#include "describingmotioninonedimension.hpp"
#include "electricchargesandfields.hpp"
#include "electriccircuits.hpp"
#include "electriccurrent.hpp"
#include "electromagneticinduction.hpp"
#include "fluidmechanics.hpp"
#include "gaussslaw.hpp"
#include "gravity.hpp"
#include "kinematics.hpp"
#include "linearmomentumandthecentreofmass.hpp"
#include "newtonslaws.hpp"
#include "potentialenergyandconservationofenergy.hpp"
#include "rotationaldynamics.hpp"
#include "rotationalenergyandmomentum.hpp"
#include "simpleharmonicmotion.hpp"
#include "themagneticforce.hpp"
#include "thetheoryofspecialrelativity.hpp"
#include "waves.hpp"
#include "workandenergy.hpp"

int main() {
  std::cout << "Running test suite..." << std::endl;
  int pass = 0, fail = 0;
  std::cout << "Testing ApplyingNewtonSLaws_eqn_6_1__F... ";
  try {
    auto r = ApplyingNewtonSLaws_eqn_6_1__F(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ApplyingNewtonSLaws_eqn_6_1__F_N... ";
  try {
    auto r = ApplyingNewtonSLaws_eqn_6_1__F_N(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ApplyingNewtonSLaws_eqn_6_2__a_2... ";
  try {
    auto r = ApplyingNewtonSLaws_eqn_6_2__a_2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ApplyingNewtonSLaws_eqn_6_2__m... ";
  try {
    auto r = ApplyingNewtonSLaws_eqn_6_2__m(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ApplyingNewtonSLaws_eqn_6_2__x... ";
  try {
    auto r = ApplyingNewtonSLaws_eqn_6_2__x(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_1__t... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_1__t(1.0, 1.0, 1.0, 1.0,
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
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_1__v... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_1__v(1.0, 1.0, 1.0, 1.0,
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
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_1__x_1... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_1__x_1(1.0, 1.0, 1.0,
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
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_1__x_2... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_1__x_2(1.0, 1.0, 1.0,
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
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_1__y_1... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_1__y_1(1.0, 1.0, 1.0,
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
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_1__y_2... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_1__y_2(1.0, 1.0, 1.0,
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
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_3__v... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_3__v(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_3__vx... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_3__vx(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_3__vy... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_3__vy(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_3__x... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_3__x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_3__y... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_3__y(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_4__a... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_4__a(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_4__ax... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_4__ax(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_4__ay... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_4__ay(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_4__x... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_4__x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_4__y... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_4__y(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_5__a... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_5__a(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_5__t... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_5__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_5__v... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_5__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_5__v_0... ";
  try {
    auto r = DescribingMotionInMultipleDimensions_eqn_4_5__v_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_6__a... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_6__a(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_6__r... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_6__r(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_6__r_0... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_6__r_0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_6__t... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_6__t(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInMultipleDimensions_eqn_4_6__v_0... ";
  try {
    auto r =
        DescribingMotionInMultipleDimensions_eqn_4_6__v_0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_1__t... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_1__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_1__v_x... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_1__v_x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_1__x... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_1__x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_1__x_0... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_1__x_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_14__a... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_14__a(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_14__a_A... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_14__a_A(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_15__t... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_15__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_15__v... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_15__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_15__v_A... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_15__v_A(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_15__v_B... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_15__v_B(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_2__ax... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_2__ax(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_2__t... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_2__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_2__v... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_2__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_2__v_0x... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_2__v_0x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_3__ax... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_3__ax(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_3__t... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_3__t(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_3__v_0x... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_3__v_0x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_3__x... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_3__x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_3__x_0... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_3__x_0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_4__a... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_4__a(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_4__t... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_4__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_4__v... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_4__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_4__v_0... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_4__v_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_6__t... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_6__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_6__v_0... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_6__v_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_6__x... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_6__x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing DescribingMotionInOneDimension_eqn_3_6__x_0... ";
  try {
    auto r = DescribingMotionInOneDimension_eqn_3_6__x_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricChargesAndFields_eqn_16_3__F_g... ";
  try {
    auto r = ElectricChargesAndFields_eqn_16_3__F_g(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricChargesAndFields_eqn_16_3__G... ";
  try {
    auto r = ElectricChargesAndFields_eqn_16_3__G(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricChargesAndFields_eqn_16_3__m_e... ";
  try {
    auto r = ElectricChargesAndFields_eqn_16_3__m_e(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricChargesAndFields_eqn_16_3__m_p... ";
  try {
    auto r = ElectricChargesAndFields_eqn_16_3__m_p(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricChargesAndFields_eqn_16_3__r... ";
  try {
    auto r = ElectricChargesAndFields_eqn_16_3__r(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_3__I1... ";
  try {
    auto r = ElectricCircuits_eqn_20_3__I1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_3__I2... ";
  try {
    auto r = ElectricCircuits_eqn_20_3__I2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_3__I3... ";
  try {
    auto r = ElectricCircuits_eqn_20_3__I3(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_4__I... ";
  try {
    auto r = ElectricCircuits_eqn_20_4__I(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_4__Reff... ";
  try {
    auto r = ElectricCircuits_eqn_20_4__Reff(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_4__Vvoltmeter... ";
  try {
    auto r = ElectricCircuits_eqn_20_4__Vvoltmeter(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_5__C... ";
  try {
    auto r = ElectricCircuits_eqn_20_5__C(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_5__IR... ";
  try {
    auto r = ElectricCircuits_eqn_20_5__IR(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_5__Q... ";
  try {
    auto r = ElectricCircuits_eqn_20_5__Q(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCircuits_eqn_20_5__V... ";
  try {
    auto r = ElectricCircuits_eqn_20_5__V(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCurrent_eqn_19_1__I... ";
  try {
    auto r = ElectricCurrent_eqn_19_1__I(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCurrent_eqn_19_1__Q... ";
  try {
    auto r = ElectricCurrent_eqn_19_1__Q(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCurrent_eqn_19_1__t... ";
  try {
    auto r = ElectricCurrent_eqn_19_1__t(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCurrent_eqn_19_6__R1... ";
  try {
    auto r = ElectricCurrent_eqn_19_6__R1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCurrent_eqn_19_6__R2... ";
  try {
    auto r = ElectricCurrent_eqn_19_6__R2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectricCurrent_eqn_19_6__Reff... ";
  try {
    auto r = ElectricCurrent_eqn_19_6__Reff(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectromagneticInduction_eqn_23_1__B0... ";
  try {
    auto r = ElectromagneticInduction_eqn_23_1__B0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectromagneticInduction_eqn_23_1__E... ";
  try {
    auto r = ElectromagneticInduction_eqn_23_1__E(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectromagneticInduction_eqn_23_1__R... ";
  try {
    auto r = ElectromagneticInduction_eqn_23_1__R(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectromagneticInduction_eqn_23_1__V... ";
  try {
    auto r = ElectromagneticInduction_eqn_23_1__V(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectromagneticInduction_eqn_23_1__a... ";
  try {
    auto r = ElectromagneticInduction_eqn_23_1__a(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectromagneticInduction_eqn_23_1__dt... ";
  try {
    auto r = ElectromagneticInduction_eqn_23_1__dt(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing ElectromagneticInduction_eqn_23_1__r... ";
  try {
    auto r = ElectromagneticInduction_eqn_23_1__r(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidMechanics_eqn_15_8__R... ";
  try {
    auto r = FluidMechanics_eqn_15_8__R(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidMechanics_eqn_15_8__R1... ";
  try {
    auto r = FluidMechanics_eqn_15_8__R1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing FluidMechanics_eqn_15_8__R2... ";
  try {
    auto r = FluidMechanics_eqn_15_8__R2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing GaussSLaw_eqn_17_2__E... ";
  try {
    auto r = GaussSLaw_eqn_17_2__E(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing GaussSLaw_eqn_17_2__R... ";
  try {
    auto r = GaussSLaw_eqn_17_2__R(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing GaussSLaw_eqn_17_3__E... ";
  try {
    auto r = GaussSLaw_eqn_17_3__E(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing GaussSLaw_eqn_17_3__Q... ";
  try {
    auto r = GaussSLaw_eqn_17_3__Q(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing GaussSLaw_eqn_17_3__R... ";
  try {
    auto r = GaussSLaw_eqn_17_3__R(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_2__G... ";
  try {
    auto r = Gravity_eqn_9_2__G(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_2__M... ";
  try {
    auto r = Gravity_eqn_9_2__M(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_2__R... ";
  try {
    auto r = Gravity_eqn_9_2__R(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_2__T... ";
  try {
    auto r = Gravity_eqn_9_2__T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_2__pi... ";
  try {
    auto r = Gravity_eqn_9_2__pi(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_3__G... ";
  try {
    auto r = Gravity_eqn_9_3__G(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_3__M... ";
  try {
    auto r = Gravity_eqn_9_3__M(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_3__R... ";
  try {
    auto r = Gravity_eqn_9_3__R(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_3__a... ";
  try {
    auto r = Gravity_eqn_9_3__a(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_3__m... ";
  try {
    auto r = Gravity_eqn_9_3__m(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_4__G... ";
  try {
    auto r = Gravity_eqn_9_4__G(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_4__M... ";
  try {
    auto r = Gravity_eqn_9_4__M(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_4__U... ";
  try {
    auto r = Gravity_eqn_9_4__U(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_4__m... ";
  try {
    auto r = Gravity_eqn_9_4__m(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Gravity_eqn_9_4__r... ";
  try {
    auto r = Gravity_eqn_9_4__r(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_1__a... ";
  try {
    auto r = Kinematics_eqn_1_1__a(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_1__t... ";
  try {
    auto r = Kinematics_eqn_1_1__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_1__v... ";
  try {
    auto r = Kinematics_eqn_1_1__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_1__v0... ";
  try {
    auto r = Kinematics_eqn_1_1__v0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_2__a... ";
  try {
    auto r = Kinematics_eqn_1_2__a(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_2__t... ";
  try {
    auto r = Kinematics_eqn_1_2__t(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_2__v0... ";
  try {
    auto r = Kinematics_eqn_1_2__v0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_2__x... ";
  try {
    auto r = Kinematics_eqn_1_2__x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_2__x0... ";
  try {
    auto r = Kinematics_eqn_1_2__x0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_3__a... ";
  try {
    auto r = Kinematics_eqn_1_3__a(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_3__dx... ";
  try {
    auto r = Kinematics_eqn_1_3__dx(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_3__v... ";
  try {
    auto r = Kinematics_eqn_1_3__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Kinematics_eqn_1_3__v0... ";
  try {
    auto r = Kinematics_eqn_1_3__v0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_1__m... ";
  try {
    auto r = LinearMomentumAndTheCentreOfMass_eqn_10_1__m(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_1__p... ";
  try {
    auto r = LinearMomentumAndTheCentreOfMass_eqn_10_1__p(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_1__v... ";
  try {
    auto r = LinearMomentumAndTheCentreOfMass_eqn_10_1__v(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_11__F_ext... ";
  try {
    auto r = LinearMomentumAndTheCentreOfMass_eqn_10_11__F_ext(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_11__M... ";
  try {
    auto r = LinearMomentumAndTheCentreOfMass_eqn_10_11__M(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing LinearMomentumAndTheCentreOfMass_eqn_10_11__a_CM... ";
  try {
    auto r = LinearMomentumAndTheCentreOfMass_eqn_10_11__a_CM(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing NewtonSLaws_eqn_5_1__F... ";
  try {
    auto r = NewtonSLaws_eqn_5_1__F(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing NewtonSLaws_eqn_5_1__a... ";
  try {
    auto r = NewtonSLaws_eqn_5_1__a(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing NewtonSLaws_eqn_5_1__m... ";
  try {
    auto r = NewtonSLaws_eqn_5_1__m(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing NewtonSLaws_eqn_5_2__F... ";
  try {
    auto r = NewtonSLaws_eqn_5_2__F(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing NewtonSLaws_eqn_5_2__F_action... ";
  try {
    auto r = NewtonSLaws_eqn_5_2__F_action(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing NewtonSLaws_eqn_5_2__F_reaction... ";
  try {
    auto r = NewtonSLaws_eqn_5_2__F_reaction(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing NewtonSLaws_eqn_5_2__x... ";
  try {
    auto r = NewtonSLaws_eqn_5_2__x(1.0);
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
      << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_6__K_A... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_6__K_A(1.0, 1.0);
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
      << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_6__K_B... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_6__K_B(1.0, 1.0);
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
      << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_6__W_net... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_6__W_net(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__A... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__A(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__B... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__B(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__K... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__K(1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__L... ";
  try {
    auto r =
        PotentialEnergyAndConservationOfEnergy_eqn_8_9__L(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__U... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__U(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__W... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__W(1.0, 1.0, 1.0);
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
      << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_C... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_C(1.0, 1.0);
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
      << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_NC... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_NC(1.0, 1.0);
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
      << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_net... ";
  try {
    auto r = PotentialEnergyAndConservationOfEnergy_eqn_8_9__W_net(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__g... ";
  try {
    auto r =
        PotentialEnergyAndConservationOfEnergy_eqn_8_9__g(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__m... ";
  try {
    auto r =
        PotentialEnergyAndConservationOfEnergy_eqn_8_9__m(1.0, 1.0, 1.0, 1.0);
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
      << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__v_x... ";
  try {
    auto r =
        PotentialEnergyAndConservationOfEnergy_eqn_8_9__v_x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing PotentialEnergyAndConservationOfEnergy_eqn_8_9__x... ";
  try {
    auto r =
        PotentialEnergyAndConservationOfEnergy_eqn_8_9__x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_10__ICM... ";
  try {
    auto r = RotationalDynamics_eqn_11_10__ICM(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_10__Ih... ";
  try {
    auto r = RotationalDynamics_eqn_11_10__Ih(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_10__M... ";
  try {
    auto r = RotationalDynamics_eqn_11_10__M(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_10__h... ";
  try {
    auto r = RotationalDynamics_eqn_11_10__h(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_8__I... ";
  try {
    auto r = RotationalDynamics_eqn_11_8__I(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_8__i... ";
  try {
    auto r = RotationalDynamics_eqn_11_8__i(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_8__m... ";
  try {
    auto r = RotationalDynamics_eqn_11_8__m(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalDynamics_eqn_11_8__r... ";
  try {
    auto r = RotationalDynamics_eqn_11_8__r(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalEnergyAndMomentum_eqn_12_11__I... ";
  try {
    auto r = RotationalEnergyAndMomentum_eqn_12_11__I(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalEnergyAndMomentum_eqn_12_11__L... ";
  try {
    auto r = RotationalEnergyAndMomentum_eqn_12_11__L(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing RotationalEnergyAndMomentum_eqn_12_11__v... ";
  try {
    auto r = RotationalEnergyAndMomentum_eqn_12_11__v(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_1__E... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_1__E(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_1__m... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_1__m(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_1__v... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_1__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_1__x... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_1__x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_2__d2x... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_2__d2x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_2__dt2... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_2__dt2(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_2__m... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_2__m(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_2__x... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_2__x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_2__x0... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_2__x0(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_2__x1... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_2__x1(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing SimpleHarmonicMotion_eqn_13_2__x2... ";
  try {
    auto r = SimpleHarmonicMotion_eqn_13_2__x2(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheMagneticForce_eqn_21_4__B... ";
  try {
    auto r = TheMagneticForce_eqn_21_4__B(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheMagneticForce_eqn_21_4__T... ";
  try {
    auto r = TheMagneticForce_eqn_21_4__T(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheMagneticForce_eqn_21_4__m... ";
  try {
    auto r = TheMagneticForce_eqn_21_4__m(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheMagneticForce_eqn_21_4__q... ";
  try {
    auto r = TheMagneticForce_eqn_21_4__q(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheMagneticForce_eqn_21_4__v... ";
  try {
    auto r = TheMagneticForce_eqn_21_4__v(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_11__E... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_11__E(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_11__c... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_11__c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_11__m_0... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_11__m_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_11__p... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_11__p(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__FE... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__FE(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__c... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__l... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__l(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__r... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__r(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__t... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__t(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__u_0... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__u_0(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__u_x... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__u_x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__v... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__vx... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__vx(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_2__x... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_2__x(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_6__c... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_6__c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_6__u... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_6__u(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_6__ux... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_6__ux(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_6__v... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_6__v(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_9__K... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_9__K(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_9__c... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_9__c(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_9__m_0... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_9__m_0(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing TheTheoryOfSpecialRelativity_eqn_24_9__u... ";
  try {
    auto r = TheTheoryOfSpecialRelativity_eqn_24_9__u(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Waves_eqn_14_11__fn... ";
  try {
    auto r = Waves_eqn_14_11__fn(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Waves_eqn_14_11__n... ";
  try {
    auto r = Waves_eqn_14_11__n(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing Waves_eqn_14_11__v... ";
  try {
    auto r = Waves_eqn_14_11__v(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_1__F... ";
  try {
    auto r = WorkAndEnergy_eqn_7_1__F(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_1__W... ";
  try {
    auto r = WorkAndEnergy_eqn_7_1__W(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_1__d... ";
  try {
    auto r = WorkAndEnergy_eqn_7_1__d(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_2__F1... ";
  try {
    auto r = WorkAndEnergy_eqn_7_2__F1(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_2__F2... ";
  try {
    auto r = WorkAndEnergy_eqn_7_2__F2(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_2__F3... ";
  try {
    auto r = WorkAndEnergy_eqn_7_2__F3(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_2__W_tot... ";
  try {
    auto r = WorkAndEnergy_eqn_7_2__W_tot(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_2__x... ";
  try {
    auto r = WorkAndEnergy_eqn_7_2__x(1.0, 1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_6__F_... ";
  try {
    auto r = WorkAndEnergy_eqn_7_6__F_(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_6__W_f... ";
  try {
    auto r = WorkAndEnergy_eqn_7_6__W_f(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_6__d... ";
  try {
    auto r = WorkAndEnergy_eqn_7_6__d(1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_7__W_F... ";
  try {
    auto r = WorkAndEnergy_eqn_7_7__W_F(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_7__W_f... ";
  try {
    auto r = WorkAndEnergy_eqn_7_7__W_f(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_7__W_g... ";
  try {
    auto r = WorkAndEnergy_eqn_7_7__W_g(1.0, 1.0, 1.0);
    std::cout << "OK (" << r.size() << " results)" << std::endl;
    pass++;
  } catch (std::exception &e) {
    std::cout << "FAIL: " << e.what() << std::endl;
    fail++;
  } catch (...) {
    std::cout << "FAIL" << std::endl;
    fail++;
  }
  std::cout << "Testing WorkAndEnergy_eqn_7_7__W_net... ";
  try {
    auto r = WorkAndEnergy_eqn_7_7__W_net(1.0, 1.0, 1.0);
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