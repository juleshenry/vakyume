from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_36(C=None, C_0=None, F_t=None, **kwargs):
        return

    @staticmethod
    def eqn_2_36__C(C_0: float, F_t: float):
        # [.pyeqn] C = C_0 * F_t
        result = []
        C = C_0*F_t
        result.append(C)
        return result

    @staticmethod
    def eqn_2_36__C_0(C: float, F_t: float):
        # [.pyeqn] C = C_0 * F_t
        result = []
        C_0 = C/F_t
        result.append(C_0)
        return result

    @staticmethod
    def eqn_2_36__F_t(C: float, C_0: float):
        # [.pyeqn] C = C_0 * F_t
        result = []
        F_t = C/C_0
        result.append(F_t)
        return result

