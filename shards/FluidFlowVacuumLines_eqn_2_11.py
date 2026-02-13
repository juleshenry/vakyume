from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_11(C_L=None, C_T=None, F_p=None, **kwargs):
        return

    @staticmethod
    def eqn_2_11__C_L(C_T: float, F_p: float):
        # [.pyeqn] C_T = C_L * F_p
        result = []
        C_L = C_T/F_p
        result.append(C_L)
        return result

    @staticmethod
    def eqn_2_11__C_T(C_L: float, F_p: float):
        # [.pyeqn] C_T = C_L * F_p
        result = []
        C_T = C_L*F_p
        result.append(C_T)
        return result

    @staticmethod
    def eqn_2_11__F_p(C_L: float, C_T: float):
        # [.pyeqn] C_T = C_L * F_p
        result = []
        F_p = C_T/C_L
        result.append(F_p)
        return result

