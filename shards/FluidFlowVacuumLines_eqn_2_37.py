from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_37(A=None, C=None, F_t=None, M=None, T=None, **kwargs):
        return

    @staticmethod
    def eqn_2_37__A(C: float, F_t: float, M: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        A = 0.000681714375311032*C**2*M/(F_t*T)
        result.append(A)
        return result

    @staticmethod
    def eqn_2_37__C(A: float, F_t: float, M: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        C = 38.3*sqrt(A*F_t*T/M)
        result.append(C)
        return result

    @staticmethod
    def eqn_2_37__F_t(A: float, C: float, M: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        F_t = 0.000681714375311032*C**2*M/(A*T)
        result.append(F_t)
        return result

    @staticmethod
    def eqn_2_37__M(A: float, C: float, F_t: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        M = 1466.89*A*F_t*T/C**2
        result.append(M)
        return result

    @staticmethod
    def eqn_2_37__T(A: float, C: float, F_t: float, M: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        T = 0.000681714375311032*C**2*M/(A*F_t)
        result.append(T)
        return result

