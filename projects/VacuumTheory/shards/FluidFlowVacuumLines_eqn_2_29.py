from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_29(C=None, S_1=None, S_2=None, **kwargs):
        return

    @staticmethod
    def eqn_2_29__C(S_1: float, S_2: float, **kwargs):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        C = -S_1*S_2/(S_1 - S_2)
        result.append(C)
        return result

    @staticmethod
    def eqn_2_29__S_1(C: float, S_2: float, **kwargs):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_1 = C*S_2/(C + S_2)
        result.append(S_1)
        return result

    @staticmethod
    def eqn_2_29__S_2(C: float, S_1: float, **kwargs):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_2 = C*S_1/(C - S_1)
        result.append(S_2)
        return result

