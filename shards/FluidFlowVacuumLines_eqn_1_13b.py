from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_1_13b(P: float = None, p_a: float = None, y_a: float = None, **kwargs):
        return

    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a / y_a
        result.append(P)
        return result

    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P * y_a
        result.append(p_a)
        return result

    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        y_a = p_a / P
        result.append(y_a)
        return result




