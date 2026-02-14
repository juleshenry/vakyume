from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_1(P=None, p_i=None, y_i=None, **kwargs):
        return

    @staticmethod
    def eqn_7_1__P(p_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_1__p_i(P: float, y_i: float, **kwargs):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_1__y_i(P: float, p_i: float, **kwargs):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return result

