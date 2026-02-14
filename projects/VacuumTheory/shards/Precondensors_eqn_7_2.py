from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_2(P_i_0=None, p_i=None, x_i=None, **kwargs):
        return

    @staticmethod
    def eqn_7_2__P_i_0(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i/x_i
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_2__p_i(P_i_0: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        p_i = P_i_0*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_2__x_i(P_i_0: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        x_i = p_i/P_i_0
        result.append(x_i)
        return result

