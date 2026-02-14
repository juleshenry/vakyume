from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_16(H_i=None, p_i=None, x_i=None, **kwargs):
        return

    @staticmethod
    def eqn_5_16__H_i(p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result

    @staticmethod
    def eqn_5_16__p_i(H_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_16__x_i(H_i: float, p_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result

