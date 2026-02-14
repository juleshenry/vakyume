from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_4a(P=None, p_c=None, p_nc=None, **kwargs):
        return

    @staticmethod
    def eqn_7_4a__P(p_c: float, p_nc: float, **kwargs):
        # [.pyeqn] p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return result

    @staticmethod
    def eqn_7_4a__p_c(P: float, p_nc: float, **kwargs):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_4a__p_nc(P: float, p_c: float, **kwargs):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return result

