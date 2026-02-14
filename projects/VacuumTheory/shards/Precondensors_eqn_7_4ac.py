from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_4ac(P_c=None, n_i=None, n_nc=None, p=None, p_i=None, **kwargs):
        return

    @staticmethod
    def eqn_7_4ac__P_c(n_i: float, n_nc: float, p: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_4ac__n_i(P_c: float, n_nc: float, p: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_4ac__n_nc(P_c: float, n_i: float, p: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_4ac__p(P_c: float, n_i: float, n_nc: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result

    @staticmethod
    def eqn_7_4ac__p_i(P_c: float, n_i: float, n_nc: float, p: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result

