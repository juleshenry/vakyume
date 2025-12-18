from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np
from kwasak import kwasak_static

class Precondensors:

    @kwasak_static
    def eqn_7_04ab(P_c: float = None, p: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ab__P_c(p: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ab__p(P_c: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ab__p_i(P_c: float, p: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = 0
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04ab__p_nc(P_c: float, p: float, p_i: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = -P_c + p
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04ac(P_c: float = None, n_i: float = None, n_nc: float = None, p: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ac__P_c(n_i: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ac__n_i(P_c: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04ac__n_nc(P_c: float, n_i: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04ac__p(P_c: float, n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ac__p_i(P_c: float, n_i: float, n_nc: float, p: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result