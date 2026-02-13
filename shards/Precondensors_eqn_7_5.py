from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_5(N_i=None, N_nc=None, P=None, P_c=None, p_i=None, **kwargs):
        return

    @staticmethod
    def eqn_7_5__N_i(N_nc: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
        return result

    @staticmethod
    def eqn_7_5__N_nc(N_i: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return result

    @staticmethod
    def eqn_7_5__P(N_i: float, N_nc: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc*p_i/N_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_5__P_c(N_i: float, N_nc: float, P: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc*p_i/N_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_5__p_i(N_i: float, N_nc: float, P: float, P_c: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i*(P - P_c)/N_nc
        result.append(p_i)
        return result

