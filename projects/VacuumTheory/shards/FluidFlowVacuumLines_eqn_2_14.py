from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_14(M=None, R=None, T=None, g_c=None, k=None, v_s=None, **kwargs):
        return

    @staticmethod
    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return result

    @staticmethod
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return result

    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return result

    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return result

    @staticmethod
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return result

