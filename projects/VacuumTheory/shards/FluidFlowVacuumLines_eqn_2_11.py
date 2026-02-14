from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_11(D=None, L=None, f=None, g_c=None, h_r=None, v=None, **kwargs):
        return

    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float, **kwargs):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float, **kwargs):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float, **kwargs):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_11__g_c(D: float, L: float, f: float, h_r: float, v: float, **kwargs):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        g_c = L*f*v**2/(2*D*h_r)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float, **kwargs):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        h_r = L*f*v**2/(2*D*g_c)
        result.append(h_r)
        return result

    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float, **kwargs):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return result

