from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_17(L=None, d=None, delta_P=None, mu=None, v=None, **kwargs):
        return

    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, v: float, **kwargs):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, v: float, **kwargs):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        d = -0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        d = 0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, v: float, **kwargs):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        delta_P = 0.0345*L*mu*v/d**2
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, v: float, **kwargs):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        mu = 28.9855072463768*d**2*delta_P/(L*v)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__v(L: float, d: float, delta_P: float, mu: float, **kwargs):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
        return result

