from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_1(D=None, Re=None, mu=None, rho=None, v=None, **kwargs):
        return

    @staticmethod
    def eqn_2_1__D(Re: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_1__Re(D: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_1__mu(D: float, Re: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_1__rho(D: float, Re: float, mu: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_1__v(D: float, Re: float, mu: float, rho: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return result

