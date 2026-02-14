from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_19b(Re=None, h=None, mu=None, rho=None, v=None, w=None, **kwargs):
        return

    @staticmethod
    def eqn_2_19b__Re(h: float, mu: float, rho: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        Re = 2*h*rho*v*w/(mu*(h + w))
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_19b__h(Re: float, mu: float, rho: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        h = Re*mu*w/(-Re*mu + 2*rho*v*w)
        result.append(h)
        return result

    @staticmethod
    def eqn_2_19b__mu(Re: float, h: float, rho: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        mu = 2*h*rho*v*w/(Re*(h + w))
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_19b__rho(Re: float, h: float, mu: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        rho = Re*mu*(h + w)/(2*h*v*w)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_19b__v(Re: float, h: float, mu: float, rho: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        v = Re*mu*(h + w)/(2*h*rho*w)
        result.append(v)
        return result

    @staticmethod
    def eqn_2_19b__w(Re: float, h: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        w = Re*h*mu/(-Re*mu + 2*h*rho*v)
        result.append(w)
        return result

