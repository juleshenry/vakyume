from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_19a(R_ll=None, Re=None, mu=None, rho=None, v=None, **kwargs):
        return

    @staticmethod
    def eqn_2_19a__R_ll(Re: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        R_ll = Re*mu/(4*rho*v)
        result.append(R_ll)
        return result

    @staticmethod
    def eqn_2_19a__Re(R_ll: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        Re = 4*R_ll*rho*v/mu
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_19a__mu(R_ll: float, Re: float, rho: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        mu = 4*R_ll*rho*v/Re
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_19a__rho(R_ll: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        rho = Re*mu/(4*R_ll*v)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_19a__v(R_ll: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        v = Re*mu/(4*R_ll*rho)
        result.append(v)
        return result

