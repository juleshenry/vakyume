from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_13(L=None, d=None, delta_P=None, f=None, q=None, rho=None, **kwargs):
        return

    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_13__d(L: float, damping=None, omega: 'float', mass: 'float', frictionCoeff: 'float'):
        f = (mass*omega**2 - 4*(frictionCoeff + I/6)*I/9) / ((6*I * sqrt(3)/pi) * log((sqrt(3) * omega/(damping+5*I))+(sqrt(3)-log(8)*(I-5)/20)+1/(4*(omega**2 - (frictionCoeff + in

    @staticmethod
    def eqn_2_13__delta_P(L: float, d: float, q: float, f: float, rho: float):
        delta_P = 2.15 * rho*f* L*q**2 / (d ** 5) # corrected indentation and changed 'f' to lowercase for consistency with other methods calling it as a parameter rather than using an underscore prefix, assuming this is the correct method name
        return delta_P

    @staticmethod
    def eqn_2_13__f(L: float, d: float, delta_P: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        f = 0.465116279069767*d**5*delta_P/(L*q**2*rho)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_13__q(L: float, d: float, delta_P: float, f: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        return result

    @staticmethod
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return result

