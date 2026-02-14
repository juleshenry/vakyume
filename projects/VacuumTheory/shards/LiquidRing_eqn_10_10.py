from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_10(bhp=None, bhp_0=None, mu=None, rho=None, **kwargs):
        return

    @staticmethod
    def eqn_10_10__bhp(bhp_0: float, mu: float, rho: float, **kwargs):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp = 0.0005*bhp_0*(31.0*mu**(4/25)*rho**(21/25) + 1000.0)
        result.append(bhp)
        return result

    @staticmethod
    def eqn_10_10__bhp_0(bhp: float, mu: float, rho: float, **kwargs):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp_0 = 2000.0*bhp/(31.0*mu**0.16*rho**0.84 + 1000.0)
        result.append(bhp_0)
        return result

    @staticmethod
    def eqn_10_10__mu(bhp: float, bhp_0: float, rho: float, **kwargs):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        mu = -4.7751763343393e-10*I*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        mu = 4.7751763343393e-10*I*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        mu = -4.7751763343393e-10*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        mu = 4.7751763343393e-10*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        return result

    @staticmethod
    def eqn_10_10__rho(bhp: float, bhp_0: float, mu: float, **kwargs):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        # Error during Sympy solve: Sympy solve failed
        def func(rho):
            # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
            return eval("(bhp_0 * (0.5 + 0.0155 * x ** 0.84 * mu ** 0.16)) - (bhp)".replace('x', str(rho)))
        raise UnsolvedException("Pending LLM/Manual Repair")

