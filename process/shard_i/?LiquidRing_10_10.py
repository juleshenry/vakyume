from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np
from kwasak import kwasak_static

class LiquidRing:

    @kwasak_static
    def eqn_10_10(bhp: float = None, bhp_0: float = None, mu: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_10__bhp(bhp_0: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp = 0.0005*bhp_0*(31.0*mu**(4/25)*rho**(21/25) + 1000.0)
        result.append(bhp)
        return result

    @staticmethod
    def eqn_10_10__bhp_0(bhp: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp_0 = 2000.0*bhp/(31.0*mu**0.16*rho**0.84 + 1000.0)
        result.append(bhp_0)
        return result

    @staticmethod
    def eqn_10_10__mu(bhp: float, bhp_0: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        mu = -204374584201.104*I*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = 204374584201.104*I*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = -204374584201.104*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = 204374584201.104*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        return result

    @staticmethod
    def eqn_10_10__rho(bhp: float, bhp_0: float, mu: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        # [Sympy Failover]
        pass # Ollama offline
