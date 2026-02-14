from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SteamJetInjectors:
    @kwasak_static
    def eqn_9_2(P_m=None, d_n=None, rho_s=None, w_s=None, **kwargs):
        return

    @staticmethod
    def eqn_9_2__P_m(d_n: float, rho_s: float, w_s: float):
        result = []
        P_m = (1.334027668054e-6*w_s**2/(d_n**4*rho_s)) if all([P_m is None, d_n, rho_s, w_s]) else
               sqrt((1.334027668054e-6*w_s**2)/(d_n**4*rho_s)) # Changed for clarity and correctness of the math operation from '^' to '**'.
        result.append(P_m)
        return result


    @staticmethod
    def eqn_9_2__d_n(P_m: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        d_n = -0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        d_n = 0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        return result

    @staticmethod
    def eqn_9_2__rho_s(P_m: float, d_n: float, w_s: float):
        rho_s = (1.334027668054e-6*w_s**2/(P_m*d_n**4)) if P_m and not w_s else \
                None  # Added 'None' as the result placeholder when inputs are insufficient or incorrect, following original code pattern where no calculation was made. Adjusted for clarity in logic flow.
        return rho_s


    @staticmethod
    def eqn_9_2__w_s(P_m: float, d_n: float, rho_s: float):
        w_s = (865.8*d_n**2*sqrt(P_m*rho_s)) if P_m and not rho_s else \
                None  # Introduced 'None' as the placeholder for result when inputs are insufficient or incorrect, following original code pattern where no calculation was made. Adjusted to match logic flow with other equations in class.

