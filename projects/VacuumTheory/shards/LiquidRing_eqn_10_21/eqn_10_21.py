from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_21__P_cap import eqn_10_21__P
from .eqn_10_21__P_d_cap import eqn_10_21__P_d
from .eqn_10_21__P_prime_cap import eqn_10_21__P_prime

class LiquidRing:
    eqn_10_21__P = eqn_10_21__P
    eqn_10_21__P_d = eqn_10_21__P_d
    eqn_10_21__P_prime = eqn_10_21__P_prime

    @kwasak_static
    def eqn_10_21(self, P=None, P_d=None, P_prime=None):
        return
