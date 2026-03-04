from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_21__P import eqn_10_21__P
from .eqn_10_21__P_d import eqn_10_21__P_d
from .eqn_10_21__P_prime import eqn_10_21__P_prime

class LiquidRing:
    eqn_10_21__P = eqn_10_21__P
    eqn_10_21__P_d = eqn_10_21__P_d
    eqn_10_21__P_prime = eqn_10_21__P_prime

    @kwasak
    def eqn_10_21(self, P=None, P_d=None, P_prime=None):
        """
        P_prime := pseudo suction pressure
        P_d := actual pump discharge pressure
        """
        return
