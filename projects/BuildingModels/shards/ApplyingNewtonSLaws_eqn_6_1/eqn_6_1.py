from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_6_1__F_cap import eqn_6_1__F
from .eqn_6_1__F_cap_N_cap import eqn_6_1__F_N


class ApplyingNewtonSLaws:
    eqn_6_1__F = eqn_6_1__F
    eqn_6_1__F_N = eqn_6_1__F_N

    @kwasak
    def eqn_6_1(self, F=None, F_N=None):
        """
        F := force
        F_N := normal force
        """
        return
