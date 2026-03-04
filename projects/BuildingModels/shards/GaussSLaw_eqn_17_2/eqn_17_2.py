from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_17_2__E_cap import eqn_17_2__E
from .eqn_17_2__R_cap import eqn_17_2__R


class GaussSLaw:
    eqn_17_2__E = eqn_17_2__E
    eqn_17_2__R = eqn_17_2__R

    @kwasak
    def eqn_17_2(self, E=None, R=None):
        """
        E := electric field
        R := radius
        """
        return
