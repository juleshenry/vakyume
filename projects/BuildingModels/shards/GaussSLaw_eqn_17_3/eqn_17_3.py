from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_17_3__E_cap import eqn_17_3__E
from .eqn_17_3__Q_cap import eqn_17_3__Q
from .eqn_17_3__R_cap import eqn_17_3__R


class GaussSLaw:
    eqn_17_3__E = eqn_17_3__E
    eqn_17_3__Q = eqn_17_3__Q
    eqn_17_3__R = eqn_17_3__R

    @kwasak
    def eqn_17_3(self, E=None, Q=None, R=None):
        """
        Q := charge
        R := radius
        """
        return
