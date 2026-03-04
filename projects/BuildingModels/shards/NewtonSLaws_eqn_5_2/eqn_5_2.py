from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_5_2__F_cap import eqn_5_2__F
from .eqn_5_2__x import eqn_5_2__x


class NewtonSLaws:
    eqn_5_2__F = eqn_5_2__F
    eqn_5_2__x = eqn_5_2__x

    @kwasak
    def eqn_5_2(self, F=None, x=None):
        """
        F := force
        F_action := action force
        F_reaction := reaction force
        F := force
        k := spring constant
        x := displacement from rest length
        """
        return
