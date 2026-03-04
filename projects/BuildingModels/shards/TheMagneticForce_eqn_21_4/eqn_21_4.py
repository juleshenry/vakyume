from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_21_4__B_cap import eqn_21_4__B
from .eqn_21_4__T_cap import eqn_21_4__T
from .eqn_21_4__m import eqn_21_4__m
from .eqn_21_4__q import eqn_21_4__q
from .eqn_21_4__v import eqn_21_4__v


class TheMagneticForce:
    eqn_21_4__B = eqn_21_4__B
    eqn_21_4__T = eqn_21_4__T
    eqn_21_4__m = eqn_21_4__m
    eqn_21_4__q = eqn_21_4__q
    eqn_21_4__v = eqn_21_4__v

    @kwasak
    def eqn_21_4(self, B=None, T=None, m=None, q=None, v=None):
        """
        m := mass
        q := charge
        B := magnetic field
        T := period
        """
        return
