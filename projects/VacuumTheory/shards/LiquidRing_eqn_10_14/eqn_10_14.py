from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_14__T_cap_c import eqn_10_14__T_c
from .eqn_10_14__T_cap_s import eqn_10_14__T_s

class LiquidRing:
    eqn_10_14__T_c = eqn_10_14__T_c
    eqn_10_14__T_s = eqn_10_14__T_s

    @kwasak
    def eqn_10_14(self, T_c=None, T_s=None):
        return
