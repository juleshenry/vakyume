from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_1__D_r_cap import eqn_10_1__D_r
from .eqn_10_1__sig_R_cap import eqn_10_1__sig_R
from .eqn_10_1__w import eqn_10_1__w

class LiquidRing:
    eqn_10_1__D_r = eqn_10_1__D_r
    eqn_10_1__sig_R = eqn_10_1__sig_R
    eqn_10_1__w = eqn_10_1__w

    @kwasak_static
    def eqn_10_1(self, D_r=None, sig_R=None, w=None):
        return
