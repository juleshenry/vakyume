from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_12__T_c_cap import eqn_10_12__T_c
from .eqn_10_12__T_s_cap import eqn_10_12__T_s

class LiquidRing:
    eqn_10_12__T_c = eqn_10_12__T_c
    eqn_10_12__T_s = eqn_10_12__T_s

    @kwasak_static
    def eqn_10_12(self, T_c=None, T_s=None):
        return
