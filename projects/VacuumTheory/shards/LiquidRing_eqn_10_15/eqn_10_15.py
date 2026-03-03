from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_15__P_cap import eqn_10_15__P
from .eqn_10_15__S_cap_T_caph import eqn_10_15__S_Th
from .eqn_10_15__S_cap_p import eqn_10_15__S_p
from .eqn_10_15__p_s import eqn_10_15__p_s


class LiquidRing:
    eqn_10_15__P = eqn_10_15__P
    eqn_10_15__S_Th = eqn_10_15__S_Th
    eqn_10_15__S_p = eqn_10_15__S_p
    eqn_10_15__p_s = eqn_10_15__p_s

    @kwasak
    def eqn_10_15(self, P=None, S_Th=None, S_p=None, p_s=None):
        return
