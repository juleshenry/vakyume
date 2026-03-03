from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_16__P_cap import eqn_10_16__P
from .eqn_10_16__S_cap_0 import eqn_10_16__S_0
from .eqn_10_16__S_cap_T_caph import eqn_10_16__S_Th
from .eqn_10_16__p_0 import eqn_10_16__p_0


class LiquidRing:
    eqn_10_16__P = eqn_10_16__P
    eqn_10_16__S_0 = eqn_10_16__S_0
    eqn_10_16__S_Th = eqn_10_16__S_Th
    eqn_10_16__p_0 = eqn_10_16__p_0

    @kwasak
    def eqn_10_16(self, P=None, S_0=None, S_Th=None, p_0=None):
        return
