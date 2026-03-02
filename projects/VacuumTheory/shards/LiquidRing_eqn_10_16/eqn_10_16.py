from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_16__P_cap import eqn_10_16__P
from .eqn_10_16__S_0_cap import eqn_10_16__S_0
from .eqn_10_16__S_Th_cap import eqn_10_16__S_Th
from .eqn_10_16__p_0 import eqn_10_16__p_0

class LiquidRing:
    eqn_10_16__P = eqn_10_16__P
    eqn_10_16__S_0 = eqn_10_16__S_0
    eqn_10_16__S_Th = eqn_10_16__S_Th
    eqn_10_16__p_0 = eqn_10_16__p_0

    @kwasak_static
    def eqn_10_16(self, P=None, S_0=None, S_Th=None, p_0=None):
        return
