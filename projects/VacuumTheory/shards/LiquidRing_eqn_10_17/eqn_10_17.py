from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_10_17__P_cap import eqn_10_17__P
from .eqn_10_17__S_0_cap import eqn_10_17__S_0
from .eqn_10_17__S_Th_cap import eqn_10_17__S_Th
from .eqn_10_17__p_0 import eqn_10_17__p_0
from .eqn_10_17__p_s import eqn_10_17__p_s

class LiquidRing:
    eqn_10_17__P = staticmethod(eqn_10_17__P)
    eqn_10_17__S_0 = staticmethod(eqn_10_17__S_0)
    eqn_10_17__S_Th = staticmethod(eqn_10_17__S_Th)
    eqn_10_17__p_0 = staticmethod(eqn_10_17__p_0)
    eqn_10_17__p_s = staticmethod(eqn_10_17__p_s)

    @kwasak_static
    def eqn_10_17(P=None, S_0=None, S_Th=None, p_0=None, p_s=None, **kwargs):
        return
