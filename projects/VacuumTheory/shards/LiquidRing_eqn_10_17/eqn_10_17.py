from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_17__P_cap import eqn_10_17__P
from .eqn_10_17__S_0_cap import eqn_10_17__S_0
from .eqn_10_17__S_Th_cap import eqn_10_17__S_Th
from .eqn_10_17__p_0 import eqn_10_17__p_0
from .eqn_10_17__p_s import eqn_10_17__p_s

class LiquidRing:
    eqn_10_17__P = eqn_10_17__P
    eqn_10_17__S_0 = eqn_10_17__S_0
    eqn_10_17__S_Th = eqn_10_17__S_Th
    eqn_10_17__p_0 = eqn_10_17__p_0
    eqn_10_17__p_s = eqn_10_17__p_s

    @kwasak
    def eqn_10_17(self, P=None, S_0=None, S_Th=None, p_0=None, p_s=None):
        return
