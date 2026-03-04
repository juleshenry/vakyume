from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_20__P_cap import eqn_10_20__P
from .eqn_10_20__S_0_cap import eqn_10_20__S_0
from .eqn_10_20__S_p_cap import eqn_10_20__S_p
from .eqn_10_20__T_e_cap import eqn_10_20__T_e
from .eqn_10_20__T_i_cap import eqn_10_20__T_i
from .eqn_10_20__p_0 import eqn_10_20__p_0
from .eqn_10_20__p_c import eqn_10_20__p_c
from .eqn_10_20__p_s import eqn_10_20__p_s

class LiquidRing:
    eqn_10_20__P = eqn_10_20__P
    eqn_10_20__S_0 = eqn_10_20__S_0
    eqn_10_20__S_p = eqn_10_20__S_p
    eqn_10_20__T_e = eqn_10_20__T_e
    eqn_10_20__T_i = eqn_10_20__T_i
    eqn_10_20__p_0 = eqn_10_20__p_0
    eqn_10_20__p_c = eqn_10_20__p_c
    eqn_10_20__p_s = eqn_10_20__p_s

    @kwasak
    def eqn_10_20(self, P=None, S_0=None, S_p=None, T_e=None, T_i=None, p_0=None, p_c=None, p_s=None):
        return
