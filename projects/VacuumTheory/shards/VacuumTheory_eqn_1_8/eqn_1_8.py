from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_1_8__M_cap import eqn_1_8__M
from .eqn_1_8__P_cap import eqn_1_8__P
from .eqn_1_8__R_cap import eqn_1_8__R
from .eqn_1_8__T_cap import eqn_1_8__T
from .eqn_1_8__V_cap import eqn_1_8__V
from .eqn_1_8__m import eqn_1_8__m


class VacuumTheory:
    eqn_1_8__M = eqn_1_8__M
    eqn_1_8__P = eqn_1_8__P
    eqn_1_8__R = eqn_1_8__R
    eqn_1_8__T = eqn_1_8__T
    eqn_1_8__V = eqn_1_8__V
    eqn_1_8__m = eqn_1_8__m

    @kwasak
    def eqn_1_8(self, M=None, P=None, R=None, T=None, V=None, m=None):
        return
