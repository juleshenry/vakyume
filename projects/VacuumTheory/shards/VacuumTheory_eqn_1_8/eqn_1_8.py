from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_8__M_cap import eqn_1_8__M
from .eqn_1_8__P_cap import eqn_1_8__P
from .eqn_1_8__R_cap import eqn_1_8__R
from .eqn_1_8__T_cap import eqn_1_8__T
from .eqn_1_8__V_cap import eqn_1_8__V
from .eqn_1_8__m import eqn_1_8__m

class VacuumTheory:
    eqn_1_8__M = staticmethod(eqn_1_8__M)
    eqn_1_8__P = staticmethod(eqn_1_8__P)
    eqn_1_8__R = staticmethod(eqn_1_8__R)
    eqn_1_8__T = staticmethod(eqn_1_8__T)
    eqn_1_8__V = staticmethod(eqn_1_8__V)
    eqn_1_8__m = staticmethod(eqn_1_8__m)

    @kwasak_static
    def eqn_1_8(M=None, P=None, R=None, T=None, V=None, m=None, **kwargs):
        return
