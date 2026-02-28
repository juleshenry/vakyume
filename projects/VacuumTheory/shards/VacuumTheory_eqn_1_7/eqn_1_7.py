from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_7__R_cap import eqn_1_7__R
from .eqn_1_7__T_cap import eqn_1_7__T
from .eqn_1_7__V_cap import eqn_1_7__V
from .eqn_1_7__n import eqn_1_7__n
from .eqn_1_7__p import eqn_1_7__p

class VacuumTheory:
    eqn_1_7__R = staticmethod(eqn_1_7__R)
    eqn_1_7__T = staticmethod(eqn_1_7__T)
    eqn_1_7__V = staticmethod(eqn_1_7__V)
    eqn_1_7__n = staticmethod(eqn_1_7__n)
    eqn_1_7__p = staticmethod(eqn_1_7__p)

    @kwasak_static
    def eqn_1_7(R=None, T=None, V=None, n=None, p=None, **kwargs):
        return
