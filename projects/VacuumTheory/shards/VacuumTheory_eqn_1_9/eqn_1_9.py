from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_9__M_cap import eqn_1_9__M
from .eqn_1_9__P_cap import eqn_1_9__P
from .eqn_1_9__R_cap import eqn_1_9__R
from .eqn_1_9__T_cap import eqn_1_9__T
from .eqn_1_9__rho import eqn_1_9__rho

class VacuumTheory:
    eqn_1_9__M = staticmethod(eqn_1_9__M)
    eqn_1_9__P = staticmethod(eqn_1_9__P)
    eqn_1_9__R = staticmethod(eqn_1_9__R)
    eqn_1_9__T = staticmethod(eqn_1_9__T)
    eqn_1_9__rho = staticmethod(eqn_1_9__rho)

    @kwasak_static
    def eqn_1_9(M=None, P=None, R=None, T=None, rho=None, **kwargs):
        return
