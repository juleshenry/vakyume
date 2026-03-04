from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_1_9__M_cap import eqn_1_9__M
from .eqn_1_9__P_cap import eqn_1_9__P
from .eqn_1_9__R_cap import eqn_1_9__R
from .eqn_1_9__T_cap import eqn_1_9__T
from .eqn_1_9__rho import eqn_1_9__rho

class VacuumTheory:
    eqn_1_9__M = eqn_1_9__M
    eqn_1_9__P = eqn_1_9__P
    eqn_1_9__R = eqn_1_9__R
    eqn_1_9__T = eqn_1_9__T
    eqn_1_9__rho = eqn_1_9__rho

    @kwasak
    def eqn_1_9(self, M=None, P=None, R=None, T=None, rho=None):
        return
