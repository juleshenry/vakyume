from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_1_11__M_cap import eqn_1_11__M
from .eqn_1_11__P_cap import eqn_1_11__P
from .eqn_1_11__T_cap import eqn_1_11__T
from .eqn_1_11__W_cap import eqn_1_11__W
from .eqn_1_11__q import eqn_1_11__q


class VacuumTheory:
    eqn_1_11__M = eqn_1_11__M
    eqn_1_11__P = eqn_1_11__P
    eqn_1_11__T = eqn_1_11__T
    eqn_1_11__W = eqn_1_11__W
    eqn_1_11__q = eqn_1_11__q

    @kwasak
    def eqn_1_11(self, M=None, P=None, T=None, W=None, q=None):
        return
