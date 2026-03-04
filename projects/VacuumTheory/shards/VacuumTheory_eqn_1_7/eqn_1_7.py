from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_1_7__R import eqn_1_7__R
from .eqn_1_7__T import eqn_1_7__T
from .eqn_1_7__V import eqn_1_7__V
from .eqn_1_7__n import eqn_1_7__n
from .eqn_1_7__p import eqn_1_7__p


class VacuumTheory:
    eqn_1_7__R = eqn_1_7__R
    eqn_1_7__T = eqn_1_7__T
    eqn_1_7__V = eqn_1_7__V
    eqn_1_7__n = eqn_1_7__n
    eqn_1_7__p = eqn_1_7__p

    @kwasak
    def eqn_1_7(self, R=None, T=None, V=None, n=None, p=None):
        return
