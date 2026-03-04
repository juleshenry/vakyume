from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_1_8__M import eqn_1_8__M
from .eqn_1_8__P import eqn_1_8__P
from .eqn_1_8__R import eqn_1_8__R
from .eqn_1_8__T import eqn_1_8__T
from .eqn_1_8__V import eqn_1_8__V
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
