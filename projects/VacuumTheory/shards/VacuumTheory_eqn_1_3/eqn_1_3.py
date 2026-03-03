from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_1_3__T_cap import eqn_1_3__T
from .eqn_1_3__k import eqn_1_3__k
from .eqn_1_3__m import eqn_1_3__m
from .eqn_1_3__v import eqn_1_3__v

class VacuumTheory:
    eqn_1_3__T = eqn_1_3__T
    eqn_1_3__k = eqn_1_3__k
    eqn_1_3__m = eqn_1_3__m
    eqn_1_3__v = eqn_1_3__v

    @kwasak
    def eqn_1_3(self, T=None, k=None, m=None, v=None):
        """
        k:= boltzmann constant
        kboltz:= 1.38e-16
        avogad:= 6.02e23
        """
        return
