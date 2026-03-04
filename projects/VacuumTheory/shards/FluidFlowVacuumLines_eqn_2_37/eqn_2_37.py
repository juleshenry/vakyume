from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_37__A import eqn_2_37__A
from .eqn_2_37__C import eqn_2_37__C
from .eqn_2_37__F_t import eqn_2_37__F_t
from .eqn_2_37__M import eqn_2_37__M
from .eqn_2_37__T import eqn_2_37__T


class FluidFlowVacuumLines:
    eqn_2_37__A = eqn_2_37__A
    eqn_2_37__C = eqn_2_37__C
    eqn_2_37__F_t = eqn_2_37__F_t
    eqn_2_37__M = eqn_2_37__M
    eqn_2_37__T = eqn_2_37__T

    @kwasak
    def eqn_2_37(self, A=None, C=None, F_t=None, M=None, T=None):
        """
        F_t:= 1, for an aperture
        """
        return
