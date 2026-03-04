from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_10_9__T_c import eqn_10_9__T_c
from .eqn_10_9__T_s import eqn_10_9__T_s
from .eqn_10_9__delta_T import eqn_10_9__delta_T


class LiquidRing:
    eqn_10_9__T_c = eqn_10_9__T_c
    eqn_10_9__T_s = eqn_10_9__T_s
    eqn_10_9__delta_T = eqn_10_9__delta_T

    @kwasak
    def eqn_10_9(self, T_c=None, T_s=None, delta_T=None):
        return
