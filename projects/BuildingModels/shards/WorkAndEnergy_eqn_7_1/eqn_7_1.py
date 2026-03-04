from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_7_1__F_cap import eqn_7_1__F
from .eqn_7_1__W_cap import eqn_7_1__W
from .eqn_7_1__d import eqn_7_1__d


class WorkAndEnergy:
    eqn_7_1__F = eqn_7_1__F
    eqn_7_1__W = eqn_7_1__W
    eqn_7_1__d = eqn_7_1__d

    @kwasak
    def eqn_7_1(self, F=None, W=None, d=None):
        """
        W := work
        F := force
        d := displacement
        x_0 := initial position
        x_1 := final position
        """
        return
