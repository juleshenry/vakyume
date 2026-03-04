from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_15__U import eqn_7_15__U
from .eqn_7_15__sum_R import eqn_7_15__sum_R

class Precondensors:
    eqn_7_15__U = eqn_7_15__U
    eqn_7_15__sum_R = eqn_7_15__sum_R

    @kwasak
    def eqn_7_15(self, U=None, sum_R=None):
        return
