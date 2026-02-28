from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_15__U_cap import eqn_7_15__U
from .eqn_7_15__sum_R_cap import eqn_7_15__sum_R

class Precondensors:
    eqn_7_15__U = staticmethod(eqn_7_15__U)
    eqn_7_15__sum_R = staticmethod(eqn_7_15__sum_R)

    @kwasak_static
    def eqn_7_15(U=None, sum_R=None, **kwargs):
        return
