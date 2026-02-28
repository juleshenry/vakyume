from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_1__P_cap import eqn_7_1__P
from .eqn_7_1__p_i import eqn_7_1__p_i
from .eqn_7_1__y_i import eqn_7_1__y_i

class Precondensors:
    eqn_7_1__P = staticmethod(eqn_7_1__P)
    eqn_7_1__p_i = staticmethod(eqn_7_1__p_i)
    eqn_7_1__y_i = staticmethod(eqn_7_1__y_i)

    @kwasak_static
    def eqn_7_1(P=None, p_i=None, y_i=None, **kwargs):
        return
