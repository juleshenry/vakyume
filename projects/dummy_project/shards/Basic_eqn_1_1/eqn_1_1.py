from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_1__x import eqn_1_1__x
from .eqn_1_1__y import eqn_1_1__y
from .eqn_1_1__z import eqn_1_1__z

class Basic:
    eqn_1_1__x = staticmethod(eqn_1_1__x)
    eqn_1_1__y = staticmethod(eqn_1_1__y)
    eqn_1_1__z = staticmethod(eqn_1_1__z)

    @kwasak_static
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return
