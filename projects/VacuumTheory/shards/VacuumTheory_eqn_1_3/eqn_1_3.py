from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_3__T_cap import eqn_1_3__T
from .eqn_1_3__k import eqn_1_3__k
from .eqn_1_3__m import eqn_1_3__m
from .eqn_1_3__v import eqn_1_3__v

class VacuumTheory:
    eqn_1_3__T = staticmethod(eqn_1_3__T)
    eqn_1_3__k = staticmethod(eqn_1_3__k)
    eqn_1_3__m = staticmethod(eqn_1_3__m)
    eqn_1_3__v = staticmethod(eqn_1_3__v)

    @kwasak_static
    def eqn_1_3(T=None, k=None, m=None, v=None, **kwargs):
        return
