from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_13a__n import eqn_1_13a__n
from .eqn_1_13a__n_a import eqn_1_13a__n_a
from .eqn_1_13a__y_a import eqn_1_13a__y_a

class VacuumTheory:
    eqn_1_13a__n = staticmethod(eqn_1_13a__n)
    eqn_1_13a__n_a = staticmethod(eqn_1_13a__n_a)
    eqn_1_13a__y_a = staticmethod(eqn_1_13a__y_a)

    @kwasak_static
    def eqn_1_13a(n=None, n_a=None, y_a=None, **kwargs):
        return
