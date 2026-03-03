from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_1_13b__P_cap import eqn_1_13b__P
from .eqn_1_13b__p_a import eqn_1_13b__p_a
from .eqn_1_13b__y_a import eqn_1_13b__y_a

class VacuumTheory:
    eqn_1_13b__P = eqn_1_13b__P
    eqn_1_13b__p_a = eqn_1_13b__p_a
    eqn_1_13b__y_a = eqn_1_13b__y_a

    @kwasak
    def eqn_1_13b(self, P=None, p_a=None, y_a=None):
        return
