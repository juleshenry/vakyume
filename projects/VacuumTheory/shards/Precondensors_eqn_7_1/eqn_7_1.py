from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_1__P import eqn_7_1__P
from .eqn_7_1__p_i import eqn_7_1__p_i
from .eqn_7_1__y_i import eqn_7_1__y_i

class Precondensors:
    eqn_7_1__P = eqn_7_1__P
    eqn_7_1__p_i = eqn_7_1__p_i
    eqn_7_1__y_i = eqn_7_1__y_i

    @kwasak
    def eqn_7_1(self, P=None, p_i=None, y_i=None):
        return
