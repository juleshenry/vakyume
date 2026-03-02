from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_1__K_i_cap import eqn_5_1__K_i
from .eqn_5_1__x_i import eqn_5_1__x_i
from .eqn_5_1__y_i import eqn_5_1__y_i

class ProcessApp1:
    eqn_5_1__K_i = eqn_5_1__K_i
    eqn_5_1__x_i = eqn_5_1__x_i
    eqn_5_1__y_i = eqn_5_1__y_i

    @kwasak_static
    def eqn_5_1(self, K_i=None, x_i=None, y_i=None):
        return
