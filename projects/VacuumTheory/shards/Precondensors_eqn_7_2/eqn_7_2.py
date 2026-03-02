from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_2__P_i_0_cap import eqn_7_2__P_i_0
from .eqn_7_2__p_i import eqn_7_2__p_i
from .eqn_7_2__x_i import eqn_7_2__x_i

class Precondensors:
    eqn_7_2__P_i_0 = eqn_7_2__P_i_0
    eqn_7_2__p_i = eqn_7_2__p_i
    eqn_7_2__x_i = eqn_7_2__x_i

    @kwasak_static
    def eqn_7_2(self, P_i_0=None, p_i=None, x_i=None):
        return
