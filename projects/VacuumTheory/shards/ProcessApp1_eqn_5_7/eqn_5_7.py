from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_7__P_cap import eqn_5_7__P
from .eqn_5_7__P_0_i_cap import eqn_5_7__P_0_i
from .eqn_5_7__gamma_i import eqn_5_7__gamma_i
from .eqn_5_7__x_i import eqn_5_7__x_i
from .eqn_5_7__y_i import eqn_5_7__y_i

class ProcessApp1:
    eqn_5_7__P = eqn_5_7__P
    eqn_5_7__P_0_i = eqn_5_7__P_0_i
    eqn_5_7__gamma_i = eqn_5_7__gamma_i
    eqn_5_7__x_i = eqn_5_7__x_i
    eqn_5_7__y_i = eqn_5_7__y_i

    @kwasak_static
    def eqn_5_7(self, P=None, P_0_i=None, gamma_i=None, x_i=None, y_i=None):
        return
