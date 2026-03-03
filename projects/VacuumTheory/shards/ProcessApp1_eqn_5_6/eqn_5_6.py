from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_6__P_cap_0_i import eqn_5_6__P_0_i
from .eqn_5_6__gamma_i import eqn_5_6__gamma_i
from .eqn_5_6__p_i import eqn_5_6__p_i
from .eqn_5_6__x_i import eqn_5_6__x_i


class ProcessApp1:
    eqn_5_6__P_0_i = eqn_5_6__P_0_i
    eqn_5_6__gamma_i = eqn_5_6__gamma_i
    eqn_5_6__p_i = eqn_5_6__p_i
    eqn_5_6__x_i = eqn_5_6__x_i

    @kwasak
    def eqn_5_6(self, P_0_i=None, gamma_i=None, p_i=None, x_i=None):
        return
