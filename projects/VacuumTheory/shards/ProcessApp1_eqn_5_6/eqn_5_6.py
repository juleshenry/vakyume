from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_6__P_0_i_cap import eqn_5_6__P_0_i
from .eqn_5_6__gamma_i import eqn_5_6__gamma_i
from .eqn_5_6__p_i import eqn_5_6__p_i
from .eqn_5_6__x_i import eqn_5_6__x_i

class ProcessApp1:
    eqn_5_6__P_0_i = eqn_5_6__P_0_i
    eqn_5_6__gamma_i = eqn_5_6__gamma_i
    eqn_5_6__p_i = eqn_5_6__p_i
    eqn_5_6__x_i = eqn_5_6__x_i

    @kwasak_static
    def eqn_5_6(self, P_0_i=None, gamma_i=None, p_i=None, x_i=None):
        return
