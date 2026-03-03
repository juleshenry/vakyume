from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_6__M_cap import eqn_7_6__M
from .eqn_7_6__P_cap import eqn_7_6__P
from .eqn_7_6__P_i_0_cap import eqn_7_6__P_i_0
from .eqn_7_6__W_air_cap import eqn_7_6__W_air
from .eqn_7_6__W_i_cap import eqn_7_6__W_i
from .eqn_7_6__p_c import eqn_7_6__p_c
from .eqn_7_6__x_i import eqn_7_6__x_i

class Precondensors:
    eqn_7_6__M = eqn_7_6__M
    eqn_7_6__P = eqn_7_6__P
    eqn_7_6__P_i_0 = eqn_7_6__P_i_0
    eqn_7_6__W_air = eqn_7_6__W_air
    eqn_7_6__W_i = eqn_7_6__W_i
    eqn_7_6__p_c = eqn_7_6__p_c
    eqn_7_6__x_i = eqn_7_6__x_i

    @kwasak_static
    def eqn_7_6(self, M=None, P=None, P_i_0=None, W_air=None, W_i=None, p_c=None, x_i=None):
        return
