from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_5_4__P_cap import eqn_5_4__P
from .eqn_5_4__P_cap_0_i import eqn_5_4__P_0_i
from .eqn_5_4__x_i import eqn_5_4__x_i
from .eqn_5_4__y_i import eqn_5_4__y_i


class ProcessApp1:
    eqn_5_4__P = eqn_5_4__P
    eqn_5_4__P_0_i = eqn_5_4__P_0_i
    eqn_5_4__x_i = eqn_5_4__x_i
    eqn_5_4__y_i = eqn_5_4__y_i

    @kwasak
    def eqn_5_4(self, P=None, P_0_i=None, x_i=None, y_i=None):
        return
