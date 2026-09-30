from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_3__H_cap_1 import eqn_3_3__H_1
from .eqn_3_3__H_cap_2 import eqn_3_3__H_2
from .eqn_3_3__P_cap import eqn_3_3__P
from .eqn_3_3__P_cap_P_cap import eqn_3_3__P_P


class PressMgmt:
    eqn_3_3__H_1 = eqn_3_3__H_1
    eqn_3_3__H_2 = eqn_3_3__H_2
    eqn_3_3__P = eqn_3_3__P
    eqn_3_3__P_P = eqn_3_3__P_P

    @kwasak
    def eqn_3_3(self, H_1=None, H_2=None, P=None, P_P=None):
        return
