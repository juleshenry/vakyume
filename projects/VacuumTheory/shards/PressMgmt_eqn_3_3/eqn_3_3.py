from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_3__H_1_cap import eqn_3_3__H_1
from .eqn_3_3__H_2_cap import eqn_3_3__H_2
from .eqn_3_3__P_cap import eqn_3_3__P
from .eqn_3_3__P_P_cap import eqn_3_3__P_P

class PressMgmt:
    eqn_3_3__H_1 = eqn_3_3__H_1
    eqn_3_3__H_2 = eqn_3_3__H_2
    eqn_3_3__P = eqn_3_3__P
    eqn_3_3__P_P = eqn_3_3__P_P

    @kwasak_static
    def eqn_3_3(self, H_1=None, H_2=None, P=None, P_P=None):
        return
