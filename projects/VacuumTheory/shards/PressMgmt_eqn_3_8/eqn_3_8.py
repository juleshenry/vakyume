from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_3_8__A_C_cap import eqn_3_8__A_C
from .eqn_3_8__H_2_cap import eqn_3_8__H_2
from .eqn_3_8__V_P_cap import eqn_3_8__V_P

class PressMgmt:
    eqn_3_8__A_C = eqn_3_8__A_C
    eqn_3_8__H_2 = eqn_3_8__H_2
    eqn_3_8__V_P = eqn_3_8__V_P

    @kwasak
    def eqn_3_8(self, A_C=None, H_2=None, V_P=None):
        return
