from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_11__A_C_cap import eqn_3_11__A_C
from .eqn_3_11__H_2_cap import eqn_3_11__H_2
from .eqn_3_11__P_cap import eqn_3_11__P
from .eqn_3_11__V_cap import eqn_3_11__V

class PressMgmt:
    eqn_3_11__A_C = eqn_3_11__A_C
    eqn_3_11__H_2 = eqn_3_11__H_2
    eqn_3_11__P = eqn_3_11__P
    eqn_3_11__V = eqn_3_11__V

    @kwasak_static
    def eqn_3_11(self, A_C=None, H_2=None, P=None, V=None):
        return
