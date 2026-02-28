from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_9__A_C_cap import eqn_3_9__A_C
from .eqn_3_9__H_1_cap import eqn_3_9__H_1
from .eqn_3_9__H_2_cap import eqn_3_9__H_2
from .eqn_3_9__P_cap import eqn_3_9__P
from .eqn_3_9__V_cap import eqn_3_9__V

class PressMgmt:
    eqn_3_9__A_C = staticmethod(eqn_3_9__A_C)
    eqn_3_9__H_1 = staticmethod(eqn_3_9__H_1)
    eqn_3_9__H_2 = staticmethod(eqn_3_9__H_2)
    eqn_3_9__P = staticmethod(eqn_3_9__P)
    eqn_3_9__V = staticmethod(eqn_3_9__V)

    @kwasak_static
    def eqn_3_9(A_C=None, H_1=None, H_2=None, P=None, V=None, **kwargs):
        return
