from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_3_9__A_C_cap import eqn_3_9__A_C
from .eqn_3_9__H_1_cap import eqn_3_9__H_1
from .eqn_3_9__H_2_cap import eqn_3_9__H_2
from .eqn_3_9__P_cap import eqn_3_9__P
from .eqn_3_9__V_cap import eqn_3_9__V


class PressMgmt:
    eqn_3_9__A_C = eqn_3_9__A_C
    eqn_3_9__H_1 = eqn_3_9__H_1
    eqn_3_9__H_2 = eqn_3_9__H_2
    eqn_3_9__P = eqn_3_9__P
    eqn_3_9__V = eqn_3_9__V

    @kwasak
    def eqn_3_9(self, A_C=None, H_1=None, H_2=None, P=None, V=None):
        return
