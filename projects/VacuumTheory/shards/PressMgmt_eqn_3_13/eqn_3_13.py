from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_13__H_cap_1 import eqn_3_13__H_1
from .eqn_3_13__H_cap_2 import eqn_3_13__H_2
from .eqn_3_13__K_capA_capP_capP_capA_cap_2 import eqn_3_13__KAPPA_2
from .eqn_3_13__P_cap import eqn_3_13__P

class PressMgmt:
    eqn_3_13__H_1 = eqn_3_13__H_1
    eqn_3_13__H_2 = eqn_3_13__H_2
    eqn_3_13__KAPPA_2 = eqn_3_13__KAPPA_2
    eqn_3_13__P = eqn_3_13__P

    @kwasak
    def eqn_3_13(self, H_1=None, H_2=None, KAPPA_2=None, P=None):
        return
