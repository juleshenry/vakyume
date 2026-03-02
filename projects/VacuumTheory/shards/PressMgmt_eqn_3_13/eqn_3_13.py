from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_3_13__H_1_cap import eqn_3_13__H_1
from .eqn_3_13__H_2_cap import eqn_3_13__H_2
from .eqn_3_13__KAPPA_2_cap import eqn_3_13__KAPPA_2
from .eqn_3_13__P_cap import eqn_3_13__P

class PressMgmt:
    eqn_3_13__H_1 = eqn_3_13__H_1
    eqn_3_13__H_2 = eqn_3_13__H_2
    eqn_3_13__KAPPA_2 = eqn_3_13__KAPPA_2
    eqn_3_13__P = eqn_3_13__P

    @kwasak_static
    def eqn_3_13(self, H_1=None, H_2=None, KAPPA_2=None, P=None):
        return
