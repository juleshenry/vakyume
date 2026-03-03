from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_3_12__H_2_cap import eqn_3_12__H_2
from .eqn_3_12__KAPPA_1_cap import eqn_3_12__KAPPA_1
from .eqn_3_12__P_cap import eqn_3_12__P

class PressMgmt:
    eqn_3_12__H_2 = eqn_3_12__H_2
    eqn_3_12__KAPPA_1 = eqn_3_12__KAPPA_1
    eqn_3_12__P = eqn_3_12__P

    @kwasak
    def eqn_3_12(self, H_2=None, KAPPA_1=None, P=None):
        return
