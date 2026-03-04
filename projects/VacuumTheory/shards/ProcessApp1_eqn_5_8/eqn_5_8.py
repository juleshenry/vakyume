from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_5_8__P_cap_0_1 import eqn_5_8__P_0_1
from .eqn_5_8__P_cap_0_2 import eqn_5_8__P_0_2
from .eqn_5_8__alpha_12 import eqn_5_8__alpha_12
from .eqn_5_8__gamma_1 import eqn_5_8__gamma_1
from .eqn_5_8__gamma_2 import eqn_5_8__gamma_2

class ProcessApp1:
    eqn_5_8__P_0_1 = eqn_5_8__P_0_1
    eqn_5_8__P_0_2 = eqn_5_8__P_0_2
    eqn_5_8__alpha_12 = eqn_5_8__alpha_12
    eqn_5_8__gamma_1 = eqn_5_8__gamma_1
    eqn_5_8__gamma_2 = eqn_5_8__gamma_2

    @kwasak
    def eqn_5_8(self, P_0_1=None, P_0_2=None, alpha_12=None, gamma_1=None, gamma_2=None):
        return
