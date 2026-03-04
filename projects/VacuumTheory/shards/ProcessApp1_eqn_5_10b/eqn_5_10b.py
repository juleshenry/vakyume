from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_5_10b__L_cap_0 import eqn_5_10b__L_0
from .eqn_5_10b__R_cap import eqn_5_10b__R
from .eqn_5_10b__V_cap_1 import eqn_5_10b__V_1

class ProcessApp1:
    eqn_5_10b__L_0 = eqn_5_10b__L_0
    eqn_5_10b__R = eqn_5_10b__R
    eqn_5_10b__V_1 = eqn_5_10b__V_1

    @kwasak
    def eqn_5_10b(self, L_0=None, R=None, V_1=None):
        return
