from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_10b__L_0_cap import eqn_5_10b__L_0
from .eqn_5_10b__R_cap import eqn_5_10b__R
from .eqn_5_10b__V_1_cap import eqn_5_10b__V_1

class ProcessApp1:
    eqn_5_10b__L_0 = eqn_5_10b__L_0
    eqn_5_10b__R = eqn_5_10b__R
    eqn_5_10b__V_1 = eqn_5_10b__V_1

    @kwasak_static
    def eqn_5_10b(self, L_0=None, R=None, V_1=None):
        return
