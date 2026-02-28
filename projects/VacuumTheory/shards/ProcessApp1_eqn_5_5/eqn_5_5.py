from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_5__P_0_1_cap import eqn_5_5__P_0_1
from .eqn_5_5__P_0_2_cap import eqn_5_5__P_0_2
from .eqn_5_5__alpha_12 import eqn_5_5__alpha_12

class ProcessApp1:
    eqn_5_5__P_0_1 = staticmethod(eqn_5_5__P_0_1)
    eqn_5_5__P_0_2 = staticmethod(eqn_5_5__P_0_2)
    eqn_5_5__alpha_12 = staticmethod(eqn_5_5__alpha_12)

    @kwasak_static
    def eqn_5_5(P_0_1=None, P_0_2=None, alpha_12=None, **kwargs):
        return
