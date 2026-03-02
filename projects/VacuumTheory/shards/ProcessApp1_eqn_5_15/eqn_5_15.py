from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_15__M_1_cap import eqn_5_15__M_1
from .eqn_5_15__M_2_cap import eqn_5_15__M_2
from .eqn_5_15__P_0_1_cap import eqn_5_15__P_0_1
from .eqn_5_15__P_0_2_cap import eqn_5_15__P_0_2
from .eqn_5_15__a_M_12_cap import eqn_5_15__a_M_12

class ProcessApp1:
    eqn_5_15__M_1 = eqn_5_15__M_1
    eqn_5_15__M_2 = eqn_5_15__M_2
    eqn_5_15__P_0_1 = eqn_5_15__P_0_1
    eqn_5_15__P_0_2 = eqn_5_15__P_0_2
    eqn_5_15__a_M_12 = eqn_5_15__a_M_12

    @kwasak_static
    def eqn_5_15(self, M_1=None, M_2=None, P_0_1=None, P_0_2=None, a_M_12=None):
        return
