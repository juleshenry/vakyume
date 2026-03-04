from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_1_10__P_1 import eqn_1_10__P_1
from .eqn_1_10__P_2 import eqn_1_10__P_2
from .eqn_1_10__T_1 import eqn_1_10__T_1
from .eqn_1_10__T_2 import eqn_1_10__T_2
from .eqn_1_10__V_1 import eqn_1_10__V_1
from .eqn_1_10__V_2 import eqn_1_10__V_2

class VacuumTheory:
    eqn_1_10__P_1 = eqn_1_10__P_1
    eqn_1_10__P_2 = eqn_1_10__P_2
    eqn_1_10__T_1 = eqn_1_10__T_1
    eqn_1_10__T_2 = eqn_1_10__T_2
    eqn_1_10__V_1 = eqn_1_10__V_1
    eqn_1_10__V_2 = eqn_1_10__V_2

    @kwasak
    def eqn_1_10(self, P_1=None, P_2=None, T_1=None, T_2=None, V_1=None, V_2=None):
        return
