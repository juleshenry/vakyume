from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_9__D_cap import eqn_5_9__D
from .eqn_5_9__L_cap_0 import eqn_5_9__L_0
from .eqn_5_9__V_cap_1 import eqn_5_9__V_1


class ProcessApp1:
    eqn_5_9__D = eqn_5_9__D
    eqn_5_9__L_0 = eqn_5_9__L_0
    eqn_5_9__V_1 = eqn_5_9__V_1

    @kwasak
    def eqn_5_9(self, D=None, L_0=None, V_1=None):
        return
