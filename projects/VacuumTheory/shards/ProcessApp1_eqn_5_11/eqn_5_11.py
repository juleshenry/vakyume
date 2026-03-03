from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_11__B_cap import eqn_5_11__B
from .eqn_5_11__L_cap_N_cap import eqn_5_11__L_N
from .eqn_5_11__V_cap_0 import eqn_5_11__V_0

class ProcessApp1:
    eqn_5_11__B = eqn_5_11__B
    eqn_5_11__L_N = eqn_5_11__L_N
    eqn_5_11__V_0 = eqn_5_11__V_0

    @kwasak
    def eqn_5_11(self, B=None, L_N=None, V_0=None):
        return
