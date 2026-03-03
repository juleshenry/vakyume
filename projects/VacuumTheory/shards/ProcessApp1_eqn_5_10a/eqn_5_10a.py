from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_10a__D_cap import eqn_5_10a__D
from .eqn_5_10a__L_0_cap import eqn_5_10a__L_0
from .eqn_5_10a__V_1_cap import eqn_5_10a__V_1

class ProcessApp1:
    eqn_5_10a__D = eqn_5_10a__D
    eqn_5_10a__L_0 = eqn_5_10a__L_0
    eqn_5_10a__V_1 = eqn_5_10a__V_1

    @kwasak
    def eqn_5_10a(self, D=None, L_0=None, V_1=None):
        return
