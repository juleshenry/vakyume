from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_5_10c__D_cap import eqn_5_10c__D
from .eqn_5_10c__L_0_cap import eqn_5_10c__L_0
from .eqn_5_10c__R_cap import eqn_5_10c__R

class ProcessApp1:
    eqn_5_10c__D = eqn_5_10c__D
    eqn_5_10c__L_0 = eqn_5_10c__L_0
    eqn_5_10c__R = eqn_5_10c__R

    @kwasak
    def eqn_5_10c(self, D=None, L_0=None, R=None):
        return
