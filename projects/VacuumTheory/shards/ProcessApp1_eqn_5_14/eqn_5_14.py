from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_14__M_cap import eqn_5_14__M
from .eqn_5_14__P_0_cap import eqn_5_14__P_0
from .eqn_5_14__T_cap import eqn_5_14__T
from .eqn_5_14__W_E_cap import eqn_5_14__W_E

class ProcessApp1:
    eqn_5_14__M = eqn_5_14__M
    eqn_5_14__P_0 = eqn_5_14__P_0
    eqn_5_14__T = eqn_5_14__T
    eqn_5_14__W_E = eqn_5_14__W_E

    @kwasak_static
    def eqn_5_14(self, M=None, P_0=None, T=None, W_E=None):
        return
