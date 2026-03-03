from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_14__M_cap import eqn_2_14__M
from .eqn_2_14__R_cap import eqn_2_14__R
from .eqn_2_14__T_cap import eqn_2_14__T
from .eqn_2_14__g_c import eqn_2_14__g_c
from .eqn_2_14__k import eqn_2_14__k
from .eqn_2_14__v_s import eqn_2_14__v_s

class FluidFlowVacuumLines:
    eqn_2_14__M = eqn_2_14__M
    eqn_2_14__R = eqn_2_14__R
    eqn_2_14__T = eqn_2_14__T
    eqn_2_14__g_c = eqn_2_14__g_c
    eqn_2_14__k = eqn_2_14__k
    eqn_2_14__v_s = eqn_2_14__v_s

    @kwasak
    def eqn_2_14(self, M=None, R=None, T=None, g_c=None, k=None, v_s=None):
        return
