from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_14__M_cap import eqn_2_14__M
from .eqn_2_14__R_cap import eqn_2_14__R
from .eqn_2_14__T_cap import eqn_2_14__T
from .eqn_2_14__g_c import eqn_2_14__g_c
from .eqn_2_14__k import eqn_2_14__k
from .eqn_2_14__v_s import eqn_2_14__v_s

class FluidFlowVacuumLines:
    eqn_2_14__M = staticmethod(eqn_2_14__M)
    eqn_2_14__R = staticmethod(eqn_2_14__R)
    eqn_2_14__T = staticmethod(eqn_2_14__T)
    eqn_2_14__g_c = staticmethod(eqn_2_14__g_c)
    eqn_2_14__k = staticmethod(eqn_2_14__k)
    eqn_2_14__v_s = staticmethod(eqn_2_14__v_s)

    @kwasak_static
    def eqn_2_14(M=None, R=None, T=None, g_c=None, k=None, v_s=None, **kwargs):
        return
