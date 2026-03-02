from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_11__D_cap import eqn_2_11__D
from .eqn_2_11__L_cap import eqn_2_11__L
from .eqn_2_11__f import eqn_2_11__f
from .eqn_2_11__g_c import eqn_2_11__g_c
from .eqn_2_11__h_r import eqn_2_11__h_r
from .eqn_2_11__v import eqn_2_11__v

class FluidFlowVacuumLines:
    eqn_2_11__D = eqn_2_11__D
    eqn_2_11__L = eqn_2_11__L
    eqn_2_11__f = eqn_2_11__f
    eqn_2_11__g_c = eqn_2_11__g_c
    eqn_2_11__h_r = eqn_2_11__h_r
    eqn_2_11__v = eqn_2_11__v

    @kwasak_static
    def eqn_2_11(self, D=None, L=None, f=None, g_c=None, h_r=None, v=None):
        return
