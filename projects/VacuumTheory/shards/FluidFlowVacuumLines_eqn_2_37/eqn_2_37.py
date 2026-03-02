from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_37__A_cap import eqn_2_37__A
from .eqn_2_37__C_cap import eqn_2_37__C
from .eqn_2_37__F_t_cap import eqn_2_37__F_t
from .eqn_2_37__M_cap import eqn_2_37__M
from .eqn_2_37__T_cap import eqn_2_37__T

class FluidFlowVacuumLines:
    eqn_2_37__A = eqn_2_37__A
    eqn_2_37__C = eqn_2_37__C
    eqn_2_37__F_t = eqn_2_37__F_t
    eqn_2_37__M = eqn_2_37__M
    eqn_2_37__T = eqn_2_37__T

    @kwasak_static
    def eqn_2_37(self, A=None, C=None, F_t=None, M=None, T=None):
        return
