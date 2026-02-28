from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_9_1__A_cap import eqn_9_1__A
from .eqn_9_1__rho_s import eqn_9_1__rho_s
from .eqn_9_1__v import eqn_9_1__v
from .eqn_9_1__w_s import eqn_9_1__w_s

class SteamJetInjectors:
    eqn_9_1__A = staticmethod(eqn_9_1__A)
    eqn_9_1__rho_s = staticmethod(eqn_9_1__rho_s)
    eqn_9_1__v = staticmethod(eqn_9_1__v)
    eqn_9_1__w_s = staticmethod(eqn_9_1__w_s)

    @kwasak_static
    def eqn_9_1(A=None, rho_s=None, v=None, w_s=None, **kwargs):
        return
