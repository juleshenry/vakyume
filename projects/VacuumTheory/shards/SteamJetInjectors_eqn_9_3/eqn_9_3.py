from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_9_3__P_s_cap import eqn_9_3__P_s
from .eqn_9_3__V_cap import eqn_9_3__V
from .eqn_9_3__t_e import eqn_9_3__t_e
from .eqn_9_3__w_j import eqn_9_3__w_j

class SteamJetInjectors:
    eqn_9_3__P_s = staticmethod(eqn_9_3__P_s)
    eqn_9_3__V = staticmethod(eqn_9_3__V)
    eqn_9_3__t_e = staticmethod(eqn_9_3__t_e)
    eqn_9_3__w_j = staticmethod(eqn_9_3__w_j)

    @kwasak_static
    def eqn_9_3(P_s=None, V=None, t_e=None, w_j=None, **kwargs):
        return
