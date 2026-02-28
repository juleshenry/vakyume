from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_9_5__V_cap import eqn_9_5__V
from .eqn_9_5__r_h import eqn_9_5__r_h
from .eqn_9_5__t_h import eqn_9_5__t_h
from .eqn_9_5__w_h import eqn_9_5__w_h

class SteamJetInjectors:
    eqn_9_5__V = staticmethod(eqn_9_5__V)
    eqn_9_5__r_h = staticmethod(eqn_9_5__r_h)
    eqn_9_5__t_h = staticmethod(eqn_9_5__t_h)
    eqn_9_5__w_h = staticmethod(eqn_9_5__w_h)

    @kwasak_static
    def eqn_9_5(V=None, r_h=None, t_h=None, w_h=None, **kwargs):
        return
