from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_9_5__V_cap import eqn_9_5__V
from .eqn_9_5__r_h import eqn_9_5__r_h
from .eqn_9_5__t_h import eqn_9_5__t_h
from .eqn_9_5__w_h import eqn_9_5__w_h

class SteamJetInjectors:
    eqn_9_5__V = eqn_9_5__V
    eqn_9_5__r_h = eqn_9_5__r_h
    eqn_9_5__t_h = eqn_9_5__t_h
    eqn_9_5__w_h = eqn_9_5__w_h

    @kwasak
    def eqn_9_5(self, V=None, r_h=None, t_h=None, w_h=None):
        return
