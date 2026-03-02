from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_9_4__AEL_cap import eqn_9_4__AEL
from .eqn_9_4__SC_cap import eqn_9_4__SC
from .eqn_9_4__r import eqn_9_4__r
from .eqn_9_4__w_s import eqn_9_4__w_s

class SteamJetInjectors:
    eqn_9_4__AEL = eqn_9_4__AEL
    eqn_9_4__SC = eqn_9_4__SC
    eqn_9_4__r = eqn_9_4__r
    eqn_9_4__w_s = eqn_9_4__w_s

    @kwasak_static
    def eqn_9_4(self, AEL=None, SC=None, r=None, w_s=None):
        return
