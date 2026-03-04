from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_24_6__c import eqn_24_6__c
from .eqn_24_6__u import eqn_24_6__u
from .eqn_24_6__ux import eqn_24_6__ux
from .eqn_24_6__v import eqn_24_6__v


class TheTheoryOfSpecialRelativity:
    eqn_24_6__c = eqn_24_6__c
    eqn_24_6__u = eqn_24_6__u
    eqn_24_6__ux = eqn_24_6__ux
    eqn_24_6__v = eqn_24_6__v

    @kwasak
    def eqn_24_6(self, c=None, u=None, ux=None, v=None):
        """
        ux := relative velocity
        u := initial velocity
        v := velocity of the train
        c := speed of light
        """
        return
