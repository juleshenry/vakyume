from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_24_9__K_cap import eqn_24_9__K
from .eqn_24_9__c import eqn_24_9__c
from .eqn_24_9__m_0 import eqn_24_9__m_0
from .eqn_24_9__u import eqn_24_9__u


class TheTheoryOfSpecialRelativity:
    eqn_24_9__K = eqn_24_9__K
    eqn_24_9__c = eqn_24_9__c
    eqn_24_9__m_0 = eqn_24_9__m_0
    eqn_24_9__u = eqn_24_9__u

    @kwasak
    def eqn_24_9(self, K=None, c=None, m_0=None, u=None):
        """
        m_0 := rest mass
        c := speed of light
        """
        return
