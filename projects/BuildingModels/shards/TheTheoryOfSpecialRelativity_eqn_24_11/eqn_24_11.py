from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_24_11__E_cap import eqn_24_11__E
from .eqn_24_11__c import eqn_24_11__c
from .eqn_24_11__m_0 import eqn_24_11__m_0
from .eqn_24_11__p import eqn_24_11__p


class TheTheoryOfSpecialRelativity:
    eqn_24_11__E = eqn_24_11__E
    eqn_24_11__c = eqn_24_11__c
    eqn_24_11__m_0 = eqn_24_11__m_0
    eqn_24_11__p = eqn_24_11__p

    @kwasak
    def eqn_24_11(self, E=None, c=None, m_0=None, p=None):
        """
        E := energy
        p := momentum
        c := speed of light
        """
        return
