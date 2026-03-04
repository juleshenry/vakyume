from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_16_3__F_cap_g import eqn_16_3__F_g
from .eqn_16_3__G_cap import eqn_16_3__G
from .eqn_16_3__m_e import eqn_16_3__m_e
from .eqn_16_3__m_p import eqn_16_3__m_p
from .eqn_16_3__r import eqn_16_3__r


class ElectricChargesAndFields:
    eqn_16_3__F_g = eqn_16_3__F_g
    eqn_16_3__G = eqn_16_3__G
    eqn_16_3__m_e = eqn_16_3__m_e
    eqn_16_3__m_p = eqn_16_3__m_p
    eqn_16_3__r = eqn_16_3__r

    @kwasak
    def eqn_16_3(self, F_g=None, G=None, m_e=None, m_p=None, r=None):
        """
        G := gravitational constant
        m_e := mass of electron
        m_p := mass of proton
        r := distance between electron and proton
        """
        return
