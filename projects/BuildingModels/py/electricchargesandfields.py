from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class ElectricChargesAndFields:
    @kwasak
    def eqn_16_3(self, F_g=None, G=None, m_e=None, m_p=None, r=None):
        """
        G := gravitational constant
        m_e := mass of electron
        m_p := mass of proton
        r := distance between electron and proton
        """
        return

    def eqn_16_3__F_g(self, G: float, m_e: float, m_p: float, r: float, **kwargs):
        # F_g = G * m_e * m_p / r ** 2
        result = []
        F_g = G * m_e * m_p / r**2
        result.append(F_g)
        return result

    def eqn_16_3__G(self, F_g: float, m_e: float, m_p: float, r: float, **kwargs):
        # F_g = G * m_e * m_p / r ** 2
        result = []
        G = F_g * r**2 / (m_e * m_p)
        result.append(G)
        return result

    def eqn_16_3__m_e(self, F_g: float, G: float, m_p: float, r: float, **kwargs):
        # F_g = G * m_e * m_p / r ** 2
        result = []
        m_e = F_g * r**2 / (G * m_p)
        result.append(m_e)
        return result

    def eqn_16_3__m_p(self, F_g: float, G: float, m_e: float, r: float, **kwargs):
        # F_g = G * m_e * m_p / r ** 2
        result = []
        m_p = F_g * r**2 / (G * m_e)
        result.append(m_p)
        return result

    def eqn_16_3__r(self, F_g: float, G: float, m_e: float, m_p: float, **kwargs):
        # F_g = G * m_e * m_p / r ** 2
        result = []
        r = -sqrt(G * m_e * m_p / F_g)
        result.append(r)
        r = sqrt(G * m_e * m_p / F_g)
        result.append(r)
        return result
