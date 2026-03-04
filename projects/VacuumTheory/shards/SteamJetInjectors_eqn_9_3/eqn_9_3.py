from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_9_3__P_s import eqn_9_3__P_s
from .eqn_9_3__V import eqn_9_3__V
from .eqn_9_3__t_e import eqn_9_3__t_e
from .eqn_9_3__w_j import eqn_9_3__w_j

class SteamJetInjectors:
    eqn_9_3__P_s = eqn_9_3__P_s
    eqn_9_3__V = eqn_9_3__V
    eqn_9_3__t_e = eqn_9_3__t_e
    eqn_9_3__w_j = eqn_9_3__w_j

    @kwasak
    def eqn_9_3(self, P_s=None, V=None, t_e=None, w_j=None):
        """
        t_e := time required to evacuate system, minutes
        P_s := design suction pressure of the ejector, torr
        V := free volume of process system, ft^3
        w_j := ejector capacity, 70 deg_F basis, lb/hr
        """
        return
