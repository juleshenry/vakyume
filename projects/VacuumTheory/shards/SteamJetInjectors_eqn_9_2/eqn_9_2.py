from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_9_2__P_m_cap import eqn_9_2__P_m
from .eqn_9_2__d_n import eqn_9_2__d_n
from .eqn_9_2__rho_s import eqn_9_2__rho_s
from .eqn_9_2__w_s import eqn_9_2__w_s

class SteamJetInjectors:
    eqn_9_2__P_m = eqn_9_2__P_m
    eqn_9_2__d_n = eqn_9_2__d_n
    eqn_9_2__rho_s = eqn_9_2__rho_s
    eqn_9_2__w_s = eqn_9_2__w_s

    @kwasak
    def eqn_9_2(self, P_m=None, d_n=None, rho_s=None, w_s=None):
        return
