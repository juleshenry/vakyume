from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_9_1__A import eqn_9_1__A
from .eqn_9_1__rho_s import eqn_9_1__rho_s
from .eqn_9_1__v import eqn_9_1__v
from .eqn_9_1__w_s import eqn_9_1__w_s

class SteamJetInjectors:
    eqn_9_1__A = eqn_9_1__A
    eqn_9_1__rho_s = eqn_9_1__rho_s
    eqn_9_1__v = eqn_9_1__v
    eqn_9_1__w_s = eqn_9_1__w_s

    @kwasak
    def eqn_9_1(self, A=None, rho_s=None, v=None, w_s=None):
        """
        w_s := motive steam flow rate, lb/hr
        v:= velocity
        A:= cross sectional area, ft^2
        rhos_s := motive steam density, lb/ft^3
        """
        return
