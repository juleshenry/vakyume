from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_7_6__F_cap_ import eqn_7_6__F_
from .eqn_7_6__W_cap_f import eqn_7_6__W_f
from .eqn_7_6__d import eqn_7_6__d


class WorkAndEnergy:
    eqn_7_6__F_ = eqn_7_6__F_
    eqn_7_6__W_f = eqn_7_6__W_f
    eqn_7_6__d = eqn_7_6__d

    @kwasak
    def eqn_7_6(self, F_=None, W_f=None, d=None):
        """
        W_f := work done by force of friction
        F_k := force of kinetic friction
        d := displacement
        """
        return
