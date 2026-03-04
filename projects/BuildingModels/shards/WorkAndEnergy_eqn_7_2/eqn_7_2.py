from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_7_2__F_cap1 import eqn_7_2__F1
from .eqn_7_2__F_cap2 import eqn_7_2__F2
from .eqn_7_2__F_cap3 import eqn_7_2__F3
from .eqn_7_2__W_cap_tot import eqn_7_2__W_tot
from .eqn_7_2__x import eqn_7_2__x


class WorkAndEnergy:
    eqn_7_2__F1 = eqn_7_2__F1
    eqn_7_2__F2 = eqn_7_2__F2
    eqn_7_2__F3 = eqn_7_2__F3
    eqn_7_2__W_tot = eqn_7_2__W_tot
    eqn_7_2__x = eqn_7_2__x

    @kwasak
    def eqn_7_2(self, F1=None, F2=None, F3=None, W_tot=None, x=None):
        """
        W_tot := total work
        """
        return
