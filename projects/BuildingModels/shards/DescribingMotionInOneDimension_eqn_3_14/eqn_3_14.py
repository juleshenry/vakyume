from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_14__a import eqn_3_14__a
from .eqn_3_14__a_A_cap import eqn_3_14__a_A


class DescribingMotionInOneDimension:
    eqn_3_14__a = eqn_3_14__a
    eqn_3_14__a_A = eqn_3_14__a_A

    @kwasak
    def eqn_3_14(self, a=None, a_A=None):
        """
        a := relative acceleration
        a_A := acceleration of the passenger
        """
        return
