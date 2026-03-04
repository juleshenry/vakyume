from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_4_6__a import eqn_4_6__a
from .eqn_4_6__r import eqn_4_6__r
from .eqn_4_6__r_0 import eqn_4_6__r_0
from .eqn_4_6__t import eqn_4_6__t
from .eqn_4_6__v_0 import eqn_4_6__v_0


class DescribingMotionInMultipleDimensions:
    eqn_4_6__a = eqn_4_6__a
    eqn_4_6__r = eqn_4_6__r
    eqn_4_6__r_0 = eqn_4_6__r_0
    eqn_4_6__t = eqn_4_6__t
    eqn_4_6__v_0 = eqn_4_6__v_0

    @kwasak
    def eqn_4_6(self, a=None, r=None, r_0=None, t=None, v_0=None):
        """
        r := position
        r_0 := initial position
        v_0 := initial velocity
        a := acceleration
        t := time
        """
        return
