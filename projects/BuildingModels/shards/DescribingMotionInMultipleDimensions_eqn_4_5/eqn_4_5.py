from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_4_5__a import eqn_4_5__a
from .eqn_4_5__t import eqn_4_5__t
from .eqn_4_5__v import eqn_4_5__v
from .eqn_4_5__v_0 import eqn_4_5__v_0


class DescribingMotionInMultipleDimensions:
    eqn_4_5__a = eqn_4_5__a
    eqn_4_5__t = eqn_4_5__t
    eqn_4_5__v = eqn_4_5__v
    eqn_4_5__v_0 = eqn_4_5__v_0

    @kwasak
    def eqn_4_5(self, a=None, t=None, v=None, v_0=None):
        """
        v := velocity
        v_0 := initial velocity
        a := acceleration
        t := time
        """
        return
