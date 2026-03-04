from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_6__t import eqn_3_6__t
from .eqn_3_6__v_0 import eqn_3_6__v_0
from .eqn_3_6__x import eqn_3_6__x
from .eqn_3_6__x_0 import eqn_3_6__x_0


class DescribingMotionInOneDimension:
    eqn_3_6__t = eqn_3_6__t
    eqn_3_6__v_0 = eqn_3_6__v_0
    eqn_3_6__x = eqn_3_6__x
    eqn_3_6__x_0 = eqn_3_6__x_0

    @kwasak
    def eqn_3_6(self, t=None, v_0=None, x=None, x_0=None):
        """
        x := position
        x_0 := initial position
        v_0 := initial velocity
        a := acceleration
        t := time
        """
        return
