from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_1__t import eqn_3_1__t
from .eqn_3_1__v_x import eqn_3_1__v_x
from .eqn_3_1__x import eqn_3_1__x
from .eqn_3_1__x_0 import eqn_3_1__x_0


class DescribingMotionInOneDimension:
    eqn_3_1__t = eqn_3_1__t
    eqn_3_1__v_x = eqn_3_1__v_x
    eqn_3_1__x = eqn_3_1__x
    eqn_3_1__x_0 = eqn_3_1__x_0

    @kwasak
    def eqn_3_1(self, t=None, v_x=None, x=None, x_0=None):
        """
        x := position
        x_0 := initial position
        v_x := velocity
        t := time
        """
        return
