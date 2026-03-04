from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_3__ax import eqn_3_3__ax
from .eqn_3_3__t import eqn_3_3__t
from .eqn_3_3__v_0x import eqn_3_3__v_0x
from .eqn_3_3__x import eqn_3_3__x
from .eqn_3_3__x_0 import eqn_3_3__x_0


class DescribingMotionInOneDimension:
    eqn_3_3__ax = eqn_3_3__ax
    eqn_3_3__t = eqn_3_3__t
    eqn_3_3__v_0x = eqn_3_3__v_0x
    eqn_3_3__x = eqn_3_3__x
    eqn_3_3__x_0 = eqn_3_3__x_0

    @kwasak
    def eqn_3_3(self, ax=None, t=None, v_0x=None, x=None, x_0=None):
        """
        x := position
        x_0 := initial position
        v_0x := initial velocity
        ax := acceleration
        t := time
        """
        return
