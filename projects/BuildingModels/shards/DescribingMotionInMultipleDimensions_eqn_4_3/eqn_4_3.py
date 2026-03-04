from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_4_3__v import eqn_4_3__v
from .eqn_4_3__vx import eqn_4_3__vx
from .eqn_4_3__vy import eqn_4_3__vy
from .eqn_4_3__x import eqn_4_3__x
from .eqn_4_3__y import eqn_4_3__y


class DescribingMotionInMultipleDimensions:
    eqn_4_3__v = eqn_4_3__v
    eqn_4_3__vx = eqn_4_3__vx
    eqn_4_3__vy = eqn_4_3__vy
    eqn_4_3__x = eqn_4_3__x
    eqn_4_3__y = eqn_4_3__y

    @kwasak
    def eqn_4_3(self, v=None, vx=None, vy=None, x=None, y=None):
        """
        v := velocity
        x := x-component
        y := y-component
        t := time
        """
        return
