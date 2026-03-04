from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_4_4__a import eqn_4_4__a
from .eqn_4_4__ax import eqn_4_4__ax
from .eqn_4_4__ay import eqn_4_4__ay
from .eqn_4_4__x import eqn_4_4__x
from .eqn_4_4__y import eqn_4_4__y


class DescribingMotionInMultipleDimensions:
    eqn_4_4__a = eqn_4_4__a
    eqn_4_4__ax = eqn_4_4__ax
    eqn_4_4__ay = eqn_4_4__ay
    eqn_4_4__x = eqn_4_4__x
    eqn_4_4__y = eqn_4_4__y

    @kwasak
    def eqn_4_4(self, a=None, ax=None, ay=None, x=None, y=None):
        """
        a := acceleration
        x := x-component
        y := y-component
        t := time
        """
        return
