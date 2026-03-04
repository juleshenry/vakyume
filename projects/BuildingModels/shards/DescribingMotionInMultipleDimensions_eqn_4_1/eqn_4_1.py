from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_4_1__t import eqn_4_1__t
from .eqn_4_1__v import eqn_4_1__v
from .eqn_4_1__x_1 import eqn_4_1__x_1
from .eqn_4_1__x_2 import eqn_4_1__x_2
from .eqn_4_1__y_1 import eqn_4_1__y_1
from .eqn_4_1__y_2 import eqn_4_1__y_2


class DescribingMotionInMultipleDimensions:
    eqn_4_1__t = eqn_4_1__t
    eqn_4_1__v = eqn_4_1__v
    eqn_4_1__x_1 = eqn_4_1__x_1
    eqn_4_1__x_2 = eqn_4_1__x_2
    eqn_4_1__y_1 = eqn_4_1__y_1
    eqn_4_1__y_2 = eqn_4_1__y_2

    @kwasak
    def eqn_4_1(self, t=None, v=None, x_1=None, x_2=None, y_1=None, y_2=None):
        """
        t := time
        x_1 := initial x-coordinate
        y_1 := initial y-coordinate
        x_2 := final x-coordinate
        y_2 := final y-coordinate
        """
        return
