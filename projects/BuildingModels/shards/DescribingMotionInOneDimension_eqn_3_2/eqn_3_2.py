from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_2__ax import eqn_3_2__ax
from .eqn_3_2__t import eqn_3_2__t
from .eqn_3_2__v import eqn_3_2__v
from .eqn_3_2__v_0x import eqn_3_2__v_0x


class DescribingMotionInOneDimension:
    eqn_3_2__ax = eqn_3_2__ax
    eqn_3_2__t = eqn_3_2__t
    eqn_3_2__v = eqn_3_2__v
    eqn_3_2__v_0x = eqn_3_2__v_0x

    @kwasak
    def eqn_3_2(self, ax=None, t=None, v=None, v_0x=None):
        """
        v := velocity
        v_0x := initial velocity
        ax := acceleration
        t := time
        """
        return
