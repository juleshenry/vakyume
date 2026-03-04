from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_4__a import eqn_3_4__a
from .eqn_3_4__t import eqn_3_4__t
from .eqn_3_4__v import eqn_3_4__v
from .eqn_3_4__v_0 import eqn_3_4__v_0


class DescribingMotionInOneDimension:
    eqn_3_4__a = eqn_3_4__a
    eqn_3_4__t = eqn_3_4__t
    eqn_3_4__v = eqn_3_4__v
    eqn_3_4__v_0 = eqn_3_4__v_0

    @kwasak
    def eqn_3_4(self, a=None, t=None, v=None, v_0=None):
        """
        x := position
        t := time
        v_0 := initial velocity
        a := acceleration
        """
        return
