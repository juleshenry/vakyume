from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_3_15__t import eqn_3_15__t
from .eqn_3_15__v import eqn_3_15__v
from .eqn_3_15__v_A_cap import eqn_3_15__v_A
from .eqn_3_15__v_B_cap import eqn_3_15__v_B


class DescribingMotionInOneDimension:
    eqn_3_15__t = eqn_3_15__t
    eqn_3_15__v = eqn_3_15__v
    eqn_3_15__v_A = eqn_3_15__v_A
    eqn_3_15__v_B = eqn_3_15__v_B

    @kwasak
    def eqn_3_15(self, t=None, v=None, v_A=None, v_B=None):
        """
        v := velocity of the passenger
        v_B := velocity of the boat
        v_A := velocity of the passenger
        t := time
        """
        return
