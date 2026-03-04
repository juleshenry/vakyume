from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_11_4__p_g import eqn_11_4__p_g
from .eqn_11_4__p_s import eqn_11_4__p_s
from .eqn_11_4__p_v import eqn_11_4__p_v


class RotaryPistonVane:
    eqn_11_4__p_g = eqn_11_4__p_g
    eqn_11_4__p_s = eqn_11_4__p_s
    eqn_11_4__p_v = eqn_11_4__p_v

    @kwasak
    def eqn_11_4(self, p_g=None, p_s=None, p_v=None):
        return
