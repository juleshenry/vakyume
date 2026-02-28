from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_11_4__p_g import eqn_11_4__p_g
from .eqn_11_4__p_s import eqn_11_4__p_s
from .eqn_11_4__p_v import eqn_11_4__p_v

class RotaryPistonVane:
    eqn_11_4__p_g = staticmethod(eqn_11_4__p_g)
    eqn_11_4__p_s = staticmethod(eqn_11_4__p_s)
    eqn_11_4__p_v = staticmethod(eqn_11_4__p_v)

    @kwasak_static
    def eqn_11_4(p_g=None, p_s=None, p_v=None, **kwargs):
        return
