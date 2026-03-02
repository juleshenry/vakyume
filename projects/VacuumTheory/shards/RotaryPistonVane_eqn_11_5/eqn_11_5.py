from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_11_5__P_0_v_cap import eqn_11_5__P_0_v
from .eqn_11_5__P_D_cap import eqn_11_5__P_D
from .eqn_11_5__p_g import eqn_11_5__p_g
from .eqn_11_5__p_v_max import eqn_11_5__p_v_max

class RotaryPistonVane:
    eqn_11_5__P_0_v = eqn_11_5__P_0_v
    eqn_11_5__P_D = eqn_11_5__P_D
    eqn_11_5__p_g = eqn_11_5__p_g
    eqn_11_5__p_v_max = eqn_11_5__p_v_max

    @kwasak_static
    def eqn_11_5(self, P_0_v=None, P_D=None, p_g=None, p_v_max=None):
        return
