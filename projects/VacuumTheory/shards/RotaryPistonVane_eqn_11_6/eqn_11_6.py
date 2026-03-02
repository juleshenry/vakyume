from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_11_6__P_0_V_cap import eqn_11_6__P_0_V
from .eqn_11_6__P_D_cap import eqn_11_6__P_D
from .eqn_11_6__P_v_0_cap import eqn_11_6__P_v_0
from .eqn_11_6__S_B_cap import eqn_11_6__S_B
from .eqn_11_6__S_D_cap import eqn_11_6__S_D
from .eqn_11_6__p_b import eqn_11_6__p_b
from .eqn_11_6__p_g import eqn_11_6__p_g
from .eqn_11_6__p_v_max import eqn_11_6__p_v_max

class RotaryPistonVane:
    eqn_11_6__P_0_V = eqn_11_6__P_0_V
    eqn_11_6__P_D = eqn_11_6__P_D
    eqn_11_6__P_v_0 = eqn_11_6__P_v_0
    eqn_11_6__S_B = eqn_11_6__S_B
    eqn_11_6__S_D = eqn_11_6__S_D
    eqn_11_6__p_b = eqn_11_6__p_b
    eqn_11_6__p_g = eqn_11_6__p_g
    eqn_11_6__p_v_max = eqn_11_6__p_v_max

    @kwasak_static
    def eqn_11_6(self, P_0_V=None, P_D=None, P_v_0=None, S_B=None, S_D=None, p_b=None, p_g=None, p_v_max=None):
        return
