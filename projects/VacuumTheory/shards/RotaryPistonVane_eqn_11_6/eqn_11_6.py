from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_11_6__P_cap_0_V_cap import eqn_11_6__P_0_V
from .eqn_11_6__P_cap_D_cap import eqn_11_6__P_D
from .eqn_11_6__P_cap_v_0 import eqn_11_6__P_v_0
from .eqn_11_6__S_cap_B_cap import eqn_11_6__S_B
from .eqn_11_6__S_cap_D_cap import eqn_11_6__S_D
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

    @kwasak
    def eqn_11_6(
        self,
        P_0_V=None,
        P_D=None,
        P_v_0=None,
        S_B=None,
        S_D=None,
        p_b=None,
        p_g=None,
        p_v_max=None,
    ):
        """
        P_0_v := saturation vapor pressure of a condensable vapor
        S_B := maximum permissible gas ballast flow rate, ft^3/min
        S_D := free air displacement of the vacuum pump, ft^3/min
        p_b := partial pressure of vapor in the ballast gas, e.g. partial pressure of water vapor in ATM, torr
        """
        return
