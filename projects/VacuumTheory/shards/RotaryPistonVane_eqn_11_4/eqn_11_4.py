from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

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
        """
        p_v := partial pressure of vapor at pump suction, torr
        p_g := pressure of permanent gas at pump suction, torr
        p_s := pump suction pressure, sum of partial pressure of vapor and partial pressure of permanent gas, torr
        P_0_V := saturation pressure of vapor at pump operating temperature, torr
        P_D := pump discharge pressure, torr
        """
        return
