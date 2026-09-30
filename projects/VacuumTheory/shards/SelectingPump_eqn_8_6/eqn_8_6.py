from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_8_6__M_cap import eqn_8_6__M
from .eqn_8_6__P_cap_1 import eqn_8_6__P_1
from .eqn_8_6__P_cap_2 import eqn_8_6__P_2
from .eqn_8_6__R_cap import eqn_8_6__R
from .eqn_8_6__T_cap import eqn_8_6__T
from .eqn_8_6__adiabatic_hp import eqn_8_6__adiabatic_hp
from .eqn_8_6__k import eqn_8_6__k
from .eqn_8_6__w import eqn_8_6__w


class SelectingPump:
    eqn_8_6__M = eqn_8_6__M
    eqn_8_6__P_1 = eqn_8_6__P_1
    eqn_8_6__P_2 = eqn_8_6__P_2
    eqn_8_6__R = eqn_8_6__R
    eqn_8_6__T = eqn_8_6__T
    eqn_8_6__adiabatic_hp = eqn_8_6__adiabatic_hp
    eqn_8_6__k = eqn_8_6__k
    eqn_8_6__w = eqn_8_6__w

    @kwasak
    def eqn_8_6(
        self,
        M=None,
        P_1=None,
        P_2=None,
        R=None,
        T=None,
        adiabatic_hp=None,
        k=None,
        w=None,
    ):
        """
        deg_R:=absolute temperature
        M:=molecular weight
        R:=gas constant, 1544 ft*lb_f / (lb*mol) * deg_R
        T:= absolute temperature, deg_R
        P:= absolute pressure, torr
        """
        return
