from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_8_7__P_cap_1 import eqn_8_7__P_1
from .eqn_8_7__P_cap_2 import eqn_8_7__P_2
from .eqn_8_7__adiabatic_hp import eqn_8_7__adiabatic_hp
from .eqn_8_7__w import eqn_8_7__w


class SelectingPump:
    eqn_8_7__P_1 = eqn_8_7__P_1
    eqn_8_7__P_2 = eqn_8_7__P_2
    eqn_8_7__adiabatic_hp = eqn_8_7__adiabatic_hp
    eqn_8_7__w = eqn_8_7__w

    @kwasak
    def eqn_8_7(self, P_1=None, P_2=None, adiabatic_hp=None, w=None):
        return
