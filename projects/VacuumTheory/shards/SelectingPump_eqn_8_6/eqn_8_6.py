from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_8_6__M_cap import eqn_8_6__M
from .eqn_8_6__P_1_cap import eqn_8_6__P_1
from .eqn_8_6__P_2_cap import eqn_8_6__P_2
from .eqn_8_6__R_cap import eqn_8_6__R
from .eqn_8_6__T_cap import eqn_8_6__T
from .eqn_8_6__adiabatic_hp import eqn_8_6__adiabatic_hp
from .eqn_8_6__k import eqn_8_6__k
from .eqn_8_6__w import eqn_8_6__w

class SelectingPump:
    eqn_8_6__M = staticmethod(eqn_8_6__M)
    eqn_8_6__P_1 = staticmethod(eqn_8_6__P_1)
    eqn_8_6__P_2 = staticmethod(eqn_8_6__P_2)
    eqn_8_6__R = staticmethod(eqn_8_6__R)
    eqn_8_6__T = staticmethod(eqn_8_6__T)
    eqn_8_6__adiabatic_hp = staticmethod(eqn_8_6__adiabatic_hp)
    eqn_8_6__k = staticmethod(eqn_8_6__k)
    eqn_8_6__w = staticmethod(eqn_8_6__w)

    @kwasak_static
    def eqn_8_6(M=None, P_1=None, P_2=None, R=None, T=None, adiabatic_hp=None, k=None, w=None, **kwargs):
        return
