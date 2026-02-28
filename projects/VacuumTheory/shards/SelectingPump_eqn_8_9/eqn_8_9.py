from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_8_9__E_j_cap import eqn_8_9__E_j
from .eqn_8_9__E_m_cap import eqn_8_9__E_m
from .eqn_8_9__e import eqn_8_9__e
from .eqn_8_9__r import eqn_8_9__r
from .eqn_8_9__s import eqn_8_9__s

class SelectingPump:
    eqn_8_9__E_j = staticmethod(eqn_8_9__E_j)
    eqn_8_9__E_m = staticmethod(eqn_8_9__E_m)
    eqn_8_9__e = staticmethod(eqn_8_9__e)
    eqn_8_9__r = staticmethod(eqn_8_9__r)
    eqn_8_9__s = staticmethod(eqn_8_9__s)

    @kwasak_static
    def eqn_8_9(E_j=None, E_m=None, e=None, r=None, s=None, **kwargs):
        return
