from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_4a__P_cap import eqn_7_4a__P
from .eqn_7_4a__p_c import eqn_7_4a__p_c
from .eqn_7_4a__p_nc import eqn_7_4a__p_nc

class Precondensors:
    eqn_7_4a__P = staticmethod(eqn_7_4a__P)
    eqn_7_4a__p_c = staticmethod(eqn_7_4a__p_c)
    eqn_7_4a__p_nc = staticmethod(eqn_7_4a__p_nc)

    @kwasak_static
    def eqn_7_4a(P=None, p_c=None, p_nc=None, **kwargs):
        return
