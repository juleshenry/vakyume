from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_4ab__P_c_cap import eqn_7_4ab__P_c
from .eqn_7_4ab__p import eqn_7_4ab__p
from .eqn_7_4ab__p_i import eqn_7_4ab__p_i
from .eqn_7_4ab__p_nc import eqn_7_4ab__p_nc

class Precondensors:
    eqn_7_4ab__P_c = staticmethod(eqn_7_4ab__P_c)
    eqn_7_4ab__p = staticmethod(eqn_7_4ab__p)
    eqn_7_4ab__p_i = staticmethod(eqn_7_4ab__p_i)
    eqn_7_4ab__p_nc = staticmethod(eqn_7_4ab__p_nc)

    @kwasak_static
    def eqn_7_4ab(P_c=None, p=None, p_i=None, p_nc=None, **kwargs):
        return
