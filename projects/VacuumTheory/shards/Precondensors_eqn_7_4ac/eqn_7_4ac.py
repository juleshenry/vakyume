from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_4ac__P_c_cap import eqn_7_4ac__P_c
from .eqn_7_4ac__n_i import eqn_7_4ac__n_i
from .eqn_7_4ac__n_nc import eqn_7_4ac__n_nc
from .eqn_7_4ac__p import eqn_7_4ac__p
from .eqn_7_4ac__p_i import eqn_7_4ac__p_i

class Precondensors:
    eqn_7_4ac__P_c = staticmethod(eqn_7_4ac__P_c)
    eqn_7_4ac__n_i = staticmethod(eqn_7_4ac__n_i)
    eqn_7_4ac__n_nc = staticmethod(eqn_7_4ac__n_nc)
    eqn_7_4ac__p = staticmethod(eqn_7_4ac__p)
    eqn_7_4ac__p_i = staticmethod(eqn_7_4ac__p_i)

    @kwasak_static
    def eqn_7_4ac(P_c=None, n_i=None, n_nc=None, p=None, p_i=None, **kwargs):
        return
