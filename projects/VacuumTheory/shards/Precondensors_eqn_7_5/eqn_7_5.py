from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_5__N_i_cap import eqn_7_5__N_i
from .eqn_7_5__N_nc_cap import eqn_7_5__N_nc
from .eqn_7_5__P_cap import eqn_7_5__P
from .eqn_7_5__P_c_cap import eqn_7_5__P_c
from .eqn_7_5__p_i import eqn_7_5__p_i

class Precondensors:
    eqn_7_5__N_i = staticmethod(eqn_7_5__N_i)
    eqn_7_5__N_nc = staticmethod(eqn_7_5__N_nc)
    eqn_7_5__P = staticmethod(eqn_7_5__P)
    eqn_7_5__P_c = staticmethod(eqn_7_5__P_c)
    eqn_7_5__p_i = staticmethod(eqn_7_5__p_i)

    @kwasak_static
    def eqn_7_5(N_i=None, N_nc=None, P=None, P_c=None, p_i=None, **kwargs):
        return
