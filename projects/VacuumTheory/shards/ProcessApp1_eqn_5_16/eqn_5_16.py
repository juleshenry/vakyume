from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_16__H_i_cap import eqn_5_16__H_i
from .eqn_5_16__p_i import eqn_5_16__p_i
from .eqn_5_16__x_i import eqn_5_16__x_i

class ProcessApp1:
    eqn_5_16__H_i = staticmethod(eqn_5_16__H_i)
    eqn_5_16__p_i = staticmethod(eqn_5_16__p_i)
    eqn_5_16__x_i = staticmethod(eqn_5_16__x_i)

    @kwasak_static
    def eqn_5_16(H_i=None, p_i=None, x_i=None, **kwargs):
        return
