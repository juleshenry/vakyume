from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_5_16__H_i_cap import eqn_5_16__H_i
from .eqn_5_16__p_i import eqn_5_16__p_i
from .eqn_5_16__x_i import eqn_5_16__x_i

class ProcessApp1:
    eqn_5_16__H_i = eqn_5_16__H_i
    eqn_5_16__p_i = eqn_5_16__p_i
    eqn_5_16__x_i = eqn_5_16__x_i

    @kwasak_static
    def eqn_5_16(self, H_i=None, p_i=None, x_i=None):
        return
