from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_4ac__P_c import eqn_7_4ac__P_c
from .eqn_7_4ac__n_i import eqn_7_4ac__n_i
from .eqn_7_4ac__n_nc import eqn_7_4ac__n_nc
from .eqn_7_4ac__p import eqn_7_4ac__p
from .eqn_7_4ac__p_i import eqn_7_4ac__p_i


class Precondensors:
    eqn_7_4ac__P_c = eqn_7_4ac__P_c
    eqn_7_4ac__n_i = eqn_7_4ac__n_i
    eqn_7_4ac__n_nc = eqn_7_4ac__n_nc
    eqn_7_4ac__p = eqn_7_4ac__p
    eqn_7_4ac__p_i = eqn_7_4ac__p_i

    @kwasak
    def eqn_7_4ac(self, P_c=None, n_i=None, n_nc=None, p=None, p_i=None):
        return
