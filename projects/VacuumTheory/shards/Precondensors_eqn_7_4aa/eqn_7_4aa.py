from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_4aa__n_i import eqn_7_4aa__n_i
from .eqn_7_4aa__n_nc import eqn_7_4aa__n_nc
from .eqn_7_4aa__p_i import eqn_7_4aa__p_i
from .eqn_7_4aa__p_nc import eqn_7_4aa__p_nc

class Precondensors:
    eqn_7_4aa__n_i = staticmethod(eqn_7_4aa__n_i)
    eqn_7_4aa__n_nc = staticmethod(eqn_7_4aa__n_nc)
    eqn_7_4aa__p_i = staticmethod(eqn_7_4aa__p_i)
    eqn_7_4aa__p_nc = staticmethod(eqn_7_4aa__p_nc)

    @kwasak_static
    def eqn_7_4aa(n_i=None, n_nc=None, p_i=None, p_nc=None, **kwargs):
        return
