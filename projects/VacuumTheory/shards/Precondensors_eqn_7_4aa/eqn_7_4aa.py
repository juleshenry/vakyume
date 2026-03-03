from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_4aa__n_i import eqn_7_4aa__n_i
from .eqn_7_4aa__n_nc import eqn_7_4aa__n_nc
from .eqn_7_4aa__p_i import eqn_7_4aa__p_i
from .eqn_7_4aa__p_nc import eqn_7_4aa__p_nc


class Precondensors:
    eqn_7_4aa__n_i = eqn_7_4aa__n_i
    eqn_7_4aa__n_nc = eqn_7_4aa__n_nc
    eqn_7_4aa__p_i = eqn_7_4aa__p_i
    eqn_7_4aa__p_nc = eqn_7_4aa__p_nc

    @kwasak
    def eqn_7_4aa(self, n_i=None, n_nc=None, p_i=None, p_nc=None):
        return
