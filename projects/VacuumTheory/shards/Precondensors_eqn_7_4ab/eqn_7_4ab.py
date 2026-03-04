from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_4ab__P_c_cap import eqn_7_4ab__P_c
from .eqn_7_4ab__p import eqn_7_4ab__p
from .eqn_7_4ab__p_i import eqn_7_4ab__p_i
from .eqn_7_4ab__p_nc import eqn_7_4ab__p_nc


class Precondensors:
    eqn_7_4ab__P_c = eqn_7_4ab__P_c
    eqn_7_4ab__p = eqn_7_4ab__p
    eqn_7_4ab__p_i = eqn_7_4ab__p_i
    eqn_7_4ab__p_nc = eqn_7_4ab__p_nc

    @kwasak
    def eqn_7_4ab(self, P_c=None, p=None, p_i=None, p_nc=None):
        return
