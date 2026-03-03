from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_4a__P_cap import eqn_7_4a__P
from .eqn_7_4a__p_c import eqn_7_4a__p_c
from .eqn_7_4a__p_nc import eqn_7_4a__p_nc


class Precondensors:
    eqn_7_4a__P = eqn_7_4a__P
    eqn_7_4a__p_c = eqn_7_4a__p_c
    eqn_7_4a__p_nc = eqn_7_4a__p_nc

    @kwasak
    def eqn_7_4a(self, P=None, p_c=None, p_nc=None):
        return
