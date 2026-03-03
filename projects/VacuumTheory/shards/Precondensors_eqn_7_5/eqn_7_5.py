from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_5__N_cap_i import eqn_7_5__N_i
from .eqn_7_5__N_cap_nc import eqn_7_5__N_nc
from .eqn_7_5__P_cap import eqn_7_5__P
from .eqn_7_5__P_cap_c import eqn_7_5__P_c
from .eqn_7_5__p_i import eqn_7_5__p_i

class Precondensors:
    eqn_7_5__N_i = eqn_7_5__N_i
    eqn_7_5__N_nc = eqn_7_5__N_nc
    eqn_7_5__P = eqn_7_5__P
    eqn_7_5__P_c = eqn_7_5__P_c
    eqn_7_5__p_i = eqn_7_5__p_i

    @kwasak
    def eqn_7_5(self, N_i=None, N_nc=None, P=None, P_c=None, p_i=None):
        return
