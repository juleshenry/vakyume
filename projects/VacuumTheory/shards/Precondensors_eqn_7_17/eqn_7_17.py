from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_17__R_0_cap import eqn_7_17__R_0
from .eqn_7_17__R_nc_cap import eqn_7_17__R_nc
from .eqn_7_17__h_c import eqn_7_17__h_c


class Precondensors:
    eqn_7_17__R_0 = eqn_7_17__R_0
    eqn_7_17__R_nc = eqn_7_17__R_nc
    eqn_7_17__h_c = eqn_7_17__h_c

    @kwasak
    def eqn_7_17(self, R_0=None, R_nc=None, h_c=None):
        return
