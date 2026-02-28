from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_17__R_0_cap import eqn_7_17__R_0
from .eqn_7_17__R_nc_cap import eqn_7_17__R_nc
from .eqn_7_17__h_c import eqn_7_17__h_c

class Precondensors:
    eqn_7_17__R_0 = staticmethod(eqn_7_17__R_0)
    eqn_7_17__R_nc = staticmethod(eqn_7_17__R_nc)
    eqn_7_17__h_c = staticmethod(eqn_7_17__h_c)

    @kwasak_static
    def eqn_7_17(R_0=None, R_nc=None, h_c=None, **kwargs):
        return
