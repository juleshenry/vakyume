from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_18b__R_ll_cap import eqn_2_18b__R_ll
from .eqn_2_18b__h import eqn_2_18b__h
from .eqn_2_18b__w import eqn_2_18b__w

class FluidFlowVacuumLines:
    eqn_2_18b__R_ll = staticmethod(eqn_2_18b__R_ll)
    eqn_2_18b__h = staticmethod(eqn_2_18b__h)
    eqn_2_18b__w = staticmethod(eqn_2_18b__w)

    @kwasak_static
    def eqn_2_18b(R_ll=None, h=None, w=None, **kwargs):
        return
