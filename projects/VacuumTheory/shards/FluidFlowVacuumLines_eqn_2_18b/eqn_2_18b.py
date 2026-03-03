from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_18b__R_cap_ll import eqn_2_18b__R_ll
from .eqn_2_18b__h import eqn_2_18b__h
from .eqn_2_18b__w import eqn_2_18b__w

class FluidFlowVacuumLines:
    eqn_2_18b__R_ll = eqn_2_18b__R_ll
    eqn_2_18b__h = eqn_2_18b__h
    eqn_2_18b__w = eqn_2_18b__w

    @kwasak
    def eqn_2_18b(self, R_ll=None, h=None, w=None):
        return
