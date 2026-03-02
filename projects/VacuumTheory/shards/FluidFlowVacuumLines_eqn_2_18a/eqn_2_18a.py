from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_18a__D_eq_cap import eqn_2_18a__D_eq
from .eqn_2_18a__R_ll_cap import eqn_2_18a__R_ll

class FluidFlowVacuumLines:
    eqn_2_18a__D_eq = eqn_2_18a__D_eq
    eqn_2_18a__R_ll = eqn_2_18a__R_ll

    @kwasak_static
    def eqn_2_18a(self, D_eq=None, R_ll=None):
        return
