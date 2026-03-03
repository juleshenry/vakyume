from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_1_12__Total_P_cap import eqn_1_12__Total_P
from .eqn_1_12__sum_partial_pressures import eqn_1_12__sum_partial_pressures

class VacuumTheory:
    eqn_1_12__Total_P = eqn_1_12__Total_P
    eqn_1_12__sum_partial_pressures = eqn_1_12__sum_partial_pressures

    @kwasak_static
    def eqn_1_12(self, Total_P=None, sum_partial_pressures=None):
        return
