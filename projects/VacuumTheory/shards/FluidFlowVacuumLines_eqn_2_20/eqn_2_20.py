from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_20__L_cap import eqn_2_20__L
from .eqn_2_20__sum_equivalent_length import eqn_2_20__sum_equivalent_length
from .eqn_2_20__sum_pipe import eqn_2_20__sum_pipe

class FluidFlowVacuumLines:
    eqn_2_20__L = eqn_2_20__L
    eqn_2_20__sum_equivalent_length = eqn_2_20__sum_equivalent_length
    eqn_2_20__sum_pipe = eqn_2_20__sum_pipe

    @kwasak_static
    def eqn_2_20(self, L=None, sum_equivalent_length=None, sum_pipe=None):
        return
