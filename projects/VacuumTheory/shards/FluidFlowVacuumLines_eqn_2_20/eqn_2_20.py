from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_20__L_cap import eqn_2_20__L
from .eqn_2_20__sum_equivalent_length import eqn_2_20__sum_equivalent_length
from .eqn_2_20__sum_pipe import eqn_2_20__sum_pipe

class FluidFlowVacuumLines:
    eqn_2_20__L = eqn_2_20__L
    eqn_2_20__sum_equivalent_length = eqn_2_20__sum_equivalent_length
    eqn_2_20__sum_pipe = eqn_2_20__sum_pipe

    @kwasak
    def eqn_2_20(self, L=None, sum_equivalent_length=None, sum_pipe=None):
        """
        L:= laminar flow
        """
        return
