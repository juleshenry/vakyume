from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_29__C_cap import eqn_2_29__C
from .eqn_2_29__S_cap_1 import eqn_2_29__S_1
from .eqn_2_29__S_cap_2 import eqn_2_29__S_2


class FluidFlowVacuumLines:
    eqn_2_29__C = eqn_2_29__C
    eqn_2_29__S_1 = eqn_2_29__S_1
    eqn_2_29__S_2 = eqn_2_29__S_2

    @kwasak
    def eqn_2_29(self, C=None, S_1=None, S_2=None):
        return
