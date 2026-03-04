from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_36__C_cap import eqn_2_36__C
from .eqn_2_36__C_cap_0 import eqn_2_36__C_0
from .eqn_2_36__F_cap_t import eqn_2_36__F_t

class FluidFlowVacuumLines:
    eqn_2_36__C = eqn_2_36__C
    eqn_2_36__C_0 = eqn_2_36__C_0
    eqn_2_36__F_t = eqn_2_36__F_t

    @kwasak
    def eqn_2_36(self, C=None, C_0=None, F_t=None):
        """
        C_0:=conductance thin walled aperture
        F_t:=transmission prob. for component
        """
        return
