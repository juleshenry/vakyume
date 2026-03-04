from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_15_8__R_cap import eqn_15_8__R
from .eqn_15_8__R_cap1 import eqn_15_8__R1
from .eqn_15_8__R_cap2 import eqn_15_8__R2


class FluidMechanics:
    eqn_15_8__R = eqn_15_8__R
    eqn_15_8__R1 = eqn_15_8__R1
    eqn_15_8__R2 = eqn_15_8__R2

    @kwasak
    def eqn_15_8(self, R=None, R1=None, R2=None):
        """
        R := effective resistance
        R1 := resistance of house 1
        R2 := resistance of house 2
        P := pressure difference
        Q := flow rate
        """
        return
