from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class FluidMechanics:
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

    def eqn_15_8__R(self, R1: float, R2: float, **kwargs):
        # R = R1 + R2
        result = []
        R = R1 + R2
        result.append(R)
        return result

    def eqn_15_8__R1(self, R: float, R2: float, **kwargs):
        # R = R1 + R2
        result = []
        R1 = R - R2
        result.append(R1)
        return result

    def eqn_15_8__R2(self, R: float, R1: float, **kwargs):
        # R = R1 + R2
        result = []
        R2 = R - R1
        result.append(R2)
        return result
