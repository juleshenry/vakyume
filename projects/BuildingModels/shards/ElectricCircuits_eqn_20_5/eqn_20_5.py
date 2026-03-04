from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_20_5__C_cap import eqn_20_5__C
from .eqn_20_5__I_capR_cap import eqn_20_5__IR
from .eqn_20_5__Q_cap import eqn_20_5__Q
from .eqn_20_5__V_cap import eqn_20_5__V


class ElectricCircuits:
    eqn_20_5__C = eqn_20_5__C
    eqn_20_5__IR = eqn_20_5__IR
    eqn_20_5__Q = eqn_20_5__Q
    eqn_20_5__V = eqn_20_5__V

    @kwasak
    def eqn_20_5(self, C=None, IR=None, Q=None, V=None):
        """
        V := voltage
        I := current
        R := resistance
        Q := charge
        C := capacitance
        """
        return
