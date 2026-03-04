from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_20_3__I_cap1 import eqn_20_3__I1
from .eqn_20_3__I_cap2 import eqn_20_3__I2
from .eqn_20_3__I_cap3 import eqn_20_3__I3


class ElectricCircuits:
    eqn_20_3__I1 = eqn_20_3__I1
    eqn_20_3__I2 = eqn_20_3__I2
    eqn_20_3__I3 = eqn_20_3__I3

    @kwasak
    def eqn_20_3(self, I1=None, I2=None, I3=None):
        """
        I1 := current 1
        I2 := current 2
        I3 := current 3
        """
        return
