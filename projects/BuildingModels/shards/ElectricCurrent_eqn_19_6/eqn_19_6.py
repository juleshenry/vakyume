from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_19_6__R_cap1 import eqn_19_6__R1
from .eqn_19_6__R_cap2 import eqn_19_6__R2
from .eqn_19_6__R_capeff import eqn_19_6__Reff


class ElectricCurrent:
    eqn_19_6__R1 = eqn_19_6__R1
    eqn_19_6__R2 = eqn_19_6__R2
    eqn_19_6__Reff = eqn_19_6__Reff

    @kwasak
    def eqn_19_6(self, R1=None, R2=None, Reff=None):
        """
        Reff := effective resistance
        """
        return
