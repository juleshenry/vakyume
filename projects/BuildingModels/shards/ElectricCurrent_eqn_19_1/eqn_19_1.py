from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_19_1__I_cap import eqn_19_1__I
from .eqn_19_1__Q_cap import eqn_19_1__Q
from .eqn_19_1__t import eqn_19_1__t


class ElectricCurrent:
    eqn_19_1__I = eqn_19_1__I
    eqn_19_1__Q = eqn_19_1__Q
    eqn_19_1__t = eqn_19_1__t

    @kwasak
    def eqn_19_1(self, I=None, Q=None, t=None):
        """
        I := electric current
        Q := electric charge
        t := time
        """
        return
