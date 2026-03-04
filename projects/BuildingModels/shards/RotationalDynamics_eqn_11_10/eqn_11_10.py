from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_11_10__I_capC_capM_cap import eqn_11_10__ICM
from .eqn_11_10__I_caph import eqn_11_10__Ih
from .eqn_11_10__M_cap import eqn_11_10__M
from .eqn_11_10__h import eqn_11_10__h


class RotationalDynamics:
    eqn_11_10__ICM = eqn_11_10__ICM
    eqn_11_10__Ih = eqn_11_10__Ih
    eqn_11_10__M = eqn_11_10__M
    eqn_11_10__h = eqn_11_10__h

    @kwasak
    def eqn_11_10(self, ICM=None, Ih=None, M=None, h=None):
        """
        Ih := moment of inertia about a parallel axis
        ICM := moment of inertia about the centre of mass
        M := mass of the object
        h := distance from the centre of mass to the parallel axis
        """
        return
