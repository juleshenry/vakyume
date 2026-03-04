from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_10_11__F_cap_ext import eqn_10_11__F_ext
from .eqn_10_11__M_cap import eqn_10_11__M
from .eqn_10_11__a_C_capM_cap import eqn_10_11__a_CM


class LinearMomentumAndTheCentreOfMass:
    eqn_10_11__F_ext = eqn_10_11__F_ext
    eqn_10_11__M = eqn_10_11__M
    eqn_10_11__a_CM = eqn_10_11__a_CM

    @kwasak
    def eqn_10_11(self, F_ext=None, M=None, a_CM=None):
        """
        F_ext := external force
        M := total mass
        a_CM := acceleration of the centre of mass
        """
        return
