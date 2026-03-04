from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_10_1__m import eqn_10_1__m
from .eqn_10_1__p import eqn_10_1__p
from .eqn_10_1__v import eqn_10_1__v


class LinearMomentumAndTheCentreOfMass:
    eqn_10_1__m = eqn_10_1__m
    eqn_10_1__p = eqn_10_1__p
    eqn_10_1__v = eqn_10_1__v

    @kwasak
    def eqn_10_1(self, m=None, p=None, v=None):
        """
        m := mass
        v := velocity
        p := momentum
        t := time
        a := acceleration
        F := force
        """
        return
