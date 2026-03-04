from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_12_11__I_cap import eqn_12_11__I
from .eqn_12_11__L_cap import eqn_12_11__L
from .eqn_12_11__v import eqn_12_11__v


class RotationalEnergyAndMomentum:
    eqn_12_11__I = eqn_12_11__I
    eqn_12_11__L = eqn_12_11__L
    eqn_12_11__v = eqn_12_11__v

    @kwasak
    def eqn_12_11(self, I=None, L=None, v=None):
        """
        m := mass of a particle
        r := position of the particle
        v := velocity of the particle
        I := moment of inertia
        """
        return
