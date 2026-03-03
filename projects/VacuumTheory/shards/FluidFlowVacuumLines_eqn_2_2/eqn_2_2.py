from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_2__delta import eqn_2_2__delta
from .eqn_2_2__lambd import eqn_2_2__lambd
from .eqn_2_2__psi import eqn_2_2__psi


class FluidFlowVacuumLines:
    eqn_2_2__delta = eqn_2_2__delta
    eqn_2_2__lambd = eqn_2_2__lambd
    eqn_2_2__psi = eqn_2_2__psi

    @kwasak
    def eqn_2_2(self, delta=None, lambd=None, psi=None):
        """
        lambd := average mean free path , in
        delta := mol. diam , in
        psi:= mol. density molecules/in^3
        """
        return
