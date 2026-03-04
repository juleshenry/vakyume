from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_5__D_cap import eqn_2_5__D
from .eqn_2_5__L_cap import eqn_2_5__L
from .eqn_2_5__delta_P_cap import eqn_2_5__delta_P
from .eqn_2_5__mu import eqn_2_5__mu
from .eqn_2_5__q import eqn_2_5__q

class FluidFlowVacuumLines:
    eqn_2_5__D = eqn_2_5__D
    eqn_2_5__L = eqn_2_5__L
    eqn_2_5__delta_P = eqn_2_5__delta_P
    eqn_2_5__mu = eqn_2_5__mu
    eqn_2_5__q = eqn_2_5__q

    @kwasak
    def eqn_2_5(self, D=None, L=None, delta_P=None, mu=None, q=None):
        """
        q:=volumetric flow cm^3/s
        D:= pipe diam.,cm
        delta_P := upstream-downstream pressure, dyne/cm^3
        L:=length, cm
        mu:= coef. of visco., poise
        """
        return
