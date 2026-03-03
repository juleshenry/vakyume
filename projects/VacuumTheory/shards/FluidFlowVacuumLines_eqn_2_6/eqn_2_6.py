from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_6__lambd import eqn_2_6__lambd
from .eqn_2_6__mu import eqn_2_6__mu
from .eqn_2_6__rho import eqn_2_6__rho
from .eqn_2_6__v_a import eqn_2_6__v_a

class FluidFlowVacuumLines:
    eqn_2_6__lambd = eqn_2_6__lambd
    eqn_2_6__mu = eqn_2_6__mu
    eqn_2_6__rho = eqn_2_6__rho
    eqn_2_6__v_a = eqn_2_6__v_a

    @kwasak
    def eqn_2_6(self, lambd=None, mu=None, rho=None, v_a=None):
        """
        mu :=viscosity, poise
        rho:= density, g/cm^3
        lambd:= mean free path, cm
        """
        return
