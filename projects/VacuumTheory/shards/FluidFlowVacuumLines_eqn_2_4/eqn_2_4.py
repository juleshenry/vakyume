from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_4___beta import eqn_2_4___beta
from .eqn_2_4__mu import eqn_2_4__mu
from .eqn_2_4__vel_grad import eqn_2_4__vel_grad


class FluidFlowVacuumLines:
    eqn_2_4___beta = eqn_2_4___beta
    eqn_2_4__mu = eqn_2_4__mu
    eqn_2_4__vel_grad = eqn_2_4__vel_grad

    @kwasak
    def eqn_2_4(self, _beta=None, mu=None, vel_grad=None):
        """
        mu:=coefficient of viscosity
        """
        return
