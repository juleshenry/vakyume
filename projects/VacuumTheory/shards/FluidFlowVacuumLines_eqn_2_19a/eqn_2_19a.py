from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_19a__R_ll_cap import eqn_2_19a__R_ll
from .eqn_2_19a__Re_cap import eqn_2_19a__Re
from .eqn_2_19a__mu import eqn_2_19a__mu
from .eqn_2_19a__rho import eqn_2_19a__rho
from .eqn_2_19a__v import eqn_2_19a__v

class FluidFlowVacuumLines:
    eqn_2_19a__R_ll = eqn_2_19a__R_ll
    eqn_2_19a__Re = eqn_2_19a__Re
    eqn_2_19a__mu = eqn_2_19a__mu
    eqn_2_19a__rho = eqn_2_19a__rho
    eqn_2_19a__v = eqn_2_19a__v

    @kwasak
    def eqn_2_19a(self, R_ll=None, Re=None, mu=None, rho=None, v=None):
        return
