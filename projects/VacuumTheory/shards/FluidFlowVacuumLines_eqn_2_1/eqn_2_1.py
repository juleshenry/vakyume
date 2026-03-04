from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_1__D_cap import eqn_2_1__D
from .eqn_2_1__Re_cap import eqn_2_1__Re
from .eqn_2_1__mu import eqn_2_1__mu
from .eqn_2_1__rho import eqn_2_1__rho
from .eqn_2_1__v import eqn_2_1__v


class FluidFlowVacuumLines:
    eqn_2_1__D = eqn_2_1__D
    eqn_2_1__Re = eqn_2_1__Re
    eqn_2_1__mu = eqn_2_1__mu
    eqn_2_1__rho = eqn_2_1__rho
    eqn_2_1__v = eqn_2_1__v

    @kwasak
    def eqn_2_1(self, D=None, Re=None, mu=None, rho=None, v=None):
        return
