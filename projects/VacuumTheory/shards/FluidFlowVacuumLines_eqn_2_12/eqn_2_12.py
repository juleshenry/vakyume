from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_12__L_cap import eqn_2_12__L
from .eqn_2_12__d import eqn_2_12__d
from .eqn_2_12__delta_P_cap import eqn_2_12__delta_P
from .eqn_2_12__f import eqn_2_12__f
from .eqn_2_12__g import eqn_2_12__g
from .eqn_2_12__rho import eqn_2_12__rho
from .eqn_2_12__v import eqn_2_12__v

class FluidFlowVacuumLines:
    eqn_2_12__L = eqn_2_12__L
    eqn_2_12__d = eqn_2_12__d
    eqn_2_12__delta_P = eqn_2_12__delta_P
    eqn_2_12__f = eqn_2_12__f
    eqn_2_12__g = eqn_2_12__g
    eqn_2_12__rho = eqn_2_12__rho
    eqn_2_12__v = eqn_2_12__v

    @kwasak_static
    def eqn_2_12(self, L=None, d=None, delta_P=None, f=None, g=None, rho=None, v=None):
        return
