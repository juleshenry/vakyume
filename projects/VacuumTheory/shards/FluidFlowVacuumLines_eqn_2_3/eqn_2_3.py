from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_3__D_cap import eqn_2_3__D
from .eqn_2_3__kn import eqn_2_3__kn
from .eqn_2_3__lambd import eqn_2_3__lambd

class FluidFlowVacuumLines:
    eqn_2_3__D = eqn_2_3__D
    eqn_2_3__kn = eqn_2_3__kn
    eqn_2_3__lambd = eqn_2_3__lambd

    @kwasak_static
    def eqn_2_3(self, D=None, kn=None, lambd=None):
        return
