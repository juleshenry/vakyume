from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_16__Re_cap import eqn_2_16__Re
from .eqn_2_16__f import eqn_2_16__f

class FluidFlowVacuumLines:
    eqn_2_16__Re = eqn_2_16__Re
    eqn_2_16__f = eqn_2_16__f

    @kwasak_static
    def eqn_2_16(self, Re=None, f=None):
        return
