from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_15__R_cape import eqn_2_15__Re
from .eqn_2_15__f import eqn_2_15__f

class FluidFlowVacuumLines:
    eqn_2_15__Re = eqn_2_15__Re
    eqn_2_15__f = eqn_2_15__f

    @kwasak
    def eqn_2_15(self, Re=None, f=None):
        return
