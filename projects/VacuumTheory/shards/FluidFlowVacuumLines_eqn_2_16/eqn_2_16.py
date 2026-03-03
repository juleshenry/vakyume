from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_16__R_cape import eqn_2_16__Re
from .eqn_2_16__f import eqn_2_16__f

class FluidFlowVacuumLines:
    eqn_2_16__Re = eqn_2_16__Re
    eqn_2_16__f = eqn_2_16__f

    @kwasak
    def eqn_2_16(self, Re=None, f=None):
        return
