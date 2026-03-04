from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_17b__L_cap import eqn_2_17b__L
from .eqn_2_17b__d import eqn_2_17b__d
from .eqn_2_17b__delta_P_cap import eqn_2_17b__delta_P
from .eqn_2_17b__mu import eqn_2_17b__mu
from .eqn_2_17b__q import eqn_2_17b__q


class FluidFlowVacuumLines:
    eqn_2_17b__L = eqn_2_17b__L
    eqn_2_17b__d = eqn_2_17b__d
    eqn_2_17b__delta_P = eqn_2_17b__delta_P
    eqn_2_17b__mu = eqn_2_17b__mu
    eqn_2_17b__q = eqn_2_17b__q

    @kwasak
    def eqn_2_17b(self, L=None, d=None, delta_P=None, mu=None, q=None):
        return
