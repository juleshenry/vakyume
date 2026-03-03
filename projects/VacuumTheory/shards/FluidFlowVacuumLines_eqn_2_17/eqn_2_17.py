from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_17__L_cap import eqn_2_17__L
from .eqn_2_17__d import eqn_2_17__d
from .eqn_2_17__delta_P_cap import eqn_2_17__delta_P
from .eqn_2_17__mu import eqn_2_17__mu
from .eqn_2_17__v import eqn_2_17__v


class FluidFlowVacuumLines:
    eqn_2_17__L = eqn_2_17__L
    eqn_2_17__d = eqn_2_17__d
    eqn_2_17__delta_P = eqn_2_17__delta_P
    eqn_2_17__mu = eqn_2_17__mu
    eqn_2_17__v = eqn_2_17__v

    @kwasak
    def eqn_2_17(self, L=None, d=None, delta_P=None, mu=None, v=None):
        return
