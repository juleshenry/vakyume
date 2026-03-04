from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_17a__L import eqn_2_17a__L
from .eqn_2_17a__d import eqn_2_17a__d
from .eqn_2_17a__delta_P import eqn_2_17a__delta_P
from .eqn_2_17a__mu import eqn_2_17a__mu
from .eqn_2_17a__v import eqn_2_17a__v


class FluidFlowVacuumLines:
    eqn_2_17a__L = eqn_2_17a__L
    eqn_2_17a__d = eqn_2_17a__d
    eqn_2_17a__delta_P = eqn_2_17a__delta_P
    eqn_2_17a__mu = eqn_2_17a__mu
    eqn_2_17a__v = eqn_2_17a__v

    @kwasak
    def eqn_2_17a(self, L=None, d=None, delta_P=None, mu=None, v=None):
        return
