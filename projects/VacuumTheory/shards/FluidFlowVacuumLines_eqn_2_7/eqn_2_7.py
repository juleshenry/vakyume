from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_7__T_cap import eqn_2_7__T
from .eqn_2_7__k import eqn_2_7__k
from .eqn_2_7__m import eqn_2_7__m
from .eqn_2_7__v_a import eqn_2_7__v_a


class FluidFlowVacuumLines:
    eqn_2_7__T = eqn_2_7__T
    eqn_2_7__k = eqn_2_7__k
    eqn_2_7__m = eqn_2_7__m
    eqn_2_7__v_a = eqn_2_7__v_a

    @kwasak
    def eqn_2_7(self, T=None, k=None, m=None, v_a=None):
        return
