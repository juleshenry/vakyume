from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_7__T_cap import eqn_2_7__T
from .eqn_2_7__k import eqn_2_7__k
from .eqn_2_7__m import eqn_2_7__m
from .eqn_2_7__v_a import eqn_2_7__v_a

class FluidFlowVacuumLines:
    eqn_2_7__T = staticmethod(eqn_2_7__T)
    eqn_2_7__k = staticmethod(eqn_2_7__k)
    eqn_2_7__m = staticmethod(eqn_2_7__m)
    eqn_2_7__v_a = staticmethod(eqn_2_7__v_a)

    @kwasak_static
    def eqn_2_7(T=None, k=None, m=None, v_a=None, **kwargs):
        return
