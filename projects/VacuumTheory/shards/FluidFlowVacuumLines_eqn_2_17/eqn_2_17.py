from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_17__L_cap import eqn_2_17__L
from .eqn_2_17__d import eqn_2_17__d
from .eqn_2_17__delta_P_cap import eqn_2_17__delta_P
from .eqn_2_17__mu import eqn_2_17__mu
from .eqn_2_17__v import eqn_2_17__v

class FluidFlowVacuumLines:
    eqn_2_17__L = staticmethod(eqn_2_17__L)
    eqn_2_17__d = staticmethod(eqn_2_17__d)
    eqn_2_17__delta_P = staticmethod(eqn_2_17__delta_P)
    eqn_2_17__mu = staticmethod(eqn_2_17__mu)
    eqn_2_17__v = staticmethod(eqn_2_17__v)

    @kwasak_static
    def eqn_2_17(L=None, d=None, delta_P=None, mu=None, v=None, **kwargs):
        return
