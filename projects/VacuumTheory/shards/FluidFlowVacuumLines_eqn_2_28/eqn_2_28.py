from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_28__C_cap import eqn_2_28__C
from .eqn_2_28__D_cap import eqn_2_28__D
from .eqn_2_28__L_cap import eqn_2_28__L
from .eqn_2_28__P_p_cap import eqn_2_28__P_p
from .eqn_2_28__mu import eqn_2_28__mu

class FluidFlowVacuumLines:
    eqn_2_28__C = staticmethod(eqn_2_28__C)
    eqn_2_28__D = staticmethod(eqn_2_28__D)
    eqn_2_28__L = staticmethod(eqn_2_28__L)
    eqn_2_28__P_p = staticmethod(eqn_2_28__P_p)
    eqn_2_28__mu = staticmethod(eqn_2_28__mu)

    @kwasak_static
    def eqn_2_28(C=None, D=None, L=None, P_p=None, mu=None, **kwargs):
        return
