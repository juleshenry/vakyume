from cmath import log, sqrt, exp
from math import e, pi
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
    eqn_2_28__C = eqn_2_28__C
    eqn_2_28__D = eqn_2_28__D
    eqn_2_28__L = eqn_2_28__L
    eqn_2_28__P_p = eqn_2_28__P_p
    eqn_2_28__mu = eqn_2_28__mu

    @kwasak_static
    def eqn_2_28(self, C=None, D=None, L=None, P_p=None, mu=None):
        return
