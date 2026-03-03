from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_34__C_cap import eqn_2_34__C
from .eqn_2_34__C_1_cap import eqn_2_34__C_1
from .eqn_2_34__C_2_cap import eqn_2_34__C_2
from .eqn_2_34__D_cap import eqn_2_34__D
from .eqn_2_34__L_cap import eqn_2_34__L
from .eqn_2_34__P_p_cap import eqn_2_34__P_p
from .eqn_2_34__mu import eqn_2_34__mu

class FluidFlowVacuumLines:
    eqn_2_34__C = eqn_2_34__C
    eqn_2_34__C_1 = eqn_2_34__C_1
    eqn_2_34__C_2 = eqn_2_34__C_2
    eqn_2_34__D = eqn_2_34__D
    eqn_2_34__L = eqn_2_34__L
    eqn_2_34__P_p = eqn_2_34__P_p
    eqn_2_34__mu = eqn_2_34__mu

    @kwasak
    def eqn_2_34(self, C=None, C_1=None, C_2=None, D=None, L=None, P_p=None, mu=None):
        return
