from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_28__C import eqn_2_28__C
from .eqn_2_28__D import eqn_2_28__D
from .eqn_2_28__L import eqn_2_28__L
from .eqn_2_28__P_p import eqn_2_28__P_p
from .eqn_2_28__mu import eqn_2_28__mu

class FluidFlowVacuumLines:
    eqn_2_28__C = eqn_2_28__C
    eqn_2_28__D = eqn_2_28__D
    eqn_2_28__L = eqn_2_28__L
    eqn_2_28__P_p = eqn_2_28__P_p
    eqn_2_28__mu = eqn_2_28__mu

    @kwasak
    def eqn_2_28(self, C=None, D=None, L=None, P_p=None, mu=None):
        return
