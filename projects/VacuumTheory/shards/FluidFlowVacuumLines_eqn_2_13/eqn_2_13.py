from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_13__L_cap import eqn_2_13__L
from .eqn_2_13__d import eqn_2_13__d
from .eqn_2_13__delta_P_cap import eqn_2_13__delta_P
from .eqn_2_13__f import eqn_2_13__f
from .eqn_2_13__q import eqn_2_13__q
from .eqn_2_13__rho import eqn_2_13__rho


class FluidFlowVacuumLines:
    eqn_2_13__L = eqn_2_13__L
    eqn_2_13__d = eqn_2_13__d
    eqn_2_13__delta_P = eqn_2_13__delta_P
    eqn_2_13__f = eqn_2_13__f
    eqn_2_13__q = eqn_2_13__q
    eqn_2_13__rho = eqn_2_13__rho

    @kwasak
    def eqn_2_13(self, L=None, d=None, delta_P=None, f=None, q=None, rho=None):
        return
