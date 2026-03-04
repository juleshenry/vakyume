from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_2_35__C_L import eqn_2_35__C_L
from .eqn_2_35__C_T import eqn_2_35__C_T
from .eqn_2_35__F_p import eqn_2_35__F_p

class FluidFlowVacuumLines:
    eqn_2_35__C_L = eqn_2_35__C_L
    eqn_2_35__C_T = eqn_2_35__C_T
    eqn_2_35__F_p = eqn_2_35__F_p

    @kwasak
    def eqn_2_35(self, C_L=None, C_T=None, F_p=None):
        """
        F_P:= correction factor for Poiseuille's eqn from Figure
        """
        return
