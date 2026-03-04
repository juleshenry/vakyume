from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_8_1__NC import eqn_8_1__NC
from .eqn_8_1__NS import eqn_8_1__NS
from .eqn_8_1__SCON import eqn_8_1__SCON
from .eqn_8_1__installation_cost import eqn_8_1__installation_cost


class SelectingPump:
    eqn_8_1__NC = eqn_8_1__NC
    eqn_8_1__NS = eqn_8_1__NS
    eqn_8_1__SCON = eqn_8_1__SCON
    eqn_8_1__installation_cost = eqn_8_1__installation_cost

    @kwasak
    def eqn_8_1(self, NC=None, NS=None, SCON=None, installation_cost=None):
        """
        NS:= number ejector stages
        NC:= number of condensors
        SCON:=steam consumption based on 100-psig motive steam, lb/hr
        """
        return
