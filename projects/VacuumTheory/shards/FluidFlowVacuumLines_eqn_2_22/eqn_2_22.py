from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_22__P_cap_s import eqn_2_22__P_s
from .eqn_2_22__Q_cap_throughput import eqn_2_22__Q_throughput
from .eqn_2_22__S_cap_p import eqn_2_22__S_p


class FluidFlowVacuumLines:
    eqn_2_22__P_s = eqn_2_22__P_s
    eqn_2_22__Q_throughput = eqn_2_22__Q_throughput
    eqn_2_22__S_p = eqn_2_22__S_p

    @kwasak
    def eqn_2_22(self, P_s=None, Q_throughput=None, S_p=None):
        """
        Q:= through_put, sucking pressure P
        S_p:= dV / Dt
        """
        return
