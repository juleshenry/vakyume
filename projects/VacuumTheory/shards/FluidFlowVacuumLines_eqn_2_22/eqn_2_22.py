from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_22__P_s_cap import eqn_2_22__P_s
from .eqn_2_22__Q_throughput_cap import eqn_2_22__Q_throughput
from .eqn_2_22__S_p_cap import eqn_2_22__S_p

class FluidFlowVacuumLines:
    eqn_2_22__P_s = eqn_2_22__P_s
    eqn_2_22__Q_throughput = eqn_2_22__Q_throughput
    eqn_2_22__S_p = eqn_2_22__S_p

    @kwasak_static
    def eqn_2_22(self, P_s=None, Q_throughput=None, S_p=None):
        return
