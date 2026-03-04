from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_2_25__C_cap import eqn_2_25__C
from .eqn_2_25__P_cap_1 import eqn_2_25__P_1
from .eqn_2_25__P_cap_2 import eqn_2_25__P_2
from .eqn_2_25__Q_cap_throughput import eqn_2_25__Q_throughput


class FluidFlowVacuumLines:
    eqn_2_25__C = eqn_2_25__C
    eqn_2_25__P_1 = eqn_2_25__P_1
    eqn_2_25__P_2 = eqn_2_25__P_2
    eqn_2_25__Q_throughput = eqn_2_25__Q_throughput

    @kwasak
    def eqn_2_25(self, C=None, P_1=None, P_2=None, Q_throughput=None):
        return
