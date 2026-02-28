from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_2_25__C_cap import eqn_2_25__C
from .eqn_2_25__P_1_cap import eqn_2_25__P_1
from .eqn_2_25__P_2_cap import eqn_2_25__P_2
from .eqn_2_25__Q_throughput_cap import eqn_2_25__Q_throughput

class FluidFlowVacuumLines:
    eqn_2_25__C = staticmethod(eqn_2_25__C)
    eqn_2_25__P_1 = staticmethod(eqn_2_25__P_1)
    eqn_2_25__P_2 = staticmethod(eqn_2_25__P_2)
    eqn_2_25__Q_throughput = staticmethod(eqn_2_25__Q_throughput)

    @kwasak_static
    def eqn_2_25(C=None, P_1=None, P_2=None, Q_throughput=None, **kwargs):
        return
