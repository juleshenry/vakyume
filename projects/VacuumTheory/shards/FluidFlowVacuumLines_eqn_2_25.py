from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_25(C=None, P_1=None, P_2=None, Q_throughput=None, **kwargs):
        return

    @staticmethod
    def eqn_2_25__C(P_1: float, P_2: float, Q_throughput: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        C = Q_throughput/(P_1 - P_2)
        result.append(C)
        return result

    @staticmethod
    def eqn_2_25__P_1(C: float, P_2: float, Q_throughput: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        P_1 = P_2 + Q_throughput/C
        result.append(P_1)
        return result

    @staticmethod
    def eqn_2_25__P_2(C: float, P_1: float, Q_throughput: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        P_2 = P_1 - Q_throughput/C
        result.append(P_2)
        return result

    @staticmethod
    def eqn_2_25__Q_throughput(C: float, P_1: float, P_2: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        Q_throughput = C*(P_1 - P_2)
        result.append(Q_throughput)
        return result

