from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_22(P_s=None, Q_throughput=None, S_p=None, **kwargs):
        return

    @staticmethod
    def eqn_2_22__P_s(Q_throughput: float, S_p: float):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        P_s = Q_throughput/S_p
        result.append(P_s)
        return result

    @staticmethod
    def eqn_2_22__Q_throughput(P_s: float, S_p: float):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        Q_throughput = P_s*S_p
        result.append(Q_throughput)
        return result

    @staticmethod
    def eqn_2_22__S_p(P_s: float, Q_throughput: float):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        S_p = Q_throughput/P_s
        result.append(S_p)
        return result

