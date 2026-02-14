from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class RotaryPistonVane:
    @kwasak_static
    def eqn_11_1(PS=None, Q_0=None, Q_external_gas_throughput=None, V=None, dP=None, dT=None, **kwargs):
        return

    @staticmethod
    def eqn_11_1__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_1__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_1__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_1__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_1__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_1__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

