from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class RotaryPistonVane:
    @kwasak_static
    def eqn_11_01(
        PS: float = None,
        Q_0: float = None,
        Q_external_gas_throughput: float = None,
        V: float = None,
        dP: float = None,
        dT: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_11_01__PS(
        Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float
    ):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V * dP / dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_01__Q_0(
        PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float
    ):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V * dP / dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_01__Q_external_gas_throughput(
        PS: float, Q_0: float, V: float, dP: float, dT: float
    ):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V * dP / dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_01__V(
        PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float
    ):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT * (-PS + Q_0 + Q_external_gas_throughput) / dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_01__dP(
        PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float
    ):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT * (-PS + Q_0 + Q_external_gas_throughput) / V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_01__dT(
        PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float
    ):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V * dP / (-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result


