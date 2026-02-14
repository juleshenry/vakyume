from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_3(N_mfw=None, Q_gas=None, T=None, **kwargs):
        return

    @staticmethod
    def eqn_10_3__N_mfw(Q_gas: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        N_mfw = 0.108108108108108*Q_gas/T
        result.append(N_mfw)
        return result

    @staticmethod
    def eqn_10_3__Q_gas(N_mfw: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        Q_gas = 9.25*N_mfw*T
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_3__T(N_mfw: float, Q_gas: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        T = 0.108108108108108*Q_gas/N_mfw
        result.append(T)
        return result

