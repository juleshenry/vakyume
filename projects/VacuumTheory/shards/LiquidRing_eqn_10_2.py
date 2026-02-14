from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_2(PS=None, Q_gas=None, V=None, dP=None, dt=None, **kwargs):
        return

    @staticmethod
    def eqn_10_2__PS(Q_gas: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        PS = Q_gas - V*dP/dt
        result.append(PS)
        return result

    @staticmethod
    def eqn_10_2__Q_gas(PS: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        Q_gas = PS + V*dP/dt
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_2__V(PS: float, Q_gas: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        V = dt*(-PS + Q_gas)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_10_2__dP(PS: float, Q_gas: float, V: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dP = dt*(-PS + Q_gas)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_10_2__dt(PS: float, Q_gas: float, V: float, dP: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dt = -V*dP/(PS - Q_gas)
        result.append(dt)
        return result

