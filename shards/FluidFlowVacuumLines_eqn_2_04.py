from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_04(
        _beta: float = None, mu: float = None, vel_grad: float = None, **kwargs
    ):
        return

    @staticmethod
    def eqn_2_04___beta(mu: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu * vel_grad
        result.append(_beta)
        return result

    @staticmethod
    def eqn_2_04__mu(_beta: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta / vel_grad
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_04__vel_grad(_beta: float, mu: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta / mu
        result.append(vel_grad)
        return result


