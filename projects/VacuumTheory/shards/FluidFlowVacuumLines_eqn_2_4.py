from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_4(_beta=None, mu=None, vel_grad=None, **kwargs):
        return

    @staticmethod
    def eqn_2_4___beta(mu: float, vel_grad: float, **kwargs):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return result

    @staticmethod
    def eqn_2_4__mu(_beta: float, vel_grad: float, **kwargs):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_4__vel_grad(_beta: float, mu: float, **kwargs):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return result

