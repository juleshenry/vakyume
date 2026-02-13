from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_12(
        A: float = None,
        Q_condensor_heat_duty: float = None,
        U: float = None,
        del_T: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_7_12__A(Q_condensor_heat_duty: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty / (U * del_T)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_12__Q_condensor_heat_duty(A: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A * U * del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty / (A * del_T)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty / (A * U)
        result.append(del_T)
        return result


