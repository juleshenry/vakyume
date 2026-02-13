from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_14a(
        A: float = None,
        Q_condensor_heat_duty: float = None,
        U: float = None,
        del_T_LM: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_7_14a__A(Q_condensor_heat_duty: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty / (U * del_T_LM)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14a__Q_condensor_heat_duty(A: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A * U * del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14a__U(A: float, Q_condensor_heat_duty: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty / (A * del_T_LM)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14a__del_T_LM(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty / (A * U)
        result.append(del_T_LM)
        return result


