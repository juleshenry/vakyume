from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_10(
        L_c_P: float = None,
        Q_condensor_heat_duty: float = None,
        del_T: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_7_10__L_c_P(Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result.append(L_c_P)
        return result

    @staticmethod
    def eqn_7_10__Q_condensor_heat_duty(L_c_P: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500 * L_c_P * del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_10__del_T(L_c_P: float, Q_condensor_heat_duty: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty / (500 * L_c_P)
        result.append(del_T)
        return result


