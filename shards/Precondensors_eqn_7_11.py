from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_11(
        Q_condensor_heat_duty: float = None,
        U_v: float = None,
        V_c: float = None,
        del_T_LM: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_7_11__Q_condensor_heat_duty(U_v: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v * V_c * del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_11__U_v(Q_condensor_heat_duty: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty / (V_c * del_T_LM)
        result.append(U_v)
        return result

    @staticmethod
    def eqn_7_11__V_c(Q_condensor_heat_duty: float, U_v: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result.append(V_c)
        return result

    @staticmethod
    def eqn_7_11__del_T_LM(Q_condensor_heat_duty: float, U_v: float, V_c: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty / (U_v * V_c)
        result.append(del_T_LM)
        return result


