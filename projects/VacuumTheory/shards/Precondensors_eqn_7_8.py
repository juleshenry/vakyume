from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_8(L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, **kwargs):
        return

    @staticmethod
    def eqn_7_8__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_8__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_8__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_8__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return result

