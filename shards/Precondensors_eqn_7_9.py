from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_9(L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, rho=None, **kwargs):
        return

    @staticmethod
    def eqn_7_9__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_9__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_9__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_9__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return result

    @staticmethod
    def eqn_7_9__rho(L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return result

