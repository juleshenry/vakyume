from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_8(P_0_1=None, P_0_2=None, alpha_12=None, gamma_1=None, gamma_2=None, **kwargs):
        return

    @staticmethod
    def eqn_5_8__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_8__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_8__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result

    @staticmethod
    def eqn_5_8__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return result

    @staticmethod
    def eqn_5_8__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return result

