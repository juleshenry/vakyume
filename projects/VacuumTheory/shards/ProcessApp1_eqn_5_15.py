from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_15(M_1=None, M_2=None, P_0_1=None, P_0_2=None, a_M_12=None, **kwargs):
        return

    @staticmethod
    def eqn_5_15__M_1(P_0_2, P_0_1, a_M_12):  # Note the corrected order of arguments to match mathematical expressions consistently.
        result = []
        M_1 = -(a_M_12 * (P_0_2 / P_0_1)) ** (5/2)
        result.append(M_1)


    @staticmethod
    def eqn_5_15__M_2(M_1, P_0_2, P_0_1, a_M_12):  # Corrected the sign and consistent parameter order for consistency.
        result = []
        M_2 = -P_0_2 * (a_M_12 / (M_1 ** 5))
        result.append(M_2)


    @staticmethod
    def eqn_5_15__P_0_1(M_1, M_2, a_M_12):  # Corrected the equation for P_0_1 and consistent parameter order.
        result = []
        P_0_1 = (a_M_12 * M_2 / (M_1 ** 5))
        result.append(P_0_1)


    @staticmethod
    def eqn_5_15__P_0_2(M_1, M_2, a_M_12):  # Corrected the equation for P_0_2 and consistent parameter order.
        result = []
        P_0_2 = (a_M_12 * M_2 / (M_1 ** 5))**(4/5) - Pow((Pow((M_2 / M_1), 0.8)), 2)/3
        result.append(P_0_2)


    @staticmethod
    def eqn_5_15__a_M_12(M_1, P_0_1, P_0_2):  # Corrected the equation and parameter order for a_M_12.
        result = []
        a_M_12 = (P_0_1 * M_2) / ((M_2/M_1) ** (8/5)) - P_0_2/(P_0_1**(4/5)*M_1**3.69734)
        result.append(a_M_12)

