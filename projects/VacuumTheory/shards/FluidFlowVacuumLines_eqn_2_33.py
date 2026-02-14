from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_33(C_paralell=None, arithmetic_sum_C=None, **kwargs):
        return

    @staticmethod
    def eqn_2_33__C_paralell(arithmetic_sum_C: float):
        # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
        result = []
        C_paralell = 1/arithmetic_sum_C
        result.append(C_paralell)
        return result

    @staticmethod
    def eqn_2_33__arithmetic_sum_C(C_paralell: float):
        # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
        result = []
        arithmetic_sum_C = 1/C_paralell
        result.append(arithmetic_sum_C)
        return result

