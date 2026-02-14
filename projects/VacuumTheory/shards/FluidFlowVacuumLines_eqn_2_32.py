from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_32(C_series=None, geometric_sum_C=None, **kwargs):
        return

    @staticmethod
    def eqn_2_32__C_series(geometric_sum_C: float, **kwargs):
        # [.pyeqn] 1 / C_series = geometric_sum_C
        result = []
        C_series = 1/geometric_sum_C
        result.append(C_series)
        return result

    @staticmethod
    def eqn_2_32__geometric_sum_C(C_series: float, **kwargs):
        # [.pyeqn] 1 / C_series = geometric_sum_C
        result = []
        geometric_sum_C = 1/C_series
        result.append(geometric_sum_C)
        return result

