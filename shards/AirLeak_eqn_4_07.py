from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class AirLeak:
    @kwasak_static
    def eqn_4_07(
        W: float = None,
        W_T: float = None,
        sum_individual_leak_rates: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_4_07__W(W_T: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return result

    @staticmethod
    def eqn_4_07__W_T(W: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return result

    @staticmethod
    def eqn_4_07__sum_individual_leak_rates(W: float, W_T: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result


