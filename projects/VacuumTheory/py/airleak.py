from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class AirLeak:
    @kwasak
    def eqn_4_10(self, T=None, V=None, del_P=None, leakage=None, t=None):
        return

    def eqn_4_10__T(self, V: float, del_P: float, leakage: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        T = 3.127 * V * del_P / (leakage * t)
        result.append(T)
        return result

    def eqn_4_10__V(self, T: float, del_P: float, leakage: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        V = 0.319795330988168 * T * leakage * t / del_P
        result.append(V)
        return result

    def eqn_4_10__del_P(self, T: float, V: float, leakage: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        del_P = 0.319795330988168 * T * leakage * t / V
        result.append(del_P)
        return result

    def eqn_4_10__leakage(self, T: float, V: float, del_P: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        leakage = 3.127 * V * del_P / (T * t)
        result.append(leakage)
        return result

    def eqn_4_10__t(self, T: float, V: float, del_P: float, leakage: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        t = 3.127 * V * del_P / (T * leakage)
        result.append(t)
        return result

    @kwasak
    def eqn_4_7(self, W=None, W_T=None, sum_individual_leak_rates=None):
        return

    def eqn_4_7__W(self, W_T: float, sum_individual_leak_rates: float, **kwargs):
        # W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return result

    def eqn_4_7__W_T(self, W: float, sum_individual_leak_rates: float, **kwargs):
        # W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return result

    def eqn_4_7__sum_individual_leak_rates(self, W: float, W_T: float, **kwargs):
        # W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result
