from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp2:
    @kwasak_static
    def eqn_6_06(
        Q_r: float = None, delta_h_v: float = None, f_m: float = None, **kwargs
    ):
        return

    @staticmethod
    def eqn_6_06__Q_r(delta_h_v: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        Q_r = delta_h_v * f_m / 24
        result.append(Q_r)
        return result

    @staticmethod
    def eqn_6_06__delta_h_v(Q_r: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        delta_h_v = 24 * Q_r / f_m
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_06__f_m(Q_r: float, delta_h_v: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        f_m = 24 * Q_r / delta_h_v
        result.append(f_m)
        return result


