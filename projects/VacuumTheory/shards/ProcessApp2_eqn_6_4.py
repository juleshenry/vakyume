from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp2:
    @kwasak_static
    def eqn_6_4(Q_v=None, delta_h_v=None, w_v=None, **kwargs):
        return

    @staticmethod
    def eqn_6_4__Q_v(delta_h_v: float, w_v: float, **kwargs):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        Q_v = delta_h_v*w_v/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_4__delta_h_v(Q_v: float, w_v: float, **kwargs):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        delta_h_v = 12000*Q_v/w_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_4__w_v(Q_v: float, delta_h_v: float, **kwargs):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return result

