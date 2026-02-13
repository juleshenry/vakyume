from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_18b(R_ll=None, h=None, w=None, **kwargs):
        return

    @staticmethod
    def eqn_2_18b__R_ll(h: float, w: float):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        R_ll = h*w/(2*(h + w))
        result.append(R_ll)
        return result

    @staticmethod
    def eqn_2_18b__h(R_ll: float, w: float):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        h = 2*R_ll*w/(-2*R_ll + w)
        result.append(h)
        return result

    @staticmethod
    def eqn_2_18b__w(R_ll: float, h: float):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        w = 2*R_ll*h/(-2*R_ll + h)
        result.append(w)
        return result

