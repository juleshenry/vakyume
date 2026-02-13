from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_15(U=None, sum_R=None, **kwargs):
        return

    @staticmethod
    def eqn_7_15__U(sum_R: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        U = 1/sum_R
        result.append(U)
        return result

    @staticmethod
    def eqn_7_15__sum_R(U: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        sum_R = 1/U
        result.append(sum_R)
        return result

