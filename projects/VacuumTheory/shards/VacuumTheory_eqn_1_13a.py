from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class VacuumTheory:
    @kwasak_static
    def eqn_1_13a(n=None, n_a=None, y_a=None, **kwargs):
        return

    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float, **kwargs):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result

    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float, **kwargs):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result

    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float, **kwargs):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result

