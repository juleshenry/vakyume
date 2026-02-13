from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_20(L=None, sum_equivalent_length=None, sum_pipe=None, **kwargs):
        return

    @staticmethod
    def eqn_2_20__L(sum_equivalent_length: float, sum_pipe: float):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        L = sum_equivalent_length + sum_pipe
        result.append(L)
        return result

    @staticmethod
    def eqn_2_20__sum_equivalent_length(L: float, sum_pipe: float):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        sum_equivalent_length = L - sum_pipe
        result.append(sum_equivalent_length)
        return result

    @staticmethod
    def eqn_2_20__sum_pipe(L: float, sum_equivalent_length: float):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        sum_pipe = L - sum_equivalent_length
        result.append(sum_pipe)
        return result

