from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_04aa(
        n_i: float = None,
        n_nc: float = None,
        p_i: float = None,
        p_nc: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_7_04aa__n_i(n_nc: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_i = n_nc * p_i / p_nc
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04aa__n_nc(n_i: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_nc = n_i * p_nc / p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04aa__p_i(n_i: float, n_nc: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_i = n_i * p_nc / n_nc
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04aa__p_nc(n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_nc = n_nc * p_i / n_i
        result.append(p_nc)
        return result


