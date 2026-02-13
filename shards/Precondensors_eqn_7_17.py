from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_17(R_0: float = None, R_nc: float = None, h_c: float = None, **kwargs):
        return

    @staticmethod
    def eqn_7_17__R_0(R_nc: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1 / h_c
        result.append(R_0)
        return result

    @staticmethod
    def eqn_7_17__R_nc(R_0: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1 / h_c
        result.append(R_nc)
        return result

    @staticmethod
    def eqn_7_17__h_c(R_0: float, R_nc: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1 / (R_0 - R_nc)
        result.append(h_c)
        return result


