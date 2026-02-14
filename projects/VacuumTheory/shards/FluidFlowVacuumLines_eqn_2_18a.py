from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_18a(D_eq=None, R_ll=None, **kwargs):
        return

    @staticmethod
    def eqn_2_18a__D_eq(R_ll: float, **kwargs):
        # [.pyeqn] D_eq = 4 * R_ll
        result = []
        D_eq = 4*R_ll
        result.append(D_eq)
        return result

    @staticmethod
    def eqn_2_18a__R_ll(D_eq: float, **kwargs):
        # [.pyeqn] D_eq = 4 * R_ll
        result = []
        R_ll = D_eq/4
        result.append(R_ll)
        return result

