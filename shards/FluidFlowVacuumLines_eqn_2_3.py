from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_3(D=None, kn=None, lambd=None, **kwargs):
        return

    @staticmethod
    def eqn_2_3__D(kn: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return result

    @staticmethod
    def eqn_2_3__kn(D: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return result

    @staticmethod
    def eqn_2_3__lambd(D: float, kn: float):
        # [.pyeqn] kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return result

