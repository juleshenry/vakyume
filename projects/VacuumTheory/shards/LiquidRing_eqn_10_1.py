from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_1(D_r=None, sig_R=None, w=None, **kwargs):
        return

    @staticmethod
    def eqn_10_1__D_r(sig_R: float, w: float, **kwargs):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        D_r = 229.357798165138*sig_R/w
        result.append(D_r)
        return result

    @staticmethod
    def eqn_10_1__sig_R(D_r: float, w: float, **kwargs):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        sig_R = 0.00436*D_r*w
        result.append(sig_R)
        return result

    @staticmethod
    def eqn_10_1__w(D_r: float, sig_R: float, **kwargs):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        w = 229.357798165138*sig_R/D_r
        result.append(w)
        return result

