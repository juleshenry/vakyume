from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_13(HETP=None, H_p=None, N_ES=None, **kwargs):
        return

    @staticmethod
    def eqn_5_13__HETP(H_p: float, N_ES: float, **kwargs):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result

    @staticmethod
    def eqn_5_13__H_p(HETP: float, N_ES: float, **kwargs):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result

    @staticmethod
    def eqn_5_13__N_ES(HETP: float, H_p: float, **kwargs):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result

