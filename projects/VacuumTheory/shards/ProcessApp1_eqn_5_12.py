from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class ProcessApp1:
    @kwasak_static
    def eqn_5_12(Eff=None, N_ES=None, N_t=None, T=None, **kwargs):
        return

    @staticmethod
    def eqn_5_12__Eff(N_ES: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result

    @staticmethod
    def eqn_5_12__N_ES(Eff: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result

    @staticmethod
    def eqn_5_12__N_t(Eff: float, N_ES: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result

    @staticmethod
    def eqn_5_12__T(Eff: float, N_ES: float, N_t: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result

