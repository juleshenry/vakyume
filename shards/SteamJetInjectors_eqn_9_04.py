from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SteamJetInjectors:
    @kwasak_static
    def eqn_9_04(
        AEL: float = None,
        SC: float = None,
        r: float = None,
        w_s: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_9_04__AEL(SC: float, r: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        AEL = w_s / (SC * r)
        result.append(AEL)
        return result

    @staticmethod
    def eqn_9_04__SC(AEL: float, r: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        SC = w_s / (AEL * r)
        result.append(SC)
        return result

    @staticmethod
    def eqn_9_04__r(AEL: float, SC: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        r = w_s / (AEL * SC)
        result.append(r)
        return result

    @staticmethod
    def eqn_9_04__w_s(AEL: float, SC: float, r: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        w_s = AEL * SC * r
        result.append(w_s)
        return result


