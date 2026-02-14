from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_31(C=None, S_p=None, S_pump_speed=None, **kwargs):
        return

    @staticmethod
    def eqn_2_31__C(S_p: float, S_pump_speed: float):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        C = S_p*S_pump_speed/(S_p - S_pump_speed)
        result.append(C)
        return result

    @staticmethod
    def eqn_2_31__S_p(C: float, S_pump_speed: float):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_p = C*S_pump_speed/(C - S_pump_speed)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_2_31__S_pump_speed(C: float, S_p: float):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_pump_speed = C*S_p/(C + S_p)
        result.append(S_pump_speed)
        return result

