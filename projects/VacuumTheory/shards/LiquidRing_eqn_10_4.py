from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class LiquidRing:
    @kwasak_static
    def eqn_10_4(Q_gas=None, SP_1=None, SP_2=None, S_p=None, V=None, t=None, **kwargs):
        return

    @staticmethod
    def eqn_10_4__Q_gas(SP_1: float, SP_2: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        Q_gas = -(SP_1 - SP_2*exp(S_p*t/V))/(exp(S_p*t/V) - 1)
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_4__SP_1(Q_gas: float, SP_2: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_1 = Q_gas + (-Q_gas + SP_2)*exp(S_p*t/V)
        result.append(SP_1)
        return result

    @staticmethod
    def eqn_10_4__SP_2(Q_gas: float, SP_1: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_2 = (Q_gas*exp(S_p*t/V) - Q_gas + SP_1)*exp(-S_p*t/V)
        result.append(SP_2)
        return result

    @staticmethod
    def eqn_10_4__S_p(Q_gas: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        S_p = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_4__V(Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        V = S_p*t/log((Q_gas - SP_1)/(Q_gas - SP_2))
        result.append(V)
        return result

    @staticmethod
    def eqn_10_4__t(Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        t = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/S_p
        result.append(t)
        return result

