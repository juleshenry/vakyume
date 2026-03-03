from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class Rotary:
    @kwasak
    def eqn_11_2(self, Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None):
        return
    def eqn_11_2__Q(self, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q = (Q_0 + Q_external_gas_throughput - SP_1 + (-Q_0 + SP_2)*exp(S_vol_pump_speed*t/V))*exp(-S_vol_pump_speed*t/V)
        result.append(Q)
        return result
    def eqn_11_2__Q_0(self, Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q_0 = -(Q*exp(S_vol_pump_speed*t/V) - Q_external_gas_throughput + SP_1 - SP_2*exp(S_vol_pump_speed*t/V))/(exp(S_vol_pump_speed*t/V) - 1)
        result.append(Q_0)
        return result
    def eqn_11_2__Q_external_gas_throughput(self, Q: float, Q_0: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q_external_gas_throughput = -Q_0 + SP_1 + (Q + Q_0 - SP_2)*exp(S_vol_pump_speed*t/V)
        result.append(Q_external_gas_throughput)
        return result
    def eqn_11_2__SP_1(self, Q: float, Q_0: float, Q_external_gas_throughput: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        SP_1 = Q_0 + Q_external_gas_throughput + (-Q - Q_0 + SP_2)*exp(S_vol_pump_speed*t/V)
        result.append(SP_1)
        return result
    def eqn_11_2__SP_2(self, Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        SP_2 = (-Q_0 - Q_external_gas_throughput + SP_1 + (Q + Q_0)*exp(S_vol_pump_speed*t/V))*exp(-S_vol_pump_speed*t/V)
        result.append(SP_2)
        return result
    def eqn_11_2__S_vol_pump_speed(self, Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        S_vol_pump_speed = V*log((Q_0 + Q_external_gas_throughput - SP_1)/(Q + Q_0 - SP_2))/t
        result.append(S_vol_pump_speed)
        return result
    def eqn_11_2__V(self, Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        V = S_vol_pump_speed*t/log((Q_0 + Q_external_gas_throughput - SP_1)/(Q + Q_0 - SP_2))
        result.append(V)
        return result
    def eqn_11_2__t(self, Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, **kwargs):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        t = V*log((Q_0 + Q_external_gas_throughput - SP_1)/(Q + Q_0 - SP_2))/S_vol_pump_speed
        result.append(t)
        return result
