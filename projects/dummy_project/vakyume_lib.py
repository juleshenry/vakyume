from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class Basic:
    @kwasak
    def eqn_1_1(x=None, y=None, z=None, **kwargs):
        return
    def eqn_1_1__x(y: float, z: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        x = -y + z
        result.append(x)
        return result
    def eqn_1_1__y(x: float, z: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        y = -x + z
        result.append(y)
        return result
    def eqn_1_1__z(x: float, y: float, **kwargs):
        # [.pyeqn] z = x + y
        result = []
        z = x + y
        result.append(z)
        return result

class Rotary:
    @kwasak
    def eqn_11_2(Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None, **kwargs):
        return
    def eqn_11_2__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q = (Q_0 + Q_external_gas_throughput - SP_1 + (-Q_0 + SP_2)*exp(S_vol_pump_speed*t/V))*exp(-S_vol_pump_speed*t/V)
        result.append(Q)
        return result
    def eqn_11_2__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q_0 = -(Q*exp(S_vol_pump_speed*t/V) - Q_external_gas_throughput + SP_1 - SP_2*exp(S_vol_pump_speed*t/V))/(exp(S_vol_pump_speed*t/V) - 1)
        result.append(Q_0)
        return result
    def eqn_11_2__Q_external_gas_throughput(Q: float, Q_0: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q_external_gas_throughput = -Q_0 + SP_1 + (Q + Q_0 - SP_2)*exp(S_vol_pump_speed*t/V)
        result.append(Q_external_gas_throughput)
        return result
    def eqn_11_2__SP_1(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        SP_1 = Q_0 + Q_external_gas_throughput + (-Q - Q_0 + SP_2)*exp(S_vol_pump_speed*t/V)
        result.append(SP_1)
        return result
    def eqn_11_2__SP_2(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        SP_2 = (-Q_0 - Q_external_gas_throughput + SP_1 + (Q + Q_0)*exp(S_vol_pump_speed*t/V))*exp(-S_vol_pump_speed*t/V)
        result.append(SP_2)
        return result
    def eqn_11_2__S_vol_pump_speed(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        S_vol_pump_speed = V*log((Q_0 + Q_external_gas_throughput - SP_1)/(Q + Q_0 - SP_2))/t
        result.append(S_vol_pump_speed)
        return result
    def eqn_11_2__V(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        V = S_vol_pump_speed*t/log((Q_0 + Q_external_gas_throughput - SP_1)/(Q + Q_0 - SP_2))
        result.append(V)
        return result
    def eqn_11_2__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        t = V*log((Q_0 + Q_external_gas_throughput - SP_1)/(Q + Q_0 - SP_2))/S_vol_pump_speed
        result.append(t)
        return result
