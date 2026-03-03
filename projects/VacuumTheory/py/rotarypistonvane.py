from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class RotaryPistonVane:
    @kwasak
    def eqn_11_1(self, PS=None, Q_0=None, Q_external_gas_throughput=None, V=None, dP=None, dT=None):
        """
        Q_0 := throughput of gas flow due to system outgassing
        """
        return
    def eqn_11_1__PS(self, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float, **kwargs):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result
    def eqn_11_1__Q_0(self, PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float, **kwargs):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result
    def eqn_11_1__Q_external_gas_throughput(self, PS: float, Q_0: float, V: float, dP: float, dT: float, **kwargs):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result
    def eqn_11_1__V(self, PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float, **kwargs):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result
    def eqn_11_1__dP(self, PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float, **kwargs):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result
    def eqn_11_1__dT(self, PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, **kwargs):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result
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
    @kwasak
    def eqn_11_3(self, F_s=None, t=None, t_c=None):
        """
        t:= actual evacuation time
        t_c:= calculated evacuation time using Eq 10.4
        F_s:= system factor, based on operating experience
        """
        return
    def eqn_11_3__F_s(self, t: float, t_c: float, **kwargs):
        # t = t_c * F_s
        result = []
        F_s = t/t_c
        result.append(F_s)
        return result
    def eqn_11_3__t(self, F_s: float, t_c: float, **kwargs):
        # t = t_c * F_s
        result = []
        t = F_s*t_c
        result.append(t)
        return result
    def eqn_11_3__t_c(self, F_s: float, t: float, **kwargs):
        # t = t_c * F_s
        result = []
        t_c = t/F_s
        result.append(t_c)
        return result
    @kwasak
    def eqn_11_4(self, p_g=None, p_s=None, p_v=None):
        """
        p_v := partial pressure of vapor at pump suction, torr
        p_g := pressure of permanent gas at pump suction, torr
        p_s := pump suction pressure, sum of partial pressure of vapor and partial pressure of permanent gas, torr
        P_0_V := saturation pressure of vapor at pump operating temperature, torr
        P_D := pump discharge pressure, torr
        """
        return
    def eqn_11_4__p_g(self, p_s: float, p_v: float, **kwargs):
        # p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_g = p_s - p_v
        result.append(p_g)
        return result
    def eqn_11_4__p_s(self, p_g: float, p_v: float, **kwargs):
        # p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_s = p_g + p_v
        result.append(p_s)
        return result
    def eqn_11_4__p_v(self, p_g: float, p_s: float, **kwargs):
        # p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_v = 0
        result.append(p_v)
        p_v = -p_g + p_s
        result.append(p_v)
        return result
    @kwasak
    def eqn_11_5(self, P_0_v=None, P_D=None, p_g=None, p_v_max=None):
        """
        p_v_max := maximum allowable partial pressure p_v_max of the process vapor at the pump suction
        """
        return
    def eqn_11_5__P_0_v(self, P_D: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D*p_v_max/(p_g + p_v_max)
        result.append(P_0_v)
        return result
    def eqn_11_5__P_D(self, P_0_v: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v*(p_g + p_v_max)/p_v_max
        result.append(P_D)
        return result
    def eqn_11_5__p_g(self, P_0_v: float, P_D: float, p_v_max: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max*(-P_0_v + P_D)/P_0_v
        result.append(p_g)
        return result
    def eqn_11_5__p_v_max(self, P_0_v: float, P_D: float, p_g: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v*p_g/(P_0_v - P_D)
        result.append(p_v_max)
        return result
    @kwasak
    def eqn_11_6(self, P_0_V=None, P_D=None, P_v_0=None, S_B=None, S_D=None, p_b=None, p_g=None, p_v_max=None):
        """
        P_0_v := saturation vapor pressure of a condensable vapor
        S_B := maximum permissible gas ballast flow rate, ft^3/min
        S_D := free air displacement of the vacuum pump, ft^3/min
        p_b := partial pressure of vapor in the ballast gas, e.g. partial pressure of water vapor in ATM, torr
        """
        return
    def eqn_11_6__P_0_V(self, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_g - P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(P_0_V)
        return result
    def eqn_11_6__P_D(self, P_0_V: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0*S_D*(p_g + p_v_max)/(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)
        result.append(P_D)
        return result
    def eqn_11_6__P_v_0(self, P_0_V: float, P_D: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = P_D*(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)/(S_D*(p_g + p_v_max))
        result.append(P_v_0)
        return result
    def eqn_11_6__S_B(self, P_0_V: float, P_D: float, P_v_0: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = S_D*(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)/(P_D*(P_0_V - p_b))
        result.append(S_B)
        return result
    def eqn_11_6__S_D(self, P_0_V: float, P_D: float, P_v_0: float, S_B: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = P_D*S_B*(P_0_V - p_b)/(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)
        result.append(S_D)
        return result
    def eqn_11_6__p_b(self, P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (P_0_V*P_D*S_B - P_D*S_D*p_v_max + P_v_0*S_D*p_g + P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(p_b)
        return result
    def eqn_11_6__p_g(self, P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_v_max: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (-P_0_V*P_D*S_B + P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_v_max)/(P_v_0*S_D)
        result.append(p_g)
        return result
    def eqn_11_6__p_v_max(self, P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, **kwargs):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V*P_D*S_B - P_D*S_B*p_b + P_v_0*S_D*p_g)/(S_D*(P_D - P_v_0))
        result.append(p_v_max)
        return result
