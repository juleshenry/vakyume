from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class PressMgmt:
    @kwasak
    def eqn_3_1(self, Abs_Pressure=None, BarometricPressure=None, Vacuum=None):
        return
    def eqn_3_1__Abs_Pressure(self, BarometricPressure: float, Vacuum: float, **kwargs):
        # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return result
    def eqn_3_1__BarometricPressure(self, Abs_Pressure: float, Vacuum: float, **kwargs):
        # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return result
    def eqn_3_1__Vacuum(self, Abs_Pressure: float, BarometricPressure: float, **kwargs):
        # Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return result
    @kwasak
    def eqn_3_11(self, A_C=None, H_2=None, P=None, V=None):
        return
    def eqn_3_11__A_C(self, H_2: float, P: float, V: float, **kwargs):
        # P = A_C / V * (H_2) ** 2
        result = []
        A_C = P*V/H_2**2
        result.append(A_C)
        return result
    def eqn_3_11__H_2(self, A_C: float, P: float, V: float, **kwargs):
        # P = A_C / V * (H_2) ** 2
        result = []
        H_2 = -sqrt(P*V/A_C)
        result.append(H_2)
        H_2 = sqrt(P*V/A_C)
        result.append(H_2)
        return result
    def eqn_3_11__P(self, A_C: float, H_2: float, V: float, **kwargs):
        # P = A_C / V * (H_2) ** 2
        result = []
        P = A_C*H_2**2/V
        result.append(P)
        return result
    def eqn_3_11__V(self, A_C: float, H_2: float, P: float, **kwargs):
        # P = A_C / V * (H_2) ** 2
        result = []
        V = A_C*H_2**2/P
        result.append(V)
        return result
    @kwasak
    def eqn_3_12(self, H_2=None, KAPPA_1=None, P=None):
        """
        KAPPA := A_C / V, THE `GAUGE CONSTANT`
        """
        return
    def eqn_3_12__H_2(self, KAPPA_1: float, P: float, **kwargs):
        # P = KAPPA_1 * H_2 ** 2
        result = []
        H_2 = -sqrt(P/KAPPA_1)
        result.append(H_2)
        H_2 = sqrt(P/KAPPA_1)
        result.append(H_2)
        return result
    def eqn_3_12__KAPPA_1(self, H_2: float, P: float, **kwargs):
        # P = KAPPA_1 * H_2 ** 2
        result = []
        KAPPA_1 = P/H_2**2
        result.append(KAPPA_1)
        return result
    def eqn_3_12__P(self, H_2: float, KAPPA_1: float, **kwargs):
        # P = KAPPA_1 * H_2 ** 2
        result = []
        P = H_2**2*KAPPA_1
        result.append(P)
        return result
    @kwasak
    def eqn_3_13(self, H_1=None, H_2=None, KAPPA_2=None, P=None):
        return
    def eqn_3_13__H_1(self, H_2: float, KAPPA_2: float, P: float, **kwargs):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_1 = H_2 - P/KAPPA_2
        result.append(H_1)
        return result
    def eqn_3_13__H_2(self, H_1: float, KAPPA_2: float, P: float, **kwargs):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_2 = H_1 + P/KAPPA_2
        result.append(H_2)
        return result
    def eqn_3_13__KAPPA_2(self, H_1: float, H_2: float, P: float, **kwargs):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        KAPPA_2 = -P/(H_1 - H_2)
        result.append(KAPPA_2)
        return result
    def eqn_3_13__P(self, H_1: float, H_2: float, KAPPA_2: float, **kwargs):
        # P = KAPPA_2 * (H_2 - H_1)
        result = []
        P = KAPPA_2*(-H_1 + H_2)
        result.append(P)
        return result
    @kwasak
    def eqn_3_15(self, V_PMIN=None):
        """
        V_PMIN := `PRACTICAL MIN, 1982`
        """
        return
    def eqn_3_15__V_PMIN(self, **kwargs):
        # V_PMIN = 3.141592653589793 / 4
        result = []
        V_PMIN = 0.785398163397448
        result.append(V_PMIN)
        return result
    @kwasak
    def eqn_3_16(self, V_div_V_P_MAX=None):
        return
    def eqn_3_16__V_div_V_P_MAX(self, **kwargs):
        # V_div_V_P_MAX = 200000 / (3.141592653589793 / 4)
        result = []
        V_div_V_P_MAX = 254647.908947033
        result.append(V_div_V_P_MAX)
        return result
    @kwasak
    def eqn_3_17(self, P_MIN=None):
        return
    def eqn_3_17__P_MIN(self, **kwargs):
        # P_MIN = (3.141592653589793 / 4) / (200000)
        result = []
        P_MIN = 0.00000392699081698724
        result.append(P_MIN)
        return result
    @kwasak
    def eqn_3_2(self, G=None, G_C=None, H=None, P=None, rho=None):
        return
    def eqn_3_2__G(self, G_C: float, H: float, P: float, rho: float, **kwargs):
        # P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return result
    def eqn_3_2__G_C(self, G: float, H: float, P: float, rho: float, **kwargs):
        # P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return result
    def eqn_3_2__H(self, G: float, G_C: float, P: float, rho: float, **kwargs):
        # P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return result
    def eqn_3_2__P(self, G: float, G_C: float, H: float, rho: float, **kwargs):
        # P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return result
    def eqn_3_2__rho(self, G: float, G_C: float, H: float, P: float, **kwargs):
        # P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return result
    @kwasak
    def eqn_3_3(self, H_1=None, H_2=None, P=None, P_P=None):
        return
    def eqn_3_3__H_1(self, H_2: float, P: float, P_P: float, **kwargs):
        # P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return result
    def eqn_3_3__H_2(self, H_1: float, P: float, P_P: float, **kwargs):
        # P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 - P + P_P
        result.append(H_2)
        return result
    def eqn_3_3__P(self, H_1: float, H_2: float, P_P: float, **kwargs):
        # P_P - P = H_2 - H_1
        result = []
        P = H_1 - H_2 + P_P
        result.append(P)
        return result
    def eqn_3_3__P_P(self, H_1: float, H_2: float, P: float, **kwargs):
        # P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return result
    @kwasak
    def eqn_3_4(self, KAPPA=None, P=None, V=None):
        return
    def eqn_3_4__KAPPA(self, P: float, V: float, **kwargs):
        # P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return result
    def eqn_3_4__P(self, KAPPA: float, V: float, **kwargs):
        # P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return result
    def eqn_3_4__V(self, KAPPA: float, P: float, **kwargs):
        # P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return result
    @kwasak
    def eqn_3_5(self, P=None, P_P=None, V=None, V_P=None):
        return
    def eqn_3_5__P(self, P_P: float, V: float, V_P: float, **kwargs):
        # P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return result
    def eqn_3_5__P_P(self, P: float, V: float, V_P: float, **kwargs):
        # P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return result
    def eqn_3_5__V(self, P: float, P_P: float, V_P: float, **kwargs):
        # P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return result
    def eqn_3_5__V_P(self, P: float, P_P: float, V: float, **kwargs):
        # P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return result
    @kwasak
    def eqn_3_6(self, H_1=None, H_2=None, P=None, V=None, V_P=None):
        return
    def eqn_3_6__H_1(self, H_2: float, P: float, V: float, V_P: float, **kwargs):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return result
    def eqn_3_6__H_2(self, H_1: float, P: float, V: float, V_P: float, **kwargs):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P*V/V_P - P
        result.append(H_2)
        return result
    def eqn_3_6__P(self, H_1: float, H_2: float, V: float, V_P: float, **kwargs):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P*(-H_1 + H_2)/(V - V_P)
        result.append(P)
        return result
    def eqn_3_6__V(self, H_1: float, H_2: float, P: float, V_P: float, **kwargs):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return result
    def eqn_3_6__V_P(self, H_1: float, H_2: float, P: float, V: float, **kwargs):
        # P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return result
    @kwasak
    def eqn_3_8(self, A_C=None, H_2=None, V_P=None):
        return
    def eqn_3_8__A_C(self, H_2: float, V_P: float, **kwargs):
        # V_P = A_C * H_2
        result = []
        A_C = V_P/H_2
        result.append(A_C)
        return result
    def eqn_3_8__H_2(self, A_C: float, V_P: float, **kwargs):
        # V_P = A_C * H_2
        result = []
        H_2 = V_P/A_C
        result.append(H_2)
        return result
    def eqn_3_8__V_P(self, A_C: float, H_2: float, **kwargs):
        # V_P = A_C * H_2
        result = []
        V_P = A_C*H_2
        result.append(V_P)
        return result
    @kwasak
    def eqn_3_9(self, A_C=None, H_1=None, H_2=None, P=None, V=None):
        return
    def eqn_3_9__A_C(self, H_1: float, H_2: float, P: float, V: float, **kwargs):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        A_C = P*V/(H_2*(-H_1 + H_2 + P))
        result.append(A_C)
        return result
    def eqn_3_9__H_1(self, A_C: float, H_2: float, P: float, V: float, **kwargs):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_1 = H_2 + P - P*V/(A_C*H_2)
        result.append(H_1)
        return result
    def eqn_3_9__H_2(self, A_C: float, H_1: float, P: float, V: float, **kwargs):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_2 = (A_C*(H_1 - P) - sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        H_2 = (A_C*(H_1 - P) + sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        return result
    def eqn_3_9__P(self, A_C: float, H_1: float, H_2: float, V: float, **kwargs):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        P = A_C*H_2*(H_1 - H_2)/(A_C*H_2 - V)
        result.append(P)
        return result
    def eqn_3_9__V(self, A_C: float, H_1: float, H_2: float, P: float, **kwargs):
        # P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        V = A_C*H_2*(-H_1 + H_2 + P)/P
        result.append(V)
        return result
