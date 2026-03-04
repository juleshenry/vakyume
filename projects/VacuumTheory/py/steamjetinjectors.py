from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class SteamJetInjectors:
    @kwasak
    def eqn_9_1(self, A=None, rho_s=None, v=None, w_s=None):
        """
        w_s := motive steam flow rate, lb/hr
        v:= velocity
        A:= cross sectional area, ft^2
        rhos_s := motive steam density, lb/ft^3
        """
        return

    def eqn_9_1__A(self, rho_s: float, v: float, w_s: float, **kwargs):
        # w_s = v * A * rho_s
        result = []
        A = w_s / (rho_s * v)
        result.append(A)
        return result

    def eqn_9_1__rho_s(self, A: float, v: float, w_s: float, **kwargs):
        # w_s = v * A * rho_s
        result = []
        rho_s = w_s / (A * v)
        result.append(rho_s)
        return result

    def eqn_9_1__v(self, A: float, rho_s: float, w_s: float, **kwargs):
        # w_s = v * A * rho_s
        result = []
        v = w_s / (A * rho_s)
        result.append(v)
        return result

    def eqn_9_1__w_s(self, A: float, rho_s: float, v: float, **kwargs):
        # w_s = v * A * rho_s
        result = []
        w_s = A * rho_s * v
        result.append(w_s)
        return result

    @kwasak
    def eqn_9_2(self, P_m=None, d_n=None, rho_s=None, w_s=None):
        """
        d_n := nozzle throat diameter
        P_m := motive steam pressure at point 1, psia
        rhos_s := motive steam density at point 1, lb/ft^3
        """
        return

    def eqn_9_2__P_m(self, d_n: float, rho_s: float, w_s: float, **kwargs):
        # w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        P_m = 1.334027668054e-6 * w_s**2 / (d_n**4 * rho_s)
        result.append(P_m)
        return result

    def eqn_9_2__d_n(self, P_m: float, rho_s: float, w_s: float, **kwargs):
        # w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        d_n = -0.0339853079285911 * sqrt(w_s / (P_m * rho_s) ** 0.5)
        result.append(d_n)
        d_n = 0.0339853079285911 * sqrt(w_s / (P_m * rho_s) ** 0.5)
        result.append(d_n)
        return result

    def eqn_9_2__rho_s(self, P_m: float, d_n: float, w_s: float, **kwargs):
        # w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        rho_s = 1.334027668054e-6 * w_s**2 / (P_m * d_n**4)
        result.append(rho_s)
        return result

    def eqn_9_2__w_s(self, P_m: float, d_n: float, rho_s: float, **kwargs):
        # w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        w_s = 865.8 * d_n**2 * sqrt(P_m * rho_s)
        result.append(w_s)
        return result

    @kwasak
    def eqn_9_3(self, P_s=None, V=None, t_e=None, w_j=None):
        """
        t_e := time required to evacuate system, minutes
        P_s := design suction pressure of the ejector, torr
        V := free volume of process system, ft^3
        w_j := ejector capacity, 70 deg_F basis, lb/hr
        """
        return

    def eqn_9_3__P_s(self, V: float, t_e: float, w_j: float, **kwargs):
        # t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        P_s = 33.3333333333333 * (23.0 * V - 10.0 * t_e * w_j) / V
        result.append(P_s)
        return result

    def eqn_9_3__V(self, P_s: float, t_e: float, w_j: float, **kwargs):
        # t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        V = -1000.0 * t_e * w_j / (3.0 * P_s - 2300.0)
        result.append(V)
        return result

    def eqn_9_3__t_e(self, P_s: float, V: float, w_j: float, **kwargs):
        # t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        t_e = 0.001 * V * (2300.0 - 3.0 * P_s) / w_j
        result.append(t_e)
        return result

    def eqn_9_3__w_j(self, P_s: float, V: float, t_e: float, **kwargs):
        # t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        w_j = 0.001 * V * (2300.0 - 3.0 * P_s) / t_e
        result.append(w_j)
        return result

    @kwasak
    def eqn_9_4(self, AEL=None, SC=None, r=None, w_s=None):
        """
        w_s:= motive steam requirement
        r := pounds of steam required to compress 1 lb air from ejector suction pressure P_s to discharge pressure P_d
        SC := size correction factor
        """
        return

    def eqn_9_4__AEL(self, SC: float, r: float, w_s: float, **kwargs):
        # w_s = AEL * r * SC
        result = []
        AEL = w_s / (SC * r)
        result.append(AEL)
        return result

    def eqn_9_4__SC(self, AEL: float, r: float, w_s: float, **kwargs):
        # w_s = AEL * r * SC
        result = []
        SC = w_s / (AEL * r)
        result.append(SC)
        return result

    def eqn_9_4__r(self, AEL: float, SC: float, w_s: float, **kwargs):
        # w_s = AEL * r * SC
        result = []
        r = w_s / (AEL * SC)
        result.append(r)
        return result

    def eqn_9_4__w_s(self, AEL: float, SC: float, r: float, **kwargs):
        # w_s = AEL * r * SC
        result = []
        w_s = AEL * SC * r
        result.append(w_s)
        return result

    @kwasak
    def eqn_9_5(self, V=None, r_h=None, t_h=None, w_h=None):
        """
        w_h:= motive steam hogging
        r_h:=pounds of 100-psig stream required per cubic foot
        V:= process system free volume, ft^3
        t_h := time permitted for evatuation, hr
        """
        return

    def eqn_9_5__V(self, r_h: float, t_h: float, w_h: float, **kwargs):
        # w_h = r_h * V / t_h
        result = []
        V = t_h * w_h / r_h
        result.append(V)
        return result

    def eqn_9_5__r_h(self, V: float, t_h: float, w_h: float, **kwargs):
        # w_h = r_h * V / t_h
        result = []
        r_h = t_h * w_h / V
        result.append(r_h)
        return result

    def eqn_9_5__t_h(self, V: float, r_h: float, w_h: float, **kwargs):
        # w_h = r_h * V / t_h
        result = []
        t_h = V * r_h / w_h
        result.append(t_h)
        return result

    def eqn_9_5__w_h(self, V: float, r_h: float, t_h: float, **kwargs):
        # w_h = r_h * V / t_h
        result = []
        w_h = V * r_h / t_h
        result.append(w_h)
        return result
