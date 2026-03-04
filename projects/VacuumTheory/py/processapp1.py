from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class ProcessApp1:
    @kwasak
    def eqn_5_1(self, K_i=None, x_i=None, y_i=None):
        """
        K_i := volatility
        y_i := component concentration, vapor
        x_i := component concentration, liquid
        """
        return

    def eqn_5_1__K_i(self, x_i: float, y_i: float, **kwargs):
        # K_i = y_i / x_i
        result = []
        K_i = y_i / x_i
        result.append(K_i)
        return result

    def eqn_5_1__x_i(self, K_i: float, y_i: float, **kwargs):
        # K_i = y_i / x_i
        result = []
        x_i = y_i / K_i
        result.append(x_i)
        return result

    def eqn_5_1__y_i(self, K_i: float, x_i: float, **kwargs):
        # K_i = y_i / x_i
        result = []
        y_i = K_i * x_i
        result.append(y_i)
        return result

    @kwasak
    def eqn_5_10a(self, D=None, L_0=None, V_1=None):
        return

    def eqn_5_10a__D(self, L_0: float, V_1: float, **kwargs):
        # L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    def eqn_5_10a__L_0(self, D: float, V_1: float, **kwargs):
        # L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    def eqn_5_10a__V_1(self, D: float, L_0: float, **kwargs):
        # L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak
    def eqn_5_10b(self, L_0=None, R=None, V_1=None):
        return

    def eqn_5_10b__L_0(self, R: float, V_1: float, **kwargs):
        # L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R * V_1 / (R + 1)
        result.append(L_0)
        return result

    def eqn_5_10b__R(self, L_0: float, V_1: float, **kwargs):
        # L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0 / (L_0 - V_1)
        result.append(R)
        return result

    def eqn_5_10b__V_1(self, L_0: float, R: float, **kwargs):
        # L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0 / R
        result.append(V_1)
        return result

    @kwasak
    def eqn_5_10c(self, D=None, L_0=None, R=None):
        return

    def eqn_5_10c__D(self, L_0: float, R: float, **kwargs):
        # (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0 / R
        result.append(D)
        return result

    def eqn_5_10c__L_0(self, D: float, R: float, **kwargs):
        # (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D * R
        result.append(L_0)
        return result

    def eqn_5_10c__R(self, D: float, L_0: float, **kwargs):
        # (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0 / D
        result.append(R)
        return result

    @kwasak
    def eqn_5_11(self, B=None, L_N=None, V_0=None):
        return

    def eqn_5_11__B(self, L_N: float, V_0: float, **kwargs):
        # L_N / V_0 = (V_0 + B) / V_0
        result = []
        B = L_N - V_0
        result.append(B)
        return result

    def eqn_5_11__L_N(self, B: float, V_0: float, **kwargs):
        # L_N / V_0 = (V_0 + B) / V_0
        result = []
        L_N = B + V_0
        result.append(L_N)
        return result

    def eqn_5_11__V_0(self, B: float, L_N: float, **kwargs):
        # L_N / V_0 = (V_0 + B) / V_0
        result = []
        V_0 = -B + L_N
        result.append(V_0)
        return result

    @kwasak
    def eqn_5_12(self, Eff=None, N_ES=None, N_t=None, T=None):
        return

    def eqn_5_12__Eff(self, N_ES: float, N_t: float, T: float, **kwargs):
        # N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES / N_t) ** (1 / T)
        result.append(Eff)
        return result

    def eqn_5_12__N_ES(self, Eff: float, N_t: float, T: float, **kwargs):
        # N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T * N_t
        result.append(N_ES)
        return result

    def eqn_5_12__N_t(self, Eff: float, N_ES: float, T: float, **kwargs):
        # N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES / Eff**T
        result.append(N_t)
        return result

    def eqn_5_12__T(self, Eff: float, N_ES: float, N_t: float, **kwargs):
        # N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES / N_t) / log(Eff)
        result.append(T)
        return result

    @kwasak
    def eqn_5_13(self, HETP=None, H_p=None, N_ES=None):
        return

    def eqn_5_13__HETP(self, H_p: float, N_ES: float, **kwargs):
        # H_p = N_ES * HETP
        result = []
        HETP = H_p / N_ES
        result.append(HETP)
        return result

    def eqn_5_13__H_p(self, HETP: float, N_ES: float, **kwargs):
        # H_p = N_ES * HETP
        result = []
        H_p = HETP * N_ES
        result.append(H_p)
        return result

    def eqn_5_13__N_ES(self, HETP: float, H_p: float, **kwargs):
        # H_p = N_ES * HETP
        result = []
        N_ES = H_p / HETP
        result.append(N_ES)
        return result

    @kwasak
    def eqn_5_14(self, M=None, P_0=None, T=None, W_E=None):
        return

    def eqn_5_14__M(self, P_0: float, T: float, W_E: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261 * T * W_E**2 / P_0**2
        result.append(M)
        return result

    def eqn_5_14__P_0(self, M: float, T: float, W_E: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926 * W_E / sqrt(M / T)
        result.append(P_0)
        return result

    def eqn_5_14__T(self, M: float, P_0: float, W_E: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889 * M * P_0**2 / W_E**2
        result.append(T)
        return result

    def eqn_5_14__W_E(self, M: float, P_0: float, T: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583 * P_0 * sqrt(M / T)
        result.append(W_E)
        return result

    @kwasak
    def eqn_5_15(self, M_1=None, M_2=None, P_0_1=None, P_0_2=None, a_M_12=None):
        return

    def eqn_5_15__M_1(
        self, M_2: float, P_0_1: float, P_0_2: float, a_M_12: float, **kwargs
    ):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2 / (P_0_2 * a_M_12 / P_0_1) ** (5 / 2)
        result.append(M_1)
        M_1 = M_2 / (P_0_2 * a_M_12 / P_0_1) ** (5 / 2)
        result.append(M_1)
        return result

    def eqn_5_15__M_2(
        self, M_1: float, P_0_1: float, P_0_2: float, a_M_12: float, **kwargs
    ):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_2 = -M_1 * (P_0_2 * a_M_12 / P_0_1) ** (5 / 2)
        result.append(M_2)
        M_2 = M_1 * (P_0_2 * a_M_12 / P_0_1) ** (5 / 2)
        result.append(M_2)
        return result

    def eqn_5_15__P_0_1(
        self, M_1: float, M_2: float, P_0_2: float, a_M_12: float, **kwargs
    ):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2 * a_M_12 / (M_2 / M_1) ** (2 / 5)
        result.append(P_0_1)
        return result

    def eqn_5_15__P_0_2(
        self, M_1: float, M_2: float, P_0_1: float, a_M_12: float, **kwargs
    ):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1 * (M_2 / M_1) ** (2 / 5) / a_M_12
        result.append(P_0_2)
        return result

    def eqn_5_15__a_M_12(
        self, M_1: float, M_2: float, P_0_1: float, P_0_2: float, **kwargs
    ):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1 * (M_2 / M_1) ** (2 / 5) / P_0_2
        result.append(a_M_12)
        return result

    @kwasak
    def eqn_5_16(self, H_i=None, p_i=None, x_i=None):
        return

    def eqn_5_16__H_i(self, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * H_i
        result = []
        H_i = p_i / x_i
        result.append(H_i)
        return result

    def eqn_5_16__p_i(self, H_i: float, x_i: float, **kwargs):
        # p_i = x_i * H_i
        result = []
        p_i = H_i * x_i
        result.append(p_i)
        return result

    def eqn_5_16__x_i(self, H_i: float, p_i: float, **kwargs):
        # p_i = x_i * H_i
        result = []
        x_i = p_i / H_i
        result.append(x_i)
        return result

    @kwasak
    def eqn_5_17(self, H_2_1=None, H_2_3=None, H_2_mi=None, x_1=None, x_3=None):
        return

    def eqn_5_17__H_2_1(
        self, H_2_3: float, H_2_mi: float, x_1: float, x_3: float, **kwargs
    ):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3 * log(H_2_3) + log(H_2_mi)) / x_1)
        result.append(H_2_1)
        return result

    def eqn_5_17__H_2_3(
        self, H_2_1: float, H_2_mi: float, x_1: float, x_3: float, **kwargs
    ):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1 * log(H_2_1) + log(H_2_mi)) / x_3)
        result.append(H_2_3)
        return result

    def eqn_5_17__H_2_mi(
        self, H_2_1: float, H_2_3: float, x_1: float, x_3: float, **kwargs
    ):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1 * log(H_2_1) + x_3 * log(H_2_3))
        result.append(H_2_mi)
        return result

    def eqn_5_17__x_1(
        self, H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float, **kwargs
    ):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3 * log(H_2_3) + log(H_2_mi)) / log(H_2_1)
        result.append(x_1)
        return result

    def eqn_5_17__x_3(
        self, H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float, **kwargs
    ):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1 * log(H_2_1) + log(H_2_mi)) / log(H_2_3)
        result.append(x_3)
        return result

    @kwasak
    def eqn_5_2a(self, K_1=None, K_2=None, alpha_1_2=None):
        return

    def eqn_5_2a__K_1(self, K_2: float, alpha_1_2: float, **kwargs):
        # alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2 * alpha_1_2
        result.append(K_1)
        return result

    def eqn_5_2a__K_2(self, K_1: float, alpha_1_2: float, **kwargs):
        # alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1 / alpha_1_2
        result.append(K_2)
        return result

    def eqn_5_2a__alpha_1_2(self, K_1: float, K_2: float, **kwargs):
        # alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1 / K_2
        result.append(alpha_1_2)
        return result

    @kwasak
    def eqn_5_2b(self, K_1=None, K_2=None, x_1=None, x_2=None, y_1=None, y_2=None):
        return

    def eqn_5_2b__K_1(
        self, K_2: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs
    ):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2 * x_2 * y_1 / (x_1 * y_2)
        result.append(K_1)
        return result

    def eqn_5_2b__K_2(
        self, K_1: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs
    ):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1 * x_1 * y_2 / (x_2 * y_1)
        result.append(K_2)
        return result

    def eqn_5_2b__x_1(
        self, K_1: float, K_2: float, x_2: float, y_1: float, y_2: float, **kwargs
    ):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2 * x_2 * y_1 / (K_1 * y_2)
        result.append(x_1)
        return result

    def eqn_5_2b__x_2(
        self, K_1: float, K_2: float, x_1: float, y_1: float, y_2: float, **kwargs
    ):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1 * x_1 * y_2 / (K_2 * y_1)
        result.append(x_2)
        return result

    def eqn_5_2b__y_1(
        self, K_1: float, K_2: float, x_1: float, x_2: float, y_2: float, **kwargs
    ):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1 * x_1 * y_2 / (K_2 * x_2)
        result.append(y_1)
        return result

    def eqn_5_2b__y_2(
        self, K_1: float, K_2: float, x_1: float, x_2: float, y_1: float, **kwargs
    ):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2 * x_2 * y_1 / (K_1 * x_1)
        result.append(y_2)
        return result

    @kwasak
    def eqn_5_3(self, P_0_i=None, p_i=None, x_i=None):
        """
        p_i := component partial pressure
        x_i := liquid component mole fraction
        P_0_i := pure component vapor pressure at equilibrium temperature
        """
        return

    def eqn_5_3__P_0_i(self, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * P_0_i
        result = []
        P_0_i = p_i / x_i
        result.append(P_0_i)
        return result

    def eqn_5_3__p_i(self, P_0_i: float, x_i: float, **kwargs):
        # p_i = x_i * P_0_i
        result = []
        p_i = P_0_i * x_i
        result.append(p_i)
        return result

    def eqn_5_3__x_i(self, P_0_i: float, p_i: float, **kwargs):
        # p_i = x_i * P_0_i
        result = []
        x_i = p_i / P_0_i
        result.append(x_i)
        return result

    @kwasak
    def eqn_5_4(self, P=None, P_0_i=None, x_i=None, y_i=None):
        return

    def eqn_5_4__P(self, P_0_i: float, x_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        P = P_0_i * x_i / y_i
        result.append(P)
        return result

    def eqn_5_4__P_0_i(self, P: float, x_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        P_0_i = P * y_i / x_i
        result.append(P_0_i)
        return result

    def eqn_5_4__x_i(self, P: float, P_0_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        x_i = P * y_i / P_0_i
        result.append(x_i)
        return result

    def eqn_5_4__y_i(self, P: float, P_0_i: float, x_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i * x_i / P
        result.append(y_i)
        return result

    @kwasak
    def eqn_5_5(self, P_0_1=None, P_0_2=None, alpha_12=None):
        return

    def eqn_5_5__P_0_1(self, P_0_2: float, alpha_12: float, **kwargs):
        # alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2 * alpha_12
        result.append(P_0_1)
        return result

    def eqn_5_5__P_0_2(self, P_0_1: float, alpha_12: float, **kwargs):
        # alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1 / alpha_12
        result.append(P_0_2)
        return result

    def eqn_5_5__alpha_12(self, P_0_1: float, P_0_2: float, **kwargs):
        # alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1 / P_0_2
        result.append(alpha_12)
        return result

    @kwasak
    def eqn_5_6(self, P_0_i=None, gamma_i=None, p_i=None, x_i=None):
        return

    def eqn_5_6__P_0_i(self, gamma_i: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i / (gamma_i * x_i)
        result.append(P_0_i)
        return result

    def eqn_5_6__gamma_i(self, P_0_i: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i / (P_0_i * x_i)
        result.append(gamma_i)
        return result

    def eqn_5_6__p_i(self, P_0_i: float, gamma_i: float, x_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i * gamma_i * x_i
        result.append(p_i)
        return result

    def eqn_5_6__x_i(self, P_0_i: float, gamma_i: float, p_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i / (P_0_i * gamma_i)
        result.append(x_i)
        return result

    @kwasak
    def eqn_5_7(self, P=None, P_0_i=None, gamma_i=None, x_i=None, y_i=None):
        return

    def eqn_5_7__P(
        self, P_0_i: float, gamma_i: float, x_i: float, y_i: float, **kwargs
    ):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i * gamma_i * x_i / y_i
        result.append(P)
        return result

    def eqn_5_7__P_0_i(
        self, P: float, gamma_i: float, x_i: float, y_i: float, **kwargs
    ):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P * y_i / (gamma_i * x_i)
        result.append(P_0_i)
        return result

    def eqn_5_7__gamma_i(
        self, P: float, P_0_i: float, x_i: float, y_i: float, **kwargs
    ):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P * y_i / (P_0_i * x_i)
        result.append(gamma_i)
        return result

    def eqn_5_7__x_i(
        self, P: float, P_0_i: float, gamma_i: float, y_i: float, **kwargs
    ):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P * y_i / (P_0_i * gamma_i)
        result.append(x_i)
        return result

    def eqn_5_7__y_i(
        self, P: float, P_0_i: float, gamma_i: float, x_i: float, **kwargs
    ):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i * gamma_i * x_i / P
        result.append(y_i)
        return result

    @kwasak
    def eqn_5_8(
        self, P_0_1=None, P_0_2=None, alpha_12=None, gamma_1=None, gamma_2=None
    ):
        return

    def eqn_5_8__P_0_1(
        self, P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float, **kwargs
    ):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2 * alpha_12 * gamma_2 / gamma_1
        result.append(P_0_1)
        return result

    def eqn_5_8__P_0_2(
        self, P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float, **kwargs
    ):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1 * gamma_1 / (alpha_12 * gamma_2)
        result.append(P_0_2)
        return result

    def eqn_5_8__alpha_12(
        self, P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float, **kwargs
    ):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1 * gamma_1 / (P_0_2 * gamma_2)
        result.append(alpha_12)
        return result

    def eqn_5_8__gamma_1(
        self, P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float, **kwargs
    ):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2 * alpha_12 * gamma_2 / P_0_1
        result.append(gamma_1)
        return result

    def eqn_5_8__gamma_2(
        self, P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float, **kwargs
    ):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1 * gamma_1 / (P_0_2 * alpha_12)
        result.append(gamma_2)
        return result

    @kwasak
    def eqn_5_9(self, D=None, L_0=None, V_1=None):
        return

    def eqn_5_9__D(self, L_0: float, V_1: float, **kwargs):
        # L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    def eqn_5_9__L_0(self, D: float, V_1: float, **kwargs):
        # L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    def eqn_5_9__V_1(self, D: float, L_0: float, **kwargs):
        # L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result
