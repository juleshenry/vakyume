from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class VacuumTheory:
    @kwasak
    def eqn_1_10(self, P_1=None, P_2=None, T_1=None, T_2=None, V_1=None, V_2=None):
        return
    def eqn_1_10__P_1(self, P_2: float, T_1: float, T_2: float, V_1: float, V_2: float, **kwargs):
        # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return result
    def eqn_1_10__P_2(self, P_1: float, T_1: float, T_2: float, V_1: float, V_2: float, **kwargs):
        # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return result
    def eqn_1_10__T_1(self, P_1: float, P_2: float, T_2: float, V_1: float, V_2: float, **kwargs):
        # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return result
    def eqn_1_10__T_2(self, P_1: float, P_2: float, T_1: float, V_1: float, V_2: float, **kwargs):
        # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return result
    def eqn_1_10__V_1(self, P_1: float, P_2: float, T_1: float, T_2: float, V_2: float, **kwargs):
        # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return result
    def eqn_1_10__V_2(self, P_1: float, P_2: float, T_1: float, T_2: float, V_1: float, **kwargs):
        # P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return result
    @kwasak
    def eqn_1_11(self, M=None, P=None, T=None, W=None, q=None):
        """
        W := lb/hr flow
        M := molecular weight
        P := Torr
        T := R degrees temp
        """
        return
    def eqn_1_11__M(self, P: float, T: float, W: float, q: float, **kwargs):
        # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result
    def eqn_1_11__P(self, M: float, T: float, W: float, q: float, **kwargs):
        # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return result
    def eqn_1_11__T(self, M: float, P: float, W: float, q: float, **kwargs):
        # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return result
    def eqn_1_11__W(self, M: float, P: float, T: float, q: float, **kwargs):
        # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return result
    def eqn_1_11__q(self, M: float, P: float, T: float, W: float, **kwargs):
        # q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result
    @kwasak
    def eqn_1_12(self, Total_P=None, sum_partial_pressures=None):
        return
    def eqn_1_12__Total_P(self, sum_partial_pressures: float, **kwargs):
        # Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return result
    def eqn_1_12__sum_partial_pressures(self, Total_P: float, **kwargs):
        # Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return result
    @kwasak
    def eqn_1_13a(self, n=None, n_a=None, y_a=None):
        return
    def eqn_1_13a__n(self, n_a: float, y_a: float, **kwargs):
        # y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result
    def eqn_1_13a__n_a(self, n: float, y_a: float, **kwargs):
        # y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result
    def eqn_1_13a__y_a(self, n: float, n_a: float, **kwargs):
        # y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result
    @kwasak
    def eqn_1_13b(self, P=None, p_a=None, y_a=None):
        return
    def eqn_1_13b__P(self, p_a: float, y_a: float, **kwargs):
        # y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return result
    def eqn_1_13b__p_a(self, P: float, y_a: float, **kwargs):
        # y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return result
    def eqn_1_13b__y_a(self, P: float, p_a: float, **kwargs):
        # y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return result
    @kwasak
    def eqn_1_3(self, T=None, k=None, m=None, v=None):
        """
        k:= boltzmann constant
        kboltz:= 1.38e-16
        avogad:= 6.02e23
        """
        return
    def eqn_1_3__T(self, k: float, m: float, v: float, **kwargs):
        # .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return result
    def eqn_1_3__k(self, T: float, m: float, v: float, **kwargs):
        # .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return result
    def eqn_1_3__m(self, T: float, k: float, v: float, **kwargs):
        # .5 * m * v**2 = 1.5 * k * T
        result = []
        m = 3.0*T*k/v**2
        result.append(m)
        return result
    def eqn_1_3__v(self, T: float, k: float, m: float, **kwargs):
        # .5 * m * v**2 = 1.5 * k * T
        result = []
        v = -1.73205080756888*sqrt(T*k/m)
        result.append(v)
        v = 1.73205080756888*sqrt(T*k/m)
        result.append(v)
        return result
    @kwasak
    def eqn_1_7(self, R=None, T=None, V=None, n=None, p=None):
        return
    def eqn_1_7__R(self, T: float, V: float, n: float, p: float, **kwargs):
        # p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return result
    def eqn_1_7__T(self, R: float, V: float, n: float, p: float, **kwargs):
        # p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return result
    def eqn_1_7__V(self, R: float, T: float, n: float, p: float, **kwargs):
        # p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return result
    def eqn_1_7__n(self, R: float, T: float, V: float, p: float, **kwargs):
        # p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return result
    def eqn_1_7__p(self, R: float, T: float, V: float, n: float, **kwargs):
        # p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return result
    @kwasak
    def eqn_1_8(self, M=None, P=None, R=None, T=None, V=None, m=None):
        return
    def eqn_1_8__M(self, P: float, R: float, T: float, V: float, m: float, **kwargs):
        # P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return result
    def eqn_1_8__P(self, M: float, R: float, T: float, V: float, m: float, **kwargs):
        # P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return result
    def eqn_1_8__R(self, M: float, P: float, T: float, V: float, m: float, **kwargs):
        # P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return result
    def eqn_1_8__T(self, M: float, P: float, R: float, V: float, m: float, **kwargs):
        # P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return result
    def eqn_1_8__V(self, M: float, P: float, R: float, T: float, m: float, **kwargs):
        # P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return result
    def eqn_1_8__m(self, M: float, P: float, R: float, T: float, V: float, **kwargs):
        # P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return result
    @kwasak
    def eqn_1_9(self, M=None, P=None, R=None, T=None, rho=None):
        return
    def eqn_1_9__M(self, P: float, R: float, T: float, rho: float, **kwargs):
        # rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result
    def eqn_1_9__P(self, M: float, R: float, T: float, rho: float, **kwargs):
        # rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return result
    def eqn_1_9__R(self, M: float, P: float, T: float, rho: float, **kwargs):
        # rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result
    def eqn_1_9__T(self, M: float, P: float, R: float, rho: float, **kwargs):
        # rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return result
    def eqn_1_9__rho(self, M: float, P: float, R: float, T: float, **kwargs):
        # rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return result
