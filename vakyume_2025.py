from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np


class VacuumTheory:

    @kwasak_static
    def eqn_1_03(T: float = None, k: float = None, m: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_03__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_1_03__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_1_03__m(T: float, k: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        m = 3.0*T*k/v**2
        result.append(m)
        return result

    @staticmethod
    def eqn_1_03__v(T: float, k: float, m: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        v = -1.73205080756888*sqrt(T*k/m)
        result.append(v)
        v = 1.73205080756888*sqrt(T*k/m)
        result.append(v)
        return result

    @kwasak_static
    def eqn_1_07(R: float = None, T: float = None, V: float = None, n: float = None, p: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_07__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_07__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_07__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return result

    @staticmethod
    def eqn_1_07__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return result

    @staticmethod
    def eqn_1_07__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return result

    @kwasak_static
    def eqn_1_08(M: float = None, P: float = None, R: float = None, T: float = None, V: float = None, m: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_08__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_08__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_08__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_08__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_08__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return result

    @staticmethod
    def eqn_1_08__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return result

    @kwasak_static
    def eqn_1_09(M: float = None, P: float = None, R: float = None, T: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_09__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result

    @staticmethod
    def eqn_1_09__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return result

    @staticmethod
    def eqn_1_09__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_09__T(M: float, P: float, R: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_09__rho(M: float, P: float, R: float, T: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_1_10(P_1: float = None, P_2: float = None, T_1: float = None, T_2: float = None, V_1: float = None, V_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return result

    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return result

    @kwasak_static
    def eqn_1_11(M: float = None, P: float = None, T: float = None, W: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return result

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result

    @kwasak_static
    def eqn_1_12(Total_P: float = None, sum_partial_pressures: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return result

    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return result

    @kwasak_static
    def eqn_1_13a(n: float = None, n_a: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result

    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result

    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result

    @kwasak_static
    def eqn_1_13b(P: float = None, p_a: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return result

    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return result

    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return result


class FluidFlowVacuumLines:

    @kwasak_static
    def eqn_2_01(D: float = None, Re: float = None, mu: float = None, rho: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_01__D(Re: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_01__Re(D: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_01__mu(D: float, Re: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_01__rho(D: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_01__v(D: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_02(delta: float = None, lambd: float = None, psi: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_02__delta(lambd: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        delta = -0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        delta = 0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        return result

    @staticmethod
    def eqn_2_02__lambd(delta: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        lambd = 4.44288293815837*delta**2*psi
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_02__psi(delta: float, lambd: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        psi = 0.225079079039277*lambd/delta**2
        result.append(psi)
        return result

    @kwasak_static
    def eqn_2_03(D: float = None, kn: float = None, lambd: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_03__D(kn: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return result

    @staticmethod
    def eqn_2_03__kn(D: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return result

    @staticmethod
    def eqn_2_03__lambd(D: float, kn: float):
        # [.pyeqn] kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return result

    @kwasak_static
    def eqn_2_04(_beta: float = None, mu: float = None, vel_grad: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_04___beta(mu: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return result

    @staticmethod
    def eqn_2_04__mu(_beta: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_04__vel_grad(_beta: float, mu: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return result

    @kwasak_static
    def eqn_2_05(D: float = None, L: float = None, delta_P: float = None, mu: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_05__D(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        D = -2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = -2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_05__L(D: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703*D**4*delta_P/(mu*q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_05__delta_P(D: float, L: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252*L*mu*q/D**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_05__mu(D: float, L: float, delta_P: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703*D**4*delta_P/(L*q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_05__q(D: float, L: float, delta_P: float, mu: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703*D**4*delta_P/(L*mu)
        result.append(q)
        return result

    @kwasak_static
    def eqn_2_06(lambd: float = None, mu: float = None, rho: float = None, v_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_06__lambd(mu: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_06__mu(lambd: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_06__rho(lambd: float, mu: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_06__v_a(lambd: float, mu: float, rho: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return result

    @kwasak_static
    def eqn_2_07(T: float = None, k: float = None, m: float = None, v_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_07__T(k: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        T = 0.392699081698724*m*v_a**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_2_07__k(T: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        k = 0.392699081698724*m*v_a**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_2_07__m(T: float, k: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return result

    @staticmethod
    def eqn_2_07__v_a(T: float, k: float, m: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return result

    @kwasak_static
    def eqn_2_08(M: float = None, P_c: float = None, T_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return result

    @staticmethod
    def eqn_2_08__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result

    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

    @staticmethod
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_08(M: float = None, P_c: float = None, T_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return result

    @staticmethod
    def eqn_2_08__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result

    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

    @staticmethod
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_10(Suc_Pres: float = None, delta_P: float = None, oper_press: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_10__Suc_Pres(delta_P: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        Suc_Pres = -delta_P + oper_press
        result.append(Suc_Pres)
        return result

    @staticmethod
    def eqn_2_10__delta_P(Suc_Pres: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        delta_P = -Suc_Pres + oper_press
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return result

    @kwasak_static
    def eqn_2_11(D: float = None, L: float = None, f: float = None, g_c: float = None, h_r: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_11__g_c(D: float, L: float, f: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        g_c = L*f*v**2/(2*D*h_r)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        h_r = L*f*v**2/(2*D*g_c)
        result.append(h_r)
        return result

    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_12(L: float = None, d: float = None, delta_P: float = None, f: float = None, g: float = None, rho: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838*d*delta_P*g/(f*rho*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155*L*f*rho*v**2/(delta_P*g)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_12__delta_P(L: float, d: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        delta_P = 2.155*L*f*rho*v**2/(d*g)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_12__f(L: float, d: float, delta_P: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        f = 0.464037122969838*d*delta_P*g/(L*rho*v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155*L*f*rho*v**2/(d*delta_P)
        result.append(g)
        return result

    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838*d*delta_P*g/(L*f*v**2)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_12__v(L: float, d: float, delta_P: float, f: float, g: float, rho: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        v = -0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        v = 0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_13(L: float = None, d: float = None, delta_P: float = None, f: float = None, q: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_13__d(L: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        d = 1.16543402167043*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) - 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) + 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) - 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) + 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_13__delta_P(L: float, d: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        delta_P = 2.15*L*f*q**2*rho/d**5
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_13__f(L: float, d: float, delta_P: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        f = 0.465116279069767*d**5*delta_P/(L*q**2*rho)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_13__q(L: float, d: float, delta_P: float, f: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        return result

    @staticmethod
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_2_14(M: float = None, R: float = None, T: float = None, g_c: float = None, k: float = None, v_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return result

    @staticmethod
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return result

    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return result

    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return result

    @staticmethod
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return result

    @kwasak_static
    def eqn_2_15(Re: float = None, f: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_15__Re(f: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736/f**4
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_15__f(Re: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316/Re**(1/4)
        result.append(f)
        return result

    @kwasak_static
    def eqn_2_17(L: float = None, d: float = None, delta_P: float = None, mu: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        d = -0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        d = 0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        delta_P = 0.0345*L*mu*v/d**2
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        mu = 28.9855072463768*d**2*delta_P/(L*v)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__v(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_17(L: float = None, d: float = None, delta_P: float = None, mu: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        L = 9.52380952380952*d**4*delta_P/(mu*q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        d = -0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = -0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        delta_P = 0.105*L*mu*q/d**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        mu = 9.52380952380952*d**4*delta_P/(L*q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__q(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952*d**4*delta_P/(L*mu)
        result.append(q)
        return result


class PressMgmt:

    @kwasak_static
    def eqn_3_01(Abs_Pressure: float = None, BarometricPressure: float = None, Vacuum: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_01__Abs_Pressure(BarometricPressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return result

    @staticmethod
    def eqn_3_01__BarometricPressure(Abs_Pressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return result

    @staticmethod
    def eqn_3_01__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return result

    @kwasak_static
    def eqn_3_02(G: float = None, G_C: float = None, H: float = None, P: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_02__G(G_C: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return result

    @staticmethod
    def eqn_3_02__G_C(G: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return result

    @staticmethod
    def eqn_3_02__H(G: float, G_C: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return result

    @staticmethod
    def eqn_3_02__P(G: float, G_C: float, H: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_02__rho(G: float, G_C: float, H: float, P: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_3_03(H_1: float = None, H_2: float = None, P: float = None, P_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_03__H_1(H_2: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_03__H_2(H_1: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 - P + P_P
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_03__P(H_1: float, H_2: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P = H_1 - H_2 + P_P
        result.append(P)
        return result

    @staticmethod
    def eqn_3_03__P_P(H_1: float, H_2: float, P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return result

    @kwasak_static
    def eqn_3_04(KAPPA: float = None, P: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_04__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return result

    @staticmethod
    def eqn_3_04__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_04__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return result

    @kwasak_static
    def eqn_3_05(P: float = None, P_P: float = None, V: float = None, V_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_05__P(P_P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_05__P_P(P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return result

    @staticmethod
    def eqn_3_05__V(P: float, P_P: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_05__V_P(P: float, P_P: float, V: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return result

    @kwasak_static
    def eqn_3_06(H_1: float = None, H_2: float = None, P: float = None, V: float = None, V_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_06__H_1(H_2: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_06__H_2(H_1: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P*V/V_P - P
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_06__P(H_1: float, H_2: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P*(-H_1 + H_2)/(V - V_P)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_06__V(H_1: float, H_2: float, P: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_06__V_P(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return result


class AirLeak:

    @kwasak_static
    def eqn_4_07(W: float = None, W_T: float = None, sum_individual_leak_rates: float = None,**kwargs):
        return


    @staticmethod
    def eqn_4_07__W(W_T: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return result

    @staticmethod
    def eqn_4_07__W_T(W: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return result

    @staticmethod
    def eqn_4_07__sum_individual_leak_rates(W: float, W_T: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result

    @kwasak_static
    def eqn_4_10(T: float = None, V: float = None, del_P: float = None, leakage: float = None, t: float = None,**kwargs):
        return


    @staticmethod
    def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        T = 3.127*V*del_P/(leakage*t)
        result.append(T)
        return result

    @staticmethod
    def eqn_4_10__V(T: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        V = 0.319795330988168*T*leakage*t/del_P
        result.append(V)
        return result

    @staticmethod
    def eqn_4_10__del_P(T: float, V: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        del_P = 0.319795330988168*T*leakage*t/V
        result.append(del_P)
        return result

    @staticmethod
    def eqn_4_10__leakage(T: float, V: float, del_P: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        leakage = 3.127*V*del_P/(T*t)
        result.append(leakage)
        return result

    @staticmethod
    def eqn_4_10__t(T: float, V: float, del_P: float, leakage: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        t = 3.127*V*del_P/(T*leakage)
        result.append(t)
        return result


class ProcessAppI:

    @kwasak_static
    def eqn_5_01(K_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_01__K_i(x_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i/x_i
        result.append(K_i)
        return result

    @staticmethod
    def eqn_5_01__x_i(K_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_01__y_i(K_i: float, x_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_02a(K_1: float = None, K_2: float = None, alpha_1_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02a__K_1(K_2: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2*alpha_1_2
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02a__K_2(K_1: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1/alpha_1_2
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02a__alpha_1_2(K_1: float, K_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1/K_2
        result.append(alpha_1_2)
        return result

    @kwasak_static
    def eqn_5_02b(K_1: float = None, K_2: float = None, x_1: float = None, x_2: float = None, y_1: float = None, y_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02b__K_1(K_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02b__K_2(K_1: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02b__x_1(K_1: float, K_2: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_02b__x_2(K_1: float, K_2: float, x_1: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result

    @staticmethod
    def eqn_5_02b__y_1(K_1: float, K_2: float, x_1: float, x_2: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result

    @staticmethod
    def eqn_5_02b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
        result.append(y_2)
        return result

    @kwasak_static
    def eqn_5_03(P_0_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_03__P_0_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_03__p_i(P_0_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_03__x_i(P_0_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_04(P: float = None, P_0_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_04__P(P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_04__P_0_i(P: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P_0_i = P*y_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_04__x_i(P: float, P_0_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        x_i = P*y_i/P_0_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_04__y_i(P: float, P_0_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_05(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_05__P_0_1(P_0_2: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2*alpha_12
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_05__P_0_2(P_0_1: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1/alpha_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_05__alpha_12(P_0_1: float, P_0_2: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return result

    @kwasak_static
    def eqn_5_06(P_0_i: float = None, gamma_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_06__P_0_i(gamma_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_06__gamma_i(P_0_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_06__p_i(P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_06__x_i(P_0_i: float, gamma_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_07(P: float = None, P_0_i: float = None, gamma_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_07__P(P_0_i: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i*gamma_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_07__P_0_i(P: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P*y_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_07__gamma_i(P: float, P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P*y_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_07__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_07__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_08(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None, gamma_1: float = None, gamma_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_08__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_08__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_08__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result

    @staticmethod
    def eqn_5_08__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return result

    @staticmethod
    def eqn_5_08__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return result

    @kwasak_static
    def eqn_5_09(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_09__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_09__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_09__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10a(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10a__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10a__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10a__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10b(L_0: float = None, R: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10b__L_0(R: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10b__R(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return result

    @staticmethod
    def eqn_5_10b__V_1(L_0: float, R: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10c(D: float = None, L_0: float = None, R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10c__D(L_0: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10c__L_0(D: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10c__R(D: float, L_0: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
        result.append(R)
        return result

    @kwasak_static
    def eqn_5_11(B: float = None, L_N: float = None, V_0: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_11__B(L_N: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        B = L_N - V_0
        result.append(B)
        return result

    @staticmethod
    def eqn_5_11__L_N(B: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        L_N = B + V_0
        result.append(L_N)
        return result

    @staticmethod
    def eqn_5_11__V_0(B: float, L_N: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        V_0 = -B + L_N
        result.append(V_0)
        return result

    @kwasak_static
    def eqn_5_12(Eff: float = None, N_ES: float = None, N_t: float = None, T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_12__Eff(N_ES: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result

    @staticmethod
    def eqn_5_12__N_ES(Eff: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result

    @staticmethod
    def eqn_5_12__N_t(Eff: float, N_ES: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result

    @staticmethod
    def eqn_5_12__T(Eff: float, N_ES: float, N_t: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result

    @kwasak_static
    def eqn_5_13(HETP: float = None, H_p: float = None, N_ES: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_13__HETP(H_p: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result

    @staticmethod
    def eqn_5_13__H_p(HETP: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result

    @staticmethod
    def eqn_5_13__N_ES(HETP: float, H_p: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result

    @kwasak_static
    def eqn_5_14(M: float = None, P_0: float = None, T: float = None, W_E: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_14__M(P_0: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return result

    @staticmethod
    def eqn_5_14__P_0(M: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return result

    @staticmethod
    def eqn_5_14__T(M: float, P_0: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return result

    @staticmethod
    def eqn_5_14__W_E(M: float, P_0: float, T: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return result

    @kwasak_static
    def eqn_5_15(M_1: float = None, M_2: float = None, P_0_1: float = None, P_0_2: float = None, a_M_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_15__M_1(M_2: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        M_1 = M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        return result

    @staticmethod
    def eqn_5_15__M_2(M_1: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_2 = -M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        M_2 = M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        return result

    @staticmethod
    def eqn_5_15__P_0_1(M_1: float, M_2: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_15__P_0_2(M_1: float, M_2: float, P_0_1: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1*(M_2/M_1)**(2/5)/a_M_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_15__a_M_12(M_1: float, M_2: float, P_0_1: float, P_0_2: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
        result.append(a_M_12)
        return result

    @kwasak_static
    def eqn_5_16(H_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_16__H_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result

    @staticmethod
    def eqn_5_16__p_i(H_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_16__x_i(H_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_17(H_2_1: float = None, H_2_3: float = None, H_2_mi: float = None, x_1: float = None, x_3: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_17__H_2_1(H_2_3: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3*log(H_2_3) + log(H_2_mi))/x_1)
        result.append(H_2_1)
        return result

    @staticmethod
    def eqn_5_17__H_2_3(H_2_1: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1*log(H_2_1) + log(H_2_mi))/x_3)
        result.append(H_2_3)
        return result

    @staticmethod
    def eqn_5_17__H_2_mi(H_2_1: float, H_2_3: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1*log(H_2_1) + x_3*log(H_2_3))
        result.append(H_2_mi)
        return result

    @staticmethod
    def eqn_5_17__x_1(H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_17__x_3(H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return result


class ProcessAppIi:

    @kwasak_static
    def eqn_6_01(T_1: float = None, T_2: float = None, T_R: float = None, c_p: float = None, del_h_v: float = None, w_1: float = None, w_2: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_01__T_1(T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R) + del_h_v*w_v)/(c_p*w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_01__T_2(T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_2 = (T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R) + del_h_v*w_v)/(c_p*w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_01__T_R(T_1: float, T_2: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_R = (T_1*c_p*w_1 + T_2*c_p*w_2 - del_h_v*w_v)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_01__c_p(T_1: float, T_2: float, T_R: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = del_h_v*w_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_01__del_h_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/w_v
        result.append(del_h_v)
        return result

    @staticmethod
    def eqn_6_01__w_1(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (-T_2*c_p*w_2 + T_R*c_p*w_2 + del_h_v*w_v)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_01__w_2(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_2 = (-T_1*c_p*w_1 + T_R*c_p*w_1 + del_h_v*w_v)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result

    @staticmethod
    def eqn_6_01__w_v(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/del_h_v
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_02(Q_v: float = None, T_1: float = None, T_2: float = None, T_R: float = None, c_p: float = None, w_1: float = None, w_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_02__Q_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_02__T_1(Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000*Q_v + T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R))/(c_p*w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_02__T_2(Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_2 = (12000*Q_v + T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R))/(c_p*w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_02__T_R(Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_R = (-12000*Q_v + T_1*c_p*w_1 + T_2*c_p*w_2)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_02__c_p(Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        c_p = 12000*Q_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_02__w_1(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_1 = (12000*Q_v - T_2*c_p*w_2 + T_R*c_p*w_2)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_02__w_2(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000*Q_v - T_1*c_p*w_1 + T_R*c_p*w_1)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result

    @kwasak_static
    def eqn_6_04(Q_v: float = None, delta_h_v: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_04__Q_v(delta_h_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        Q_v = delta_h_v*w_v/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_04__delta_h_v(Q_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        delta_h_v = 12000*Q_v/w_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_04__w_v(Q_v: float, delta_h_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_06(Q_r: float = None, delta_h_v: float = None, f_m: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_06__Q_r(delta_h_v: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        Q_r = delta_h_v*f_m/24
        result.append(Q_r)
        return result

    @staticmethod
    def eqn_6_06__delta_h_v(Q_r: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        delta_h_v = 24*Q_r/f_m
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_06__f_m(Q_r: float, delta_h_v: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        f_m = 24*Q_r/delta_h_v
        result.append(f_m)
        return result

    @kwasak_static
    def eqn_6_07(C_1: float = None, C_2: float = None, T_1: float = None, T_2: float = None, c_p: float = None, delta_h_c: float = None, delta_h_v: float = None, m_b: float = None, m_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_07__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result

    @staticmethod
    def eqn_6_07__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result

    @staticmethod
    def eqn_6_07__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*m_v)/(c_p*m_b)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_07__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*m_v)/(c_p*m_b)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_07__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*m_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_07__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*m_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result

    @staticmethod
    def eqn_6_07__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/m_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_07__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_b = delta_h_v*m_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_07__m_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/delta_h_v
        result.append(m_v)
        return result

    @kwasak_static
    def eqn_6_08(C_1: float = None, C_2: float = None, T_1: float = None, T_2: float = None, c_p: float = None, delta_h_c: float = None, delta_h_v: float = None, delta_t: float = None, m_b: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_08__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result

    @staticmethod
    def eqn_6_08__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result

    @staticmethod
    def eqn_6_08__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_08__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_08__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*delta_t*w_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_08__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*delta_t*w_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result

    @staticmethod
    def eqn_6_08__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_t*w_v)
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_08__delta_t(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_t = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*w_v)
        result.append(delta_t)
        return result

    @staticmethod
    def eqn_6_08__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        m_b = delta_h_v*delta_t*w_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_08__w_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        w_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*delta_t)
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_09(A: float = None, dV_dt: float = None, delta_P: float = None, m: float = None, mu: float = None, r: float = None, r_M: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_09__A(dV_dt: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        A = (dV_dt*r_M - sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        A = (dV_dt*r_M + sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        return result

    @staticmethod
    def eqn_6_09__dV_dt(A: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_09__delta_P(A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        delta_P = A*dV_dt*r_M/(A**2 - dV_dt*m*mu*r)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_09__m(A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
        result.append(m)
        return result

    @staticmethod
    def eqn_6_09__mu(A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_09__r(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
        result.append(r)
        return result

    @staticmethod
    def eqn_6_09__r_M(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
        result.append(r_M)
        return result

    @kwasak_static
    def eqn_6_10(A: float = None, dV_dt: float = None, delta_P: float = None, mu: float = None, r_c: float = None, s: float = None, tau: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_10__A(dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt*delta_P**(s - 1)*mu*r_c*tau
        result.append(A)
        return result

    @staticmethod
    def eqn_6_10__dV_dt(A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A*delta_P**(1 - s)/(mu*r_c*tau)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_10__delta_P(A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_10__mu(A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A*delta_P**(1 - s)/(dV_dt*r_c*tau)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_10__r_c(A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
        result.append(r_c)
        return result

    @staticmethod
    def eqn_6_10__s(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A*delta_P/(dV_dt*mu*r_c*tau))/log(delta_P)
        result.append(s)
        return result

    @staticmethod
    def eqn_6_10__tau(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A*delta_P**(1 - s)/(dV_dt*mu*r_c)
        result.append(tau)
        return result

    @kwasak_static
    def eqn_6_11a(A_d: float = None, delta_T: float = None, delta_h_i: float = None, delta_m: float = None, h_d: float = None, m_b: float = None, t_R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_11a__A_d(delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
        result.append(A_d)
        return result

    @staticmethod
    def eqn_6_11a__delta_T(A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_6_11a__delta_h_i(A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d*delta_T*h_d*t_R/(delta_m*m_b)
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_6_11a__delta_m(A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d*delta_T*h_d*t_R/(delta_h_i*m_b)
        result.append(delta_m)
        return result

    @staticmethod
    def eqn_6_11a__h_d(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i*delta_m*m_b/(A_d*delta_T*t_R)
        result.append(h_d)
        return result

    @staticmethod
    def eqn_6_11a__m_b(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_11a__t_R(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        t_R = delta_h_i*delta_m*m_b/(A_d*delta_T*h_d)
        result.append(t_R)
        return result


class Precondensors:

    @kwasak_static
    def eqn_7_01(P: float = None, p_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_01__P(p_i: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_01__p_i(P: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_01__y_i(P: float, p_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_7_02(P_i_0: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_02__P_i_0(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i/x_i
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_02__p_i(P_i_0: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        p_i = P_i_0*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_02__x_i(P_i_0: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        x_i = p_i/P_i_0
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_03(P_i_0: float = None, epsilon_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_03__P_i_0(epsilon_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i/(epsilon_i*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_03__epsilon_i(P_i_0: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i/(P_i_0*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_03__p_i(P_i_0: float, epsilon_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0*epsilon_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_03__x_i(P_i_0: float, epsilon_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i/(P_i_0*epsilon_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_04a(P: float = None, p_c: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04a__P(p_c: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return result

    @staticmethod
    def eqn_7_04a__p_c(P: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_04a__p_nc(P: float, p_c: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04aa(n_i: float = None, n_nc: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04aa__n_i(n_nc: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_i = n_nc*p_i/p_nc
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04aa__n_nc(n_i: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_nc = n_i*p_nc/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04aa__p_i(n_i: float, n_nc: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_i = n_i*p_nc/n_nc
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04aa__p_nc(n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_nc = n_nc*p_i/n_i
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04ab(P_c: float = None, p: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ab__P_c(p: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ab__p(P_c: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ab__p_i(P_c: float, p: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = 0
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04ab__p_nc(P_c: float, p: float, p_i: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = -P_c + p
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04ac(P_c: float = None, n_i: float = None, n_nc: float = None, p: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ac__P_c(n_i: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ac__n_i(P_c: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04ac__n_nc(P_c: float, n_i: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04ac__p(P_c: float, n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ac__p_i(P_c: float, n_i: float, n_nc: float, p: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result

    @kwasak_static
    def eqn_7_05(N_i: float = None, N_nc: float = None, P: float = None, P_c: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_05__N_i(N_nc: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
        return result

    @staticmethod
    def eqn_7_05__N_nc(N_i: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return result

    @staticmethod
    def eqn_7_05__P(N_i: float, N_nc: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc*p_i/N_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_05__P_c(N_i: float, N_nc: float, P: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc*p_i/N_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_05__p_i(N_i: float, N_nc: float, P: float, P_c: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i*(P - P_c)/N_nc
        result.append(p_i)
        return result

    @kwasak_static
    def eqn_7_06(M: float = None, P: float = None, P_i_0: float = None, W_air: float = None, W_i: float = None, p_c: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_06__M(P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_06__P(M: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_06__P_i_0(M: float, P: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_06__W_air(M: float, P: float, P_i_0: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*x_i)
        result.append(W_air)
        return result

    @staticmethod
    def eqn_7_06__W_i(M: float, P: float, P_i_0: float, W_air: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*x_i/(29*(P - p_c))
        result.append(W_i)
        return result

    @staticmethod
    def eqn_7_06__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_06__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_07(M: float = None, P: float = None, P_i_0: float = None, W_air: float = None, W_i: float = None, epsilon_i: float = None, p_c: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_07__M(P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_07__P(M: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_07__P_i_0(M: float, P: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*epsilon_i*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_07__W_air(M: float, P: float, P_i_0: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*epsilon_i*x_i)
        result.append(W_air)
        return result

    @staticmethod
    def eqn_7_07__W_i(M: float, P: float, P_i_0: float, W_air: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*epsilon_i*x_i/(29*(P - p_c))
        result.append(W_i)
        return result

    @staticmethod
    def eqn_7_07__epsilon_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_07__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_07__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*epsilon_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_08(L_c: float = None, Q_condensor_heat_duty: float = None, c_p: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_08__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_08__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_08__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_08__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_09(L_c: float = None, Q_condensor_heat_duty: float = None, c_p: float = None, del_T: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_09__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_09__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_09__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_09__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return result

    @staticmethod
    def eqn_7_09__rho(L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_7_10(L_c_P: float = None, Q_condensor_heat_duty: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_10__L_c_P(Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty/(500*del_T)
        result.append(L_c_P)
        return result

    @staticmethod
    def eqn_7_10__Q_condensor_heat_duty(L_c_P: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500*L_c_P*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_10__del_T(L_c_P: float, Q_condensor_heat_duty: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(500*L_c_P)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_11(Q_condensor_heat_duty: float = None, U_v: float = None, V_c: float = None, del_T_LM: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_11__Q_condensor_heat_duty(U_v: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v*V_c*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_11__U_v(Q_condensor_heat_duty: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
        result.append(U_v)
        return result

    @staticmethod
    def eqn_7_11__V_c(Q_condensor_heat_duty: float, U_v: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty/(U_v*del_T_LM)
        result.append(V_c)
        return result

    @staticmethod
    def eqn_7_11__del_T_LM(Q_condensor_heat_duty: float, U_v: float, V_c: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(U_v*V_c)
        result.append(del_T_LM)
        return result

    @kwasak_static
    def eqn_7_12(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_12__A(Q_condensor_heat_duty: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty/(U*del_T)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_12__Q_condensor_heat_duty(A: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A*U*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty/(A*del_T)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty/(A*U)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_14a(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T_LM: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14a__A(Q_condensor_heat_duty: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty/(U*del_T_LM)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14a__Q_condensor_heat_duty(A: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A*U*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14a__U(A: float, Q_condensor_heat_duty: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty/(A*del_T_LM)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14a__del_T_LM(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(A*U)
        result.append(del_T_LM)
        return result

    @kwasak_static
    def eqn_7_14b(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T_1: float = None, del_T_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14b__A(Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14b__Q_condensor_heat_duty(A: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A*U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14b__U(A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty/(A*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14b__del_T_1(A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_1)
        return result

    @staticmethod
    def eqn_7_14b__del_T_2(A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_2)
        return result

    @kwasak_static
    def eqn_7_15(U: float = None, sum_R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_15__U(sum_R: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        U = 1/sum_R
        result.append(U)
        return result

    @staticmethod
    def eqn_7_15__sum_R(U: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        sum_R = 1/U
        result.append(sum_R)
        return result

    @kwasak_static
    def eqn_7_16(D_0: float = None, D_LM: float = None, D_i: float = None, R_f_0: float = None, R_fi: float = None, U_0: float = None, h_0: float = None, h_i: float = None, k_w: float = None, x_w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_16__D_0(D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_f_0*U_0*h_0 - U_0 + h_0)/(U_0*h_0*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_16__D_LM(D_0: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_0*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(D_LM)
        return result

    @staticmethod
    def eqn_7_16__D_i(D_0: float, D_LM: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_0*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_0*x_w + D_LM*R_f_0*U_0*h_0*k_w + D_LM*U_0*k_w - D_LM*h_0*k_w))
        result.append(D_i)
        return result

    @staticmethod
    def eqn_7_16__R_f_0(D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
        result.append(R_f_0)
        return result

    @staticmethod
    def eqn_7_16__R_fi(D_0: float, D_LM: float, D_i: float, R_f_0: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_f_0/D_0 - D_i/(D_0*h_0) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result

    @staticmethod
    def eqn_7_16__U_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_0*h_i*k_w/(D_0*D_LM*R_fi*h_0*h_i*k_w + D_0*D_LM*h_0*k_w + D_0*D_i*h_0*h_i*x_w + D_LM*D_i*R_f_0*h_0*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result

    @staticmethod
    def eqn_7_16__h_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_0 = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_f_0*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_0)
        return result

    @staticmethod
    def eqn_7_16__h_i(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_0*k_w/(D_0*D_LM*R_fi*U_0*h_0*k_w + D_0*D_i*U_0*h_0*x_w + D_LM*D_i*R_f_0*U_0*h_0*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_0*k_w)
        result.append(h_i)
        return result

    @staticmethod
    def eqn_7_16__k_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_0*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(k_w)
        return result

    @staticmethod
    def eqn_7_16__x_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_f_0*k_w/D_0 - D_LM*k_w/(D_0*h_0) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result

    @kwasak_static
    def eqn_7_17(R_0: float = None, R_nc: float = None, h_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_17__R_0(R_nc: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1/h_c
        result.append(R_0)
        return result

    @staticmethod
    def eqn_7_17__R_nc(R_0: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1/h_c
        result.append(R_nc)
        return result

    @staticmethod
    def eqn_7_17__h_c(R_0: float, R_nc: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1/(R_0 - R_nc)
        result.append(h_c)
        return result

    @kwasak_static
    def eqn_7_18(D_0: float = None, D_LM: float = None, D_i: float = None, R_fi: float = None, R_fo: float = None, R_nc: float = None, U_0: float = None, h_c: float = None, h_i: float = None, k_w: float = None, x_w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_18__D_0(D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_fo*U_0*h_c - R_nc*U_0*h_c - U_0 + h_c)/(U_0*h_c*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_18__D_LM(D_0: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_c*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(D_LM)
        return result

    @staticmethod
    def eqn_7_18__D_i(D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        # [Sympy Failover]
#     def eqn_7_18__D_i(D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i) - (1 / U_0)D_i
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np


class VacuumTheory:

    @kwasak_static
    def eqn_1_03(T: float = None, k: float = None, m: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_03__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_1_03__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_1_03__m(T: float, k: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        m = 3.0*T*k/v**2
        result.append(m)
        return result

    @staticmethod
    def eqn_1_03__v(T: float, k: float, m: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        v = -1.73205080756888*sqrt(T*k/m)
        result.append(v)
        v = 1.73205080756888*sqrt(T*k/m)
        result.append(v)
        return result

    @kwasak_static
    def eqn_1_07(R: float = None, T: float = None, V: float = None, n: float = None, p: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_07__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_07__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_07__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return result

    @staticmethod
    def eqn_1_07__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return result

    @staticmethod
    def eqn_1_07__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return result

    @kwasak_static
    def eqn_1_08(M: float = None, P: float = None, R: float = None, T: float = None, V: float = None, m: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_08__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_08__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_08__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_08__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_08__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return result

    @staticmethod
    def eqn_1_08__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return result

    @kwasak_static
    def eqn_1_09(M: float = None, P: float = None, R: float = None, T: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_09__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result

    @staticmethod
    def eqn_1_09__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return result

    @staticmethod
    def eqn_1_09__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_09__T(M: float, P: float, R: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_09__rho(M: float, P: float, R: float, T: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_1_10(P_1: float = None, P_2: float = None, T_1: float = None, T_2: float = None, V_1: float = None, V_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return result

    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return result

    @kwasak_static
    def eqn_1_11(M: float = None, P: float = None, T: float = None, W: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return result

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result

    @kwasak_static
    def eqn_1_12(Total_P: float = None, sum_partial_pressures: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return result

    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return result

    @kwasak_static
    def eqn_1_13a(n: float = None, n_a: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result

    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result

    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result

    @kwasak_static
    def eqn_1_13b(P: float = None, p_a: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return result

    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return result

    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return result


class FluidFlowVacuumLines:

    @kwasak_static
    def eqn_2_01(D: float = None, Re: float = None, mu: float = None, rho: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_01__D(Re: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_01__Re(D: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_01__mu(D: float, Re: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_01__rho(D: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_01__v(D: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_02(delta: float = None, lambd: float = None, psi: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_02__delta(lambd: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        delta = -0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        delta = 0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        return result

    @staticmethod
    def eqn_2_02__lambd(delta: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        lambd = 4.44288293815837*delta**2*psi
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_02__psi(delta: float, lambd: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        psi = 0.225079079039277*lambd/delta**2
        result.append(psi)
        return result

    @kwasak_static
    def eqn_2_03(D: float = None, kn: float = None, lambd: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_03__D(kn: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return result

    @staticmethod
    def eqn_2_03__kn(D: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return result

    @staticmethod
    def eqn_2_03__lambd(D: float, kn: float):
        # [.pyeqn] kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return result

    @kwasak_static
    def eqn_2_04(_beta: float = None, mu: float = None, vel_grad: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_04___beta(mu: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return result

    @staticmethod
    def eqn_2_04__mu(_beta: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_04__vel_grad(_beta: float, mu: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return result

    @kwasak_static
    def eqn_2_05(D: float = None, L: float = None, delta_P: float = None, mu: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_05__D(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        D = -2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = -2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_05__L(D: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703*D**4*delta_P/(mu*q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_05__delta_P(D: float, L: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252*L*mu*q/D**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_05__mu(D: float, L: float, delta_P: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703*D**4*delta_P/(L*q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_05__q(D: float, L: float, delta_P: float, mu: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703*D**4*delta_P/(L*mu)
        result.append(q)
        return result

    @kwasak_static
    def eqn_2_06(lambd: float = None, mu: float = None, rho: float = None, v_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_06__lambd(mu: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_06__mu(lambd: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_06__rho(lambd: float, mu: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_06__v_a(lambd: float, mu: float, rho: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return result

    @kwasak_static
    def eqn_2_07(T: float = None, k: float = None, m: float = None, v_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_07__T(k: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        T = 0.392699081698724*m*v_a**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_2_07__k(T: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        k = 0.392699081698724*m*v_a**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_2_07__m(T: float, k: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return result

    @staticmethod
    def eqn_2_07__v_a(T: float, k: float, m: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return result

    @kwasak_static
    def eqn_2_08(M: float = None, P_c: float = None, T_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return result

    @staticmethod
    def eqn_2_08__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result

    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

    @staticmethod
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_08(M: float = None, P_c: float = None, T_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return result

    @staticmethod
    def eqn_2_08__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result

    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

    @staticmethod
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_10(Suc_Pres: float = None, delta_P: float = None, oper_press: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_10__Suc_Pres(delta_P: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        Suc_Pres = -delta_P + oper_press
        result.append(Suc_Pres)
        return result

    @staticmethod
    def eqn_2_10__delta_P(Suc_Pres: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        delta_P = -Suc_Pres + oper_press
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return result

    @kwasak_static
    def eqn_2_11(D: float = None, L: float = None, f: float = None, g_c: float = None, h_r: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_11__g_c(D: float, L: float, f: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        g_c = L*f*v**2/(2*D*h_r)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        h_r = L*f*v**2/(2*D*g_c)
        result.append(h_r)
        return result

    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_12(L: float = None, d: float = None, delta_P: float = None, f: float = None, g: float = None, rho: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838*d*delta_P*g/(f*rho*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155*L*f*rho*v**2/(delta_P*g)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_12__delta_P(L: float, d: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        delta_P = 2.155*L*f*rho*v**2/(d*g)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_12__f(L: float, d: float, delta_P: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        f = 0.464037122969838*d*delta_P*g/(L*rho*v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155*L*f*rho*v**2/(d*delta_P)
        result.append(g)
        return result

    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838*d*delta_P*g/(L*f*v**2)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_12__v(L: float, d: float, delta_P: float, f: float, g: float, rho: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        v = -0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        v = 0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_13(L: float = None, d: float = None, delta_P: float = None, f: float = None, q: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_13__d(L: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        d = 1.16543402167043*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) - 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) + 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) - 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) + 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_13__delta_P(L: float, d: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        delta_P = 2.15*L*f*q**2*rho/d**5
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_13__f(L: float, d: float, delta_P: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        f = 0.465116279069767*d**5*delta_P/(L*q**2*rho)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_13__q(L: float, d: float, delta_P: float, f: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        return result

    @staticmethod
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_2_14(M: float = None, R: float = None, T: float = None, g_c: float = None, k: float = None, v_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return result

    @staticmethod
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return result

    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return result

    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return result

    @staticmethod
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return result

    @kwasak_static
    def eqn_2_15(Re: float = None, f: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_15__Re(f: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736/f**4
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_15__f(Re: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316/Re**(1/4)
        result.append(f)
        return result

    @kwasak_static
    def eqn_2_17(L: float = None, d: float = None, delta_P: float = None, mu: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        d = -0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        d = 0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        delta_P = 0.0345*L*mu*v/d**2
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        mu = 28.9855072463768*d**2*delta_P/(L*v)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__v(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_17(L: float = None, d: float = None, delta_P: float = None, mu: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        L = 9.52380952380952*d**4*delta_P/(mu*q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        d = -0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = -0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        delta_P = 0.105*L*mu*q/d**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        mu = 9.52380952380952*d**4*delta_P/(L*q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__q(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952*d**4*delta_P/(L*mu)
        result.append(q)
        return result


class PressMgmt:

    @kwasak_static
    def eqn_3_01(Abs_Pressure: float = None, BarometricPressure: float = None, Vacuum: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_01__Abs_Pressure(BarometricPressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return result

    @staticmethod
    def eqn_3_01__BarometricPressure(Abs_Pressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return result

    @staticmethod
    def eqn_3_01__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return result

    @kwasak_static
    def eqn_3_02(G: float = None, G_C: float = None, H: float = None, P: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_02__G(G_C: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return result

    @staticmethod
    def eqn_3_02__G_C(G: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return result

    @staticmethod
    def eqn_3_02__H(G: float, G_C: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return result

    @staticmethod
    def eqn_3_02__P(G: float, G_C: float, H: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_02__rho(G: float, G_C: float, H: float, P: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_3_03(H_1: float = None, H_2: float = None, P: float = None, P_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_03__H_1(H_2: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_03__H_2(H_1: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 - P + P_P
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_03__P(H_1: float, H_2: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P = H_1 - H_2 + P_P
        result.append(P)
        return result

    @staticmethod
    def eqn_3_03__P_P(H_1: float, H_2: float, P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return result

    @kwasak_static
    def eqn_3_04(KAPPA: float = None, P: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_04__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return result

    @staticmethod
    def eqn_3_04__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_04__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return result

    @kwasak_static
    def eqn_3_05(P: float = None, P_P: float = None, V: float = None, V_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_05__P(P_P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_05__P_P(P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return result

    @staticmethod
    def eqn_3_05__V(P: float, P_P: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_05__V_P(P: float, P_P: float, V: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return result

    @kwasak_static
    def eqn_3_06(H_1: float = None, H_2: float = None, P: float = None, V: float = None, V_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_06__H_1(H_2: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_06__H_2(H_1: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P*V/V_P - P
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_06__P(H_1: float, H_2: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P*(-H_1 + H_2)/(V - V_P)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_06__V(H_1: float, H_2: float, P: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_06__V_P(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return result


class AirLeak:

    @kwasak_static
    def eqn_4_07(W: float = None, W_T: float = None, sum_individual_leak_rates: float = None,**kwargs):
        return


    @staticmethod
    def eqn_4_07__W(W_T: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return result

    @staticmethod
    def eqn_4_07__W_T(W: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return result

    @staticmethod
    def eqn_4_07__sum_individual_leak_rates(W: float, W_T: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result

    @kwasak_static
    def eqn_4_10(T: float = None, V: float = None, del_P: float = None, leakage: float = None, t: float = None,**kwargs):
        return


    @staticmethod
    def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        T = 3.127*V*del_P/(leakage*t)
        result.append(T)
        return result

    @staticmethod
    def eqn_4_10__V(T: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        V = 0.319795330988168*T*leakage*t/del_P
        result.append(V)
        return result

    @staticmethod
    def eqn_4_10__del_P(T: float, V: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        del_P = 0.319795330988168*T*leakage*t/V
        result.append(del_P)
        return result

    @staticmethod
    def eqn_4_10__leakage(T: float, V: float, del_P: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        leakage = 3.127*V*del_P/(T*t)
        result.append(leakage)
        return result

    @staticmethod
    def eqn_4_10__t(T: float, V: float, del_P: float, leakage: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        t = 3.127*V*del_P/(T*leakage)
        result.append(t)
        return result


class ProcessAppI:

    @kwasak_static
    def eqn_5_01(K_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_01__K_i(x_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i/x_i
        result.append(K_i)
        return result

    @staticmethod
    def eqn_5_01__x_i(K_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_01__y_i(K_i: float, x_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_02a(K_1: float = None, K_2: float = None, alpha_1_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02a__K_1(K_2: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2*alpha_1_2
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02a__K_2(K_1: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1/alpha_1_2
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02a__alpha_1_2(K_1: float, K_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1/K_2
        result.append(alpha_1_2)
        return result

    @kwasak_static
    def eqn_5_02b(K_1: float = None, K_2: float = None, x_1: float = None, x_2: float = None, y_1: float = None, y_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02b__K_1(K_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02b__K_2(K_1: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02b__x_1(K_1: float, K_2: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_02b__x_2(K_1: float, K_2: float, x_1: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result

    @staticmethod
    def eqn_5_02b__y_1(K_1: float, K_2: float, x_1: float, x_2: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result

    @staticmethod
    def eqn_5_02b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
        result.append(y_2)
        return result

    @kwasak_static
    def eqn_5_03(P_0_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_03__P_0_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_03__p_i(P_0_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_03__x_i(P_0_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_04(P: float = None, P_0_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_04__P(P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_04__P_0_i(P: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P_0_i = P*y_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_04__x_i(P: float, P_0_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        x_i = P*y_i/P_0_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_04__y_i(P: float, P_0_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_05(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_05__P_0_1(P_0_2: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2*alpha_12
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_05__P_0_2(P_0_1: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1/alpha_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_05__alpha_12(P_0_1: float, P_0_2: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return result

    @kwasak_static
    def eqn_5_06(P_0_i: float = None, gamma_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_06__P_0_i(gamma_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_06__gamma_i(P_0_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_06__p_i(P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_06__x_i(P_0_i: float, gamma_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_07(P: float = None, P_0_i: float = None, gamma_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_07__P(P_0_i: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i*gamma_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_07__P_0_i(P: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P*y_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_07__gamma_i(P: float, P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P*y_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_07__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_07__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_08(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None, gamma_1: float = None, gamma_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_08__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_08__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_08__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result

    @staticmethod
    def eqn_5_08__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return result

    @staticmethod
    def eqn_5_08__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return result

    @kwasak_static
    def eqn_5_09(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_09__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_09__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_09__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10a(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10a__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10a__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10a__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10b(L_0: float = None, R: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10b__L_0(R: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10b__R(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return result

    @staticmethod
    def eqn_5_10b__V_1(L_0: float, R: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10c(D: float = None, L_0: float = None, R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10c__D(L_0: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10c__L_0(D: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10c__R(D: float, L_0: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
        result.append(R)
        return result

    @kwasak_static
    def eqn_5_11(B: float = None, L_N: float = None, V_0: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_11__B(L_N: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        B = L_N - V_0
        result.append(B)
        return result

    @staticmethod
    def eqn_5_11__L_N(B: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        L_N = B + V_0
        result.append(L_N)
        return result

    @staticmethod
    def eqn_5_11__V_0(B: float, L_N: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        V_0 = -B + L_N
        result.append(V_0)
        return result

    @kwasak_static
    def eqn_5_12(Eff: float = None, N_ES: float = None, N_t: float = None, T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_12__Eff(N_ES: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result

    @staticmethod
    def eqn_5_12__N_ES(Eff: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result

    @staticmethod
    def eqn_5_12__N_t(Eff: float, N_ES: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result

    @staticmethod
    def eqn_5_12__T(Eff: float, N_ES: float, N_t: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result

    @kwasak_static
    def eqn_5_13(HETP: float = None, H_p: float = None, N_ES: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_13__HETP(H_p: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result

    @staticmethod
    def eqn_5_13__H_p(HETP: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result

    @staticmethod
    def eqn_5_13__N_ES(HETP: float, H_p: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result

    @kwasak_static
    def eqn_5_14(M: float = None, P_0: float = None, T: float = None, W_E: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_14__M(P_0: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return result

    @staticmethod
    def eqn_5_14__P_0(M: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return result

    @staticmethod
    def eqn_5_14__T(M: float, P_0: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return result

    @staticmethod
    def eqn_5_14__W_E(M: float, P_0: float, T: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return result

    @kwasak_static
    def eqn_5_15(M_1: float = None, M_2: float = None, P_0_1: float = None, P_0_2: float = None, a_M_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_15__M_1(M_2: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        M_1 = M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        return result

    @staticmethod
    def eqn_5_15__M_2(M_1: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_2 = -M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        M_2 = M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        return result

    @staticmethod
    def eqn_5_15__P_0_1(M_1: float, M_2: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_15__P_0_2(M_1: float, M_2: float, P_0_1: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1*(M_2/M_1)**(2/5)/a_M_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_15__a_M_12(M_1: float, M_2: float, P_0_1: float, P_0_2: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
        result.append(a_M_12)
        return result

    @kwasak_static
    def eqn_5_16(H_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_16__H_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result

    @staticmethod
    def eqn_5_16__p_i(H_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_16__x_i(H_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_17(H_2_1: float = None, H_2_3: float = None, H_2_mi: float = None, x_1: float = None, x_3: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_17__H_2_1(H_2_3: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3*log(H_2_3) + log(H_2_mi))/x_1)
        result.append(H_2_1)
        return result

    @staticmethod
    def eqn_5_17__H_2_3(H_2_1: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1*log(H_2_1) + log(H_2_mi))/x_3)
        result.append(H_2_3)
        return result

    @staticmethod
    def eqn_5_17__H_2_mi(H_2_1: float, H_2_3: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1*log(H_2_1) + x_3*log(H_2_3))
        result.append(H_2_mi)
        return result

    @staticmethod
    def eqn_5_17__x_1(H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_17__x_3(H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return result


class ProcessAppIi:

    @kwasak_static
    def eqn_6_01(T_1: float = None, T_2: float = None, T_R: float = None, c_p: float = None, del_h_v: float = None, w_1: float = None, w_2: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_01__T_1(T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R) + del_h_v*w_v)/(c_p*w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_01__T_2(T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_2 = (T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R) + del_h_v*w_v)/(c_p*w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_01__T_R(T_1: float, T_2: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_R = (T_1*c_p*w_1 + T_2*c_p*w_2 - del_h_v*w_v)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_01__c_p(T_1: float, T_2: float, T_R: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = del_h_v*w_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_01__del_h_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/w_v
        result.append(del_h_v)
        return result

    @staticmethod
    def eqn_6_01__w_1(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (-T_2*c_p*w_2 + T_R*c_p*w_2 + del_h_v*w_v)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_01__w_2(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_2 = (-T_1*c_p*w_1 + T_R*c_p*w_1 + del_h_v*w_v)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result

    @staticmethod
    def eqn_6_01__w_v(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/del_h_v
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_02(Q_v: float = None, T_1: float = None, T_2: float = None, T_R: float = None, c_p: float = None, w_1: float = None, w_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_02__Q_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_02__T_1(Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000*Q_v + T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R))/(c_p*w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_02__T_2(Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_2 = (12000*Q_v + T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R))/(c_p*w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_02__T_R(Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_R = (-12000*Q_v + T_1*c_p*w_1 + T_2*c_p*w_2)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_02__c_p(Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        c_p = 12000*Q_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_02__w_1(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_1 = (12000*Q_v - T_2*c_p*w_2 + T_R*c_p*w_2)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_02__w_2(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000*Q_v - T_1*c_p*w_1 + T_R*c_p*w_1)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result

    @kwasak_static
    def eqn_6_04(Q_v: float = None, delta_h_v: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_04__Q_v(delta_h_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        Q_v = delta_h_v*w_v/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_04__delta_h_v(Q_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        delta_h_v = 12000*Q_v/w_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_04__w_v(Q_v: float, delta_h_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_06(Q_r: float = None, delta_h_v: float = None, f_m: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_06__Q_r(delta_h_v: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        Q_r = delta_h_v*f_m/24
        result.append(Q_r)
        return result

    @staticmethod
    def eqn_6_06__delta_h_v(Q_r: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        delta_h_v = 24*Q_r/f_m
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_06__f_m(Q_r: float, delta_h_v: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        f_m = 24*Q_r/delta_h_v
        result.append(f_m)
        return result

    @kwasak_static
    def eqn_6_07(C_1: float = None, C_2: float = None, T_1: float = None, T_2: float = None, c_p: float = None, delta_h_c: float = None, delta_h_v: float = None, m_b: float = None, m_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_07__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result

    @staticmethod
    def eqn_6_07__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result

    @staticmethod
    def eqn_6_07__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*m_v)/(c_p*m_b)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_07__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*m_v)/(c_p*m_b)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_07__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*m_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_07__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*m_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result

    @staticmethod
    def eqn_6_07__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/m_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_07__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_b = delta_h_v*m_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_07__m_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/delta_h_v
        result.append(m_v)
        return result

    @kwasak_static
    def eqn_6_08(C_1: float = None, C_2: float = None, T_1: float = None, T_2: float = None, c_p: float = None, delta_h_c: float = None, delta_h_v: float = None, delta_t: float = None, m_b: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_08__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result

    @staticmethod
    def eqn_6_08__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result

    @staticmethod
    def eqn_6_08__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_08__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_08__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*delta_t*w_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_08__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*delta_t*w_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result

    @staticmethod
    def eqn_6_08__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_t*w_v)
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_08__delta_t(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_t = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*w_v)
        result.append(delta_t)
        return result

    @staticmethod
    def eqn_6_08__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        m_b = delta_h_v*delta_t*w_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_08__w_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        w_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*delta_t)
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_09(A: float = None, dV_dt: float = None, delta_P: float = None, m: float = None, mu: float = None, r: float = None, r_M: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_09__A(dV_dt: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        A = (dV_dt*r_M - sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        A = (dV_dt*r_M + sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        return result

    @staticmethod
    def eqn_6_09__dV_dt(A: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_09__delta_P(A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        delta_P = A*dV_dt*r_M/(A**2 - dV_dt*m*mu*r)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_09__m(A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
        result.append(m)
        return result

    @staticmethod
    def eqn_6_09__mu(A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_09__r(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
        result.append(r)
        return result

    @staticmethod
    def eqn_6_09__r_M(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
        result.append(r_M)
        return result

    @kwasak_static
    def eqn_6_10(A: float = None, dV_dt: float = None, delta_P: float = None, mu: float = None, r_c: float = None, s: float = None, tau: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_10__A(dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt*delta_P**(s - 1)*mu*r_c*tau
        result.append(A)
        return result

    @staticmethod
    def eqn_6_10__dV_dt(A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A*delta_P**(1 - s)/(mu*r_c*tau)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_10__delta_P(A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_10__mu(A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A*delta_P**(1 - s)/(dV_dt*r_c*tau)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_10__r_c(A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
        result.append(r_c)
        return result

    @staticmethod
    def eqn_6_10__s(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A*delta_P/(dV_dt*mu*r_c*tau))/log(delta_P)
        result.append(s)
        return result

    @staticmethod
    def eqn_6_10__tau(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A*delta_P**(1 - s)/(dV_dt*mu*r_c)
        result.append(tau)
        return result

    @kwasak_static
    def eqn_6_11a(A_d: float = None, delta_T: float = None, delta_h_i: float = None, delta_m: float = None, h_d: float = None, m_b: float = None, t_R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_11a__A_d(delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
        result.append(A_d)
        return result

    @staticmethod
    def eqn_6_11a__delta_T(A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_6_11a__delta_h_i(A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d*delta_T*h_d*t_R/(delta_m*m_b)
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_6_11a__delta_m(A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d*delta_T*h_d*t_R/(delta_h_i*m_b)
        result.append(delta_m)
        return result

    @staticmethod
    def eqn_6_11a__h_d(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i*delta_m*m_b/(A_d*delta_T*t_R)
        result.append(h_d)
        return result

    @staticmethod
    def eqn_6_11a__m_b(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_11a__t_R(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        t_R = delta_h_i*delta_m*m_b/(A_d*delta_T*h_d)
        result.append(t_R)
        return result


class Precondensors:

    @kwasak_static
    def eqn_7_01(P: float = None, p_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_01__P(p_i: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_01__p_i(P: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_01__y_i(P: float, p_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_7_02(P_i_0: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_02__P_i_0(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i/x_i
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_02__p_i(P_i_0: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        p_i = P_i_0*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_02__x_i(P_i_0: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        x_i = p_i/P_i_0
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_03(P_i_0: float = None, epsilon_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_03__P_i_0(epsilon_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i/(epsilon_i*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_03__epsilon_i(P_i_0: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i/(P_i_0*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_03__p_i(P_i_0: float, epsilon_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0*epsilon_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_03__x_i(P_i_0: float, epsilon_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i/(P_i_0*epsilon_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_04a(P: float = None, p_c: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04a__P(p_c: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return result

    @staticmethod
    def eqn_7_04a__p_c(P: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_04a__p_nc(P: float, p_c: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04aa(n_i: float = None, n_nc: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04aa__n_i(n_nc: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_i = n_nc*p_i/p_nc
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04aa__n_nc(n_i: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_nc = n_i*p_nc/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04aa__p_i(n_i: float, n_nc: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_i = n_i*p_nc/n_nc
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04aa__p_nc(n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_nc = n_nc*p_i/n_i
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04ab(P_c: float = None, p: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ab__P_c(p: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ab__p(P_c: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ab__p_i(P_c: float, p: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = 0
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04ab__p_nc(P_c: float, p: float, p_i: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = -P_c + p
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04ac(P_c: float = None, n_i: float = None, n_nc: float = None, p: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ac__P_c(n_i: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ac__n_i(P_c: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04ac__n_nc(P_c: float, n_i: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04ac__p(P_c: float, n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ac__p_i(P_c: float, n_i: float, n_nc: float, p: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result

    @kwasak_static
    def eqn_7_05(N_i: float = None, N_nc: float = None, P: float = None, P_c: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_05__N_i(N_nc: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
        return result

    @staticmethod
    def eqn_7_05__N_nc(N_i: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return result

    @staticmethod
    def eqn_7_05__P(N_i: float, N_nc: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc*p_i/N_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_05__P_c(N_i: float, N_nc: float, P: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc*p_i/N_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_05__p_i(N_i: float, N_nc: float, P: float, P_c: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i*(P - P_c)/N_nc
        result.append(p_i)
        return result

    @kwasak_static
    def eqn_7_06(M: float = None, P: float = None, P_i_0: float = None, W_air: float = None, W_i: float = None, p_c: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_06__M(P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_06__P(M: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_06__P_i_0(M: float, P: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_06__W_air(M: float, P: float, P_i_0: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*x_i)
        result.append(W_air)
        return result

    @staticmethod
    def eqn_7_06__W_i(M: float, P: float, P_i_0: float, W_air: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*x_i/(29*(P - p_c))
        result.append(W_i)
        return result

    @staticmethod
    def eqn_7_06__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_06__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_07(M: float = None, P: float = None, P_i_0: float = None, W_air: float = None, W_i: float = None, epsilon_i: float = None, p_c: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_07__M(P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_07__P(M: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_07__P_i_0(M: float, P: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*epsilon_i*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_07__W_air(M: float, P: float, P_i_0: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*epsilon_i*x_i)
        result.append(W_air)
        return result

    @staticmethod
    def eqn_7_07__W_i(M: float, P: float, P_i_0: float, W_air: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*epsilon_i*x_i/(29*(P - p_c))
        result.append(W_i)
        return result

    @staticmethod
    def eqn_7_07__epsilon_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_07__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_07__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*epsilon_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_08(L_c: float = None, Q_condensor_heat_duty: float = None, c_p: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_08__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_08__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_08__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_08__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_09(L_c: float = None, Q_condensor_heat_duty: float = None, c_p: float = None, del_T: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_09__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_09__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_09__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_09__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return result

    @staticmethod
    def eqn_7_09__rho(L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_7_10(L_c_P: float = None, Q_condensor_heat_duty: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_10__L_c_P(Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty/(500*del_T)
        result.append(L_c_P)
        return result

    @staticmethod
    def eqn_7_10__Q_condensor_heat_duty(L_c_P: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500*L_c_P*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_10__del_T(L_c_P: float, Q_condensor_heat_duty: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(500*L_c_P)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_11(Q_condensor_heat_duty: float = None, U_v: float = None, V_c: float = None, del_T_LM: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_11__Q_condensor_heat_duty(U_v: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v*V_c*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_11__U_v(Q_condensor_heat_duty: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
        result.append(U_v)
        return result

    @staticmethod
    def eqn_7_11__V_c(Q_condensor_heat_duty: float, U_v: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty/(U_v*del_T_LM)
        result.append(V_c)
        return result

    @staticmethod
    def eqn_7_11__del_T_LM(Q_condensor_heat_duty: float, U_v: float, V_c: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(U_v*V_c)
        result.append(del_T_LM)
        return result

    @kwasak_static
    def eqn_7_12(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_12__A(Q_condensor_heat_duty: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty/(U*del_T)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_12__Q_condensor_heat_duty(A: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A*U*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty/(A*del_T)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty/(A*U)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_14a(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T_LM: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14a__A(Q_condensor_heat_duty: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty/(U*del_T_LM)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14a__Q_condensor_heat_duty(A: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A*U*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14a__U(A: float, Q_condensor_heat_duty: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty/(A*del_T_LM)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14a__del_T_LM(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(A*U)
        result.append(del_T_LM)
        return result

    @kwasak_static
    def eqn_7_14b(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T_1: float = None, del_T_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14b__A(Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14b__Q_condensor_heat_duty(A: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A*U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14b__U(A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty/(A*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14b__del_T_1(A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_1)
        return result

    @staticmethod
    def eqn_7_14b__del_T_2(A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_2)
        return result

    @kwasak_static
    def eqn_7_15(U: float = None, sum_R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_15__U(sum_R: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        U = 1/sum_R
        result.append(U)
        return result

    @staticmethod
    def eqn_7_15__sum_R(U: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        sum_R = 1/U
        result.append(sum_R)
        return result

    @kwasak_static
    def eqn_7_16(D_0: float = None, D_LM: float = None, D_i: float = None, R_f_0: float = None, R_fi: float = None, U_0: float = None, h_0: float = None, h_i: float = None, k_w: float = None, x_w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_16__D_0(D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_f_0*U_0*h_0 - U_0 + h_0)/(U_0*h_0*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_16__D_LM(D_0: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_0*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(D_LM)
        return result

    @staticmethod
    def eqn_7_16__D_i(D_0: float, D_LM: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_0*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_0*x_w + D_LM*R_f_0*U_0*h_0*k_w + D_LM*U_0*k_w - D_LM*h_0*k_w))
        result.append(D_i)
        return result

    @staticmethod
    def eqn_7_16__R_f_0(D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
        result.append(R_f_0)
        return result

    @staticmethod
    def eqn_7_16__R_fi(D_0: float, D_LM: float, D_i: float, R_f_0: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_f_0/D_0 - D_i/(D_0*h_0) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result

    @staticmethod
    def eqn_7_16__U_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_0*h_i*k_w/(D_0*D_LM*R_fi*h_0*h_i*k_w + D_0*D_LM*h_0*k_w + D_0*D_i*h_0*h_i*x_w + D_LM*D_i*R_f_0*h_0*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result

    @staticmethod
    def eqn_7_16__h_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_0 = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_f_0*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_0)
        return result

    @staticmethod
    def eqn_7_16__h_i(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_0*k_w/(D_0*D_LM*R_fi*U_0*h_0*k_w + D_0*D_i*U_0*h_0*x_w + D_LM*D_i*R_f_0*U_0*h_0*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_0*k_w)
        result.append(h_i)
        return result

    @staticmethod
    def eqn_7_16__k_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_0*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(k_w)
        return result

    @staticmethod
    def eqn_7_16__x_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_f_0*k_w/D_0 - D_LM*k_w/(D_0*h_0) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result

    @kwasak_static
    def eqn_7_17(R_0: float = None, R_nc: float = None, h_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_17__R_0(R_nc: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1/h_c
        result.append(R_0)
        return result

    @staticmethod
    def eqn_7_17__R_nc(R_0: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1/h_c
        result.append(R_nc)
        return result

    @staticmethod
    def eqn_7_17__h_c(R_0: float, R_nc: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1/(R_0 - R_nc)
        result.append(h_c)
        return result

    @kwasak_static
    def eqn_7_18(D_0: float = None, D_LM: float = None, D_i: float = None, R_fi: float = None, R_fo: float = None, R_nc: float = None, U_0: float = None, h_c: float = None, h_i: float = None, k_w: float = None, x_w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_18__D_0(D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_fo*U_0*h_c - R_nc*U_0*h_c - U_0 + h_c)/(U_0*h_c*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_18__D_LM(D_0: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_c*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(D_LM)
        return result

    @staticmethod
    def eqn_7_18__D_i(D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_c*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_c*x_w + D_LM*R_fo*U_0*h_c*k_w + D_LM*R_nc*U_0*h_c*k_w + D_LM*U_0*k_w - D_LM*h_c*k_w))
        result.append(D_i)
        return result

    @staticmethod
    def eqn_7_18__R_fi(D_0: float, D_LM: float, D_i: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_fo/D_0 - D_i*R_nc/D_0 - D_i/(D_0*h_c) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result

    @staticmethod
    def eqn_7_18__R_fo(D_0: float, D_LM: float, D_i: float, R_fi: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fo = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_nc - 1/h_c + 1/U_0
        result.append(R_fo)
        return result

    @staticmethod
    def eqn_7_18__R_nc(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_nc = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_fo - 1/h_c + 1/U_0
        result.append(R_nc)
        return result

    @staticmethod
    def eqn_7_18__U_0(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_c*h_i*k_w/(D_0*D_LM*R_fi*h_c*h_i*k_w + D_0*D_LM*h_c*k_w + D_0*D_i*h_c*h_i*x_w + D_LM*D_i*R_fo*h_c*h_i*k_w + D_LM*D_i*R_nc*h_c*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result

    @staticmethod
    def eqn_7_18__h_c(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_c = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_fo*U_0*h_i*k_w + D_LM*D_i*R_nc*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_c)
        return result

    @staticmethod
    def eqn_7_18__h_i(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_c*k_w/(D_0*D_LM*R_fi*U_0*h_c*k_w + D_0*D_i*U_0*h_c*x_w + D_LM*D_i*R_fo*U_0*h_c*k_w + D_LM*D_i*R_nc*U_0*h_c*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_c*k_w)
        result.append(h_i)
        return result

    @staticmethod
    def eqn_7_18__k_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_c*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(k_w)
        return result

    @staticmethod
    def eqn_7_18__x_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_fo*k_w/D_0 - D_LM*R_nc*k_w/D_0 - D_LM*k_w/(D_0*h_c) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result


class SelectingPump:

    @kwasak_static
    def eqn_8_01(NC: float = None, NS: float = None, SCON: float = None, installation_cost: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_01__NC(NS: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return result

    @staticmethod
    def eqn_8_01__NS(NC: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return result

    @staticmethod
    def eqn_8_01__SCON(NC: float, NS: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        # [Sympy Failover]
#     def eqn_8_01__SCON(NC: float, NS: float, installation_cost: float):16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35 - (installation_cost)SCON
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np


class VacuumTheory:

    @kwasak_static
    def eqn_1_03(T: float = None, k: float = None, m: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_03__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_1_03__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_1_03__m(T: float, k: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        m = 3.0*T*k/v**2
        result.append(m)
        return result

    @staticmethod
    def eqn_1_03__v(T: float, k: float, m: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        v = -1.73205080756888*sqrt(T*k/m)
        result.append(v)
        v = 1.73205080756888*sqrt(T*k/m)
        result.append(v)
        return result

    @kwasak_static
    def eqn_1_07(R: float = None, T: float = None, V: float = None, n: float = None, p: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_07__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_07__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_07__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return result

    @staticmethod
    def eqn_1_07__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return result

    @staticmethod
    def eqn_1_07__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return result

    @kwasak_static
    def eqn_1_08(M: float = None, P: float = None, R: float = None, T: float = None, V: float = None, m: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_08__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_08__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_08__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_08__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_08__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return result

    @staticmethod
    def eqn_1_08__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return result

    @kwasak_static
    def eqn_1_09(M: float = None, P: float = None, R: float = None, T: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_09__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result

    @staticmethod
    def eqn_1_09__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return result

    @staticmethod
    def eqn_1_09__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_09__T(M: float, P: float, R: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_09__rho(M: float, P: float, R: float, T: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_1_10(P_1: float = None, P_2: float = None, T_1: float = None, T_2: float = None, V_1: float = None, V_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return result

    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return result

    @kwasak_static
    def eqn_1_11(M: float = None, P: float = None, T: float = None, W: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return result

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result

    @kwasak_static
    def eqn_1_12(Total_P: float = None, sum_partial_pressures: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return result

    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return result

    @kwasak_static
    def eqn_1_13a(n: float = None, n_a: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result

    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result

    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result

    @kwasak_static
    def eqn_1_13b(P: float = None, p_a: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return result

    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return result

    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return result


class FluidFlowVacuumLines:

    @kwasak_static
    def eqn_2_01(D: float = None, Re: float = None, mu: float = None, rho: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_01__D(Re: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_01__Re(D: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_01__mu(D: float, Re: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_01__rho(D: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_01__v(D: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_02(delta: float = None, lambd: float = None, psi: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_02__delta(lambd: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        delta = -0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        delta = 0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        return result

    @staticmethod
    def eqn_2_02__lambd(delta: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        lambd = 4.44288293815837*delta**2*psi
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_02__psi(delta: float, lambd: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        psi = 0.225079079039277*lambd/delta**2
        result.append(psi)
        return result

    @kwasak_static
    def eqn_2_03(D: float = None, kn: float = None, lambd: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_03__D(kn: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return result

    @staticmethod
    def eqn_2_03__kn(D: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return result

    @staticmethod
    def eqn_2_03__lambd(D: float, kn: float):
        # [.pyeqn] kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return result

    @kwasak_static
    def eqn_2_04(_beta: float = None, mu: float = None, vel_grad: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_04___beta(mu: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return result

    @staticmethod
    def eqn_2_04__mu(_beta: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_04__vel_grad(_beta: float, mu: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return result

    @kwasak_static
    def eqn_2_05(D: float = None, L: float = None, delta_P: float = None, mu: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_05__D(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        D = -2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = -2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_05__L(D: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703*D**4*delta_P/(mu*q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_05__delta_P(D: float, L: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252*L*mu*q/D**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_05__mu(D: float, L: float, delta_P: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703*D**4*delta_P/(L*q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_05__q(D: float, L: float, delta_P: float, mu: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703*D**4*delta_P/(L*mu)
        result.append(q)
        return result

    @kwasak_static
    def eqn_2_06(lambd: float = None, mu: float = None, rho: float = None, v_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_06__lambd(mu: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_06__mu(lambd: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_06__rho(lambd: float, mu: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_06__v_a(lambd: float, mu: float, rho: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return result

    @kwasak_static
    def eqn_2_07(T: float = None, k: float = None, m: float = None, v_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_07__T(k: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        T = 0.392699081698724*m*v_a**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_2_07__k(T: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        k = 0.392699081698724*m*v_a**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_2_07__m(T: float, k: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return result

    @staticmethod
    def eqn_2_07__v_a(T: float, k: float, m: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return result

    @kwasak_static
    def eqn_2_08(M: float = None, P_c: float = None, T_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return result

    @staticmethod
    def eqn_2_08__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result

    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

    @staticmethod
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_08(M: float = None, P_c: float = None, T_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return result

    @staticmethod
    def eqn_2_08__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result

    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

    @staticmethod
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_10(Suc_Pres: float = None, delta_P: float = None, oper_press: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_10__Suc_Pres(delta_P: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        Suc_Pres = -delta_P + oper_press
        result.append(Suc_Pres)
        return result

    @staticmethod
    def eqn_2_10__delta_P(Suc_Pres: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        delta_P = -Suc_Pres + oper_press
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return result

    @kwasak_static
    def eqn_2_11(D: float = None, L: float = None, f: float = None, g_c: float = None, h_r: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
        return result

    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_11__g_c(D: float, L: float, f: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        g_c = L*f*v**2/(2*D*h_r)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        h_r = L*f*v**2/(2*D*g_c)
        result.append(h_r)
        return result

    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_12(L: float = None, d: float = None, delta_P: float = None, f: float = None, g: float = None, rho: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838*d*delta_P*g/(f*rho*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155*L*f*rho*v**2/(delta_P*g)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_12__delta_P(L: float, d: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        delta_P = 2.155*L*f*rho*v**2/(d*g)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_12__f(L: float, d: float, delta_P: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        f = 0.464037122969838*d*delta_P*g/(L*rho*v**2)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155*L*f*rho*v**2/(d*delta_P)
        result.append(g)
        return result

    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838*d*delta_P*g/(L*f*v**2)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_12__v(L: float, d: float, delta_P: float, f: float, g: float, rho: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        v = -0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        v = 0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_13(L: float = None, d: float = None, delta_P: float = None, f: float = None, q: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_13__d(L: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        d = 1.16543402167043*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) - 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) + 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) - 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) + 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_13__delta_P(L: float, d: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        delta_P = 2.15*L*f*q**2*rho/d**5
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_13__f(L: float, d: float, delta_P: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        f = 0.465116279069767*d**5*delta_P/(L*q**2*rho)
        result.append(f)
        return result

    @staticmethod
    def eqn_2_13__q(L: float, d: float, delta_P: float, f: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        return result

    @staticmethod
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_2_14(M: float = None, R: float = None, T: float = None, g_c: float = None, k: float = None, v_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return result

    @staticmethod
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return result

    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return result

    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return result

    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return result

    @staticmethod
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return result

    @kwasak_static
    def eqn_2_15(Re: float = None, f: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_15__Re(f: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736/f**4
        result.append(Re)
        return result

    @staticmethod
    def eqn_2_15__f(Re: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316/Re**(1/4)
        result.append(f)
        return result

    @kwasak_static
    def eqn_2_17(L: float = None, d: float = None, delta_P: float = None, mu: float = None, v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        d = -0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        d = 0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        delta_P = 0.0345*L*mu*v/d**2
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        mu = 28.9855072463768*d**2*delta_P/(L*v)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__v(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
        return result

    @kwasak_static
    def eqn_2_17(L: float = None, d: float = None, delta_P: float = None, mu: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        L = 9.52380952380952*d**4*delta_P/(mu*q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        d = -0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = -0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        delta_P = 0.105*L*mu*q/d**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        mu = 9.52380952380952*d**4*delta_P/(L*q)
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_17__q(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952*d**4*delta_P/(L*mu)
        result.append(q)
        return result


class PressMgmt:

    @kwasak_static
    def eqn_3_01(Abs_Pressure: float = None, BarometricPressure: float = None, Vacuum: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_01__Abs_Pressure(BarometricPressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return result

    @staticmethod
    def eqn_3_01__BarometricPressure(Abs_Pressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return result

    @staticmethod
    def eqn_3_01__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return result

    @kwasak_static
    def eqn_3_02(G: float = None, G_C: float = None, H: float = None, P: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_02__G(G_C: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return result

    @staticmethod
    def eqn_3_02__G_C(G: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return result

    @staticmethod
    def eqn_3_02__H(G: float, G_C: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return result

    @staticmethod
    def eqn_3_02__P(G: float, G_C: float, H: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_02__rho(G: float, G_C: float, H: float, P: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_3_03(H_1: float = None, H_2: float = None, P: float = None, P_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_03__H_1(H_2: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_03__H_2(H_1: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 - P + P_P
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_03__P(H_1: float, H_2: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P = H_1 - H_2 + P_P
        result.append(P)
        return result

    @staticmethod
    def eqn_3_03__P_P(H_1: float, H_2: float, P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return result

    @kwasak_static
    def eqn_3_04(KAPPA: float = None, P: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_04__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return result

    @staticmethod
    def eqn_3_04__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_04__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return result

    @kwasak_static
    def eqn_3_05(P: float = None, P_P: float = None, V: float = None, V_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_05__P(P_P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_05__P_P(P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return result

    @staticmethod
    def eqn_3_05__V(P: float, P_P: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_05__V_P(P: float, P_P: float, V: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return result

    @kwasak_static
    def eqn_3_06(H_1: float = None, H_2: float = None, P: float = None, V: float = None, V_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_06__H_1(H_2: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_06__H_2(H_1: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P*V/V_P - P
        result.append(H_2)
        return result

    @staticmethod
    def eqn_3_06__P(H_1: float, H_2: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P*(-H_1 + H_2)/(V - V_P)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_06__V(H_1: float, H_2: float, P: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_06__V_P(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return result


class AirLeak:

    @kwasak_static
    def eqn_4_07(W: float = None, W_T: float = None, sum_individual_leak_rates: float = None,**kwargs):
        return


    @staticmethod
    def eqn_4_07__W(W_T: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return result

    @staticmethod
    def eqn_4_07__W_T(W: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return result

    @staticmethod
    def eqn_4_07__sum_individual_leak_rates(W: float, W_T: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result

    @kwasak_static
    def eqn_4_10(T: float = None, V: float = None, del_P: float = None, leakage: float = None, t: float = None,**kwargs):
        return


    @staticmethod
    def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        T = 3.127*V*del_P/(leakage*t)
        result.append(T)
        return result

    @staticmethod
    def eqn_4_10__V(T: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        V = 0.319795330988168*T*leakage*t/del_P
        result.append(V)
        return result

    @staticmethod
    def eqn_4_10__del_P(T: float, V: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        del_P = 0.319795330988168*T*leakage*t/V
        result.append(del_P)
        return result

    @staticmethod
    def eqn_4_10__leakage(T: float, V: float, del_P: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        leakage = 3.127*V*del_P/(T*t)
        result.append(leakage)
        return result

    @staticmethod
    def eqn_4_10__t(T: float, V: float, del_P: float, leakage: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        t = 3.127*V*del_P/(T*leakage)
        result.append(t)
        return result


class ProcessAppI:

    @kwasak_static
    def eqn_5_01(K_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_01__K_i(x_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i/x_i
        result.append(K_i)
        return result

    @staticmethod
    def eqn_5_01__x_i(K_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_01__y_i(K_i: float, x_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_02a(K_1: float = None, K_2: float = None, alpha_1_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02a__K_1(K_2: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2*alpha_1_2
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02a__K_2(K_1: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1/alpha_1_2
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02a__alpha_1_2(K_1: float, K_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1/K_2
        result.append(alpha_1_2)
        return result

    @kwasak_static
    def eqn_5_02b(K_1: float = None, K_2: float = None, x_1: float = None, x_2: float = None, y_1: float = None, y_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02b__K_1(K_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02b__K_2(K_1: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02b__x_1(K_1: float, K_2: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_02b__x_2(K_1: float, K_2: float, x_1: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result

    @staticmethod
    def eqn_5_02b__y_1(K_1: float, K_2: float, x_1: float, x_2: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result

    @staticmethod
    def eqn_5_02b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
        result.append(y_2)
        return result

    @kwasak_static
    def eqn_5_03(P_0_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_03__P_0_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_03__p_i(P_0_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_03__x_i(P_0_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_04(P: float = None, P_0_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_04__P(P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_04__P_0_i(P: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P_0_i = P*y_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_04__x_i(P: float, P_0_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        x_i = P*y_i/P_0_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_04__y_i(P: float, P_0_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_05(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_05__P_0_1(P_0_2: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2*alpha_12
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_05__P_0_2(P_0_1: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1/alpha_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_05__alpha_12(P_0_1: float, P_0_2: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return result

    @kwasak_static
    def eqn_5_06(P_0_i: float = None, gamma_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_06__P_0_i(gamma_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_06__gamma_i(P_0_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_06__p_i(P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_06__x_i(P_0_i: float, gamma_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_07(P: float = None, P_0_i: float = None, gamma_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_07__P(P_0_i: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i*gamma_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_07__P_0_i(P: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P*y_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_07__gamma_i(P: float, P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P*y_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_07__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_07__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_08(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None, gamma_1: float = None, gamma_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_08__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_08__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_08__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result

    @staticmethod
    def eqn_5_08__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return result

    @staticmethod
    def eqn_5_08__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return result

    @kwasak_static
    def eqn_5_09(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_09__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_09__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_09__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10a(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10a__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10a__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10a__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10b(L_0: float = None, R: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10b__L_0(R: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10b__R(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return result

    @staticmethod
    def eqn_5_10b__V_1(L_0: float, R: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10c(D: float = None, L_0: float = None, R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10c__D(L_0: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10c__L_0(D: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10c__R(D: float, L_0: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
        result.append(R)
        return result

    @kwasak_static
    def eqn_5_11(B: float = None, L_N: float = None, V_0: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_11__B(L_N: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        B = L_N - V_0
        result.append(B)
        return result

    @staticmethod
    def eqn_5_11__L_N(B: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        L_N = B + V_0
        result.append(L_N)
        return result

    @staticmethod
    def eqn_5_11__V_0(B: float, L_N: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        V_0 = -B + L_N
        result.append(V_0)
        return result

    @kwasak_static
    def eqn_5_12(Eff: float = None, N_ES: float = None, N_t: float = None, T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_12__Eff(N_ES: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result

    @staticmethod
    def eqn_5_12__N_ES(Eff: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result

    @staticmethod
    def eqn_5_12__N_t(Eff: float, N_ES: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result

    @staticmethod
    def eqn_5_12__T(Eff: float, N_ES: float, N_t: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result

    @kwasak_static
    def eqn_5_13(HETP: float = None, H_p: float = None, N_ES: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_13__HETP(H_p: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result

    @staticmethod
    def eqn_5_13__H_p(HETP: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result

    @staticmethod
    def eqn_5_13__N_ES(HETP: float, H_p: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result

    @kwasak_static
    def eqn_5_14(M: float = None, P_0: float = None, T: float = None, W_E: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_14__M(P_0: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return result

    @staticmethod
    def eqn_5_14__P_0(M: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return result

    @staticmethod
    def eqn_5_14__T(M: float, P_0: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return result

    @staticmethod
    def eqn_5_14__W_E(M: float, P_0: float, T: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return result

    @kwasak_static
    def eqn_5_15(M_1: float = None, M_2: float = None, P_0_1: float = None, P_0_2: float = None, a_M_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_15__M_1(M_2: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        M_1 = M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        return result

    @staticmethod
    def eqn_5_15__M_2(M_1: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_2 = -M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        M_2 = M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        return result

    @staticmethod
    def eqn_5_15__P_0_1(M_1: float, M_2: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_15__P_0_2(M_1: float, M_2: float, P_0_1: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1*(M_2/M_1)**(2/5)/a_M_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_15__a_M_12(M_1: float, M_2: float, P_0_1: float, P_0_2: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
        result.append(a_M_12)
        return result

    @kwasak_static
    def eqn_5_16(H_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_16__H_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result

    @staticmethod
    def eqn_5_16__p_i(H_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_16__x_i(H_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_17(H_2_1: float = None, H_2_3: float = None, H_2_mi: float = None, x_1: float = None, x_3: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_17__H_2_1(H_2_3: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3*log(H_2_3) + log(H_2_mi))/x_1)
        result.append(H_2_1)
        return result

    @staticmethod
    def eqn_5_17__H_2_3(H_2_1: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1*log(H_2_1) + log(H_2_mi))/x_3)
        result.append(H_2_3)
        return result

    @staticmethod
    def eqn_5_17__H_2_mi(H_2_1: float, H_2_3: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1*log(H_2_1) + x_3*log(H_2_3))
        result.append(H_2_mi)
        return result

    @staticmethod
    def eqn_5_17__x_1(H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_17__x_3(H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return result


class ProcessAppIi:

    @kwasak_static
    def eqn_6_01(T_1: float = None, T_2: float = None, T_R: float = None, c_p: float = None, del_h_v: float = None, w_1: float = None, w_2: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_01__T_1(T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R) + del_h_v*w_v)/(c_p*w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_01__T_2(T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_2 = (T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R) + del_h_v*w_v)/(c_p*w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_01__T_R(T_1: float, T_2: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_R = (T_1*c_p*w_1 + T_2*c_p*w_2 - del_h_v*w_v)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_01__c_p(T_1: float, T_2: float, T_R: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = del_h_v*w_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_01__del_h_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/w_v
        result.append(del_h_v)
        return result

    @staticmethod
    def eqn_6_01__w_1(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (-T_2*c_p*w_2 + T_R*c_p*w_2 + del_h_v*w_v)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_01__w_2(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_2 = (-T_1*c_p*w_1 + T_R*c_p*w_1 + del_h_v*w_v)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result

    @staticmethod
    def eqn_6_01__w_v(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/del_h_v
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_02(Q_v: float = None, T_1: float = None, T_2: float = None, T_R: float = None, c_p: float = None, w_1: float = None, w_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_02__Q_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_02__T_1(Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000*Q_v + T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R))/(c_p*w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_02__T_2(Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_2 = (12000*Q_v + T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R))/(c_p*w_2)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_02__T_R(Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_R = (-12000*Q_v + T_1*c_p*w_1 + T_2*c_p*w_2)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_02__c_p(Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        c_p = 12000*Q_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_02__w_1(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_1 = (12000*Q_v - T_2*c_p*w_2 + T_R*c_p*w_2)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_02__w_2(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000*Q_v - T_1*c_p*w_1 + T_R*c_p*w_1)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result

    @kwasak_static
    def eqn_6_04(Q_v: float = None, delta_h_v: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_04__Q_v(delta_h_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        Q_v = delta_h_v*w_v/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_04__delta_h_v(Q_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        delta_h_v = 12000*Q_v/w_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_04__w_v(Q_v: float, delta_h_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_06(Q_r: float = None, delta_h_v: float = None, f_m: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_06__Q_r(delta_h_v: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        Q_r = delta_h_v*f_m/24
        result.append(Q_r)
        return result

    @staticmethod
    def eqn_6_06__delta_h_v(Q_r: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        delta_h_v = 24*Q_r/f_m
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_06__f_m(Q_r: float, delta_h_v: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        f_m = 24*Q_r/delta_h_v
        result.append(f_m)
        return result

    @kwasak_static
    def eqn_6_07(C_1: float = None, C_2: float = None, T_1: float = None, T_2: float = None, c_p: float = None, delta_h_c: float = None, delta_h_v: float = None, m_b: float = None, m_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_07__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result

    @staticmethod
    def eqn_6_07__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result

    @staticmethod
    def eqn_6_07__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*m_v)/(c_p*m_b)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_07__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*m_v)/(c_p*m_b)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_07__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*m_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_07__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*m_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result

    @staticmethod
    def eqn_6_07__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/m_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_07__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_b = delta_h_v*m_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_07__m_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/delta_h_v
        result.append(m_v)
        return result

    @kwasak_static
    def eqn_6_08(C_1: float = None, C_2: float = None, T_1: float = None, T_2: float = None, c_p: float = None, delta_h_c: float = None, delta_h_v: float = None, delta_t: float = None, m_b: float = None, w_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_08__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result

    @staticmethod
    def eqn_6_08__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result

    @staticmethod
    def eqn_6_08__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_08__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_6_08__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*delta_t*w_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_08__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*delta_t*w_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result

    @staticmethod
    def eqn_6_08__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_t*w_v)
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_08__delta_t(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_t = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*w_v)
        result.append(delta_t)
        return result

    @staticmethod
    def eqn_6_08__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        m_b = delta_h_v*delta_t*w_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_08__w_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        w_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*delta_t)
        result.append(w_v)
        return result

    @kwasak_static
    def eqn_6_09(A: float = None, dV_dt: float = None, delta_P: float = None, m: float = None, mu: float = None, r: float = None, r_M: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_09__A(dV_dt: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        A = (dV_dt*r_M - sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        A = (dV_dt*r_M + sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        return result

    @staticmethod
    def eqn_6_09__dV_dt(A: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_09__delta_P(A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        delta_P = A*dV_dt*r_M/(A**2 - dV_dt*m*mu*r)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_09__m(A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
        result.append(m)
        return result

    @staticmethod
    def eqn_6_09__mu(A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_09__r(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
        result.append(r)
        return result

    @staticmethod
    def eqn_6_09__r_M(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
        result.append(r_M)
        return result

    @kwasak_static
    def eqn_6_10(A: float = None, dV_dt: float = None, delta_P: float = None, mu: float = None, r_c: float = None, s: float = None, tau: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_10__A(dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt*delta_P**(s - 1)*mu*r_c*tau
        result.append(A)
        return result

    @staticmethod
    def eqn_6_10__dV_dt(A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A*delta_P**(1 - s)/(mu*r_c*tau)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_10__delta_P(A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_10__mu(A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A*delta_P**(1 - s)/(dV_dt*r_c*tau)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_10__r_c(A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
        result.append(r_c)
        return result

    @staticmethod
    def eqn_6_10__s(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A*delta_P/(dV_dt*mu*r_c*tau))/log(delta_P)
        result.append(s)
        return result

    @staticmethod
    def eqn_6_10__tau(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A*delta_P**(1 - s)/(dV_dt*mu*r_c)
        result.append(tau)
        return result

    @kwasak_static
    def eqn_6_11a(A_d: float = None, delta_T: float = None, delta_h_i: float = None, delta_m: float = None, h_d: float = None, m_b: float = None, t_R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_11a__A_d(delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
        result.append(A_d)
        return result

    @staticmethod
    def eqn_6_11a__delta_T(A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_6_11a__delta_h_i(A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d*delta_T*h_d*t_R/(delta_m*m_b)
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_6_11a__delta_m(A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d*delta_T*h_d*t_R/(delta_h_i*m_b)
        result.append(delta_m)
        return result

    @staticmethod
    def eqn_6_11a__h_d(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i*delta_m*m_b/(A_d*delta_T*t_R)
        result.append(h_d)
        return result

    @staticmethod
    def eqn_6_11a__m_b(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_11a__t_R(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        t_R = delta_h_i*delta_m*m_b/(A_d*delta_T*h_d)
        result.append(t_R)
        return result


class Precondensors:

    @kwasak_static
    def eqn_7_01(P: float = None, p_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_01__P(p_i: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_01__p_i(P: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_01__y_i(P: float, p_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_7_02(P_i_0: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_02__P_i_0(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i/x_i
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_02__p_i(P_i_0: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        p_i = P_i_0*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_02__x_i(P_i_0: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        x_i = p_i/P_i_0
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_03(P_i_0: float = None, epsilon_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_03__P_i_0(epsilon_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i/(epsilon_i*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_03__epsilon_i(P_i_0: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i/(P_i_0*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_03__p_i(P_i_0: float, epsilon_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0*epsilon_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_03__x_i(P_i_0: float, epsilon_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i/(P_i_0*epsilon_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_04a(P: float = None, p_c: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04a__P(p_c: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return result

    @staticmethod
    def eqn_7_04a__p_c(P: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_04a__p_nc(P: float, p_c: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04aa(n_i: float = None, n_nc: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04aa__n_i(n_nc: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_i = n_nc*p_i/p_nc
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04aa__n_nc(n_i: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_nc = n_i*p_nc/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04aa__p_i(n_i: float, n_nc: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_i = n_i*p_nc/n_nc
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04aa__p_nc(n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_nc = n_nc*p_i/n_i
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04ab(P_c: float = None, p: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ab__P_c(p: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ab__p(P_c: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ab__p_i(P_c: float, p: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = 0
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04ab__p_nc(P_c: float, p: float, p_i: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = -P_c + p
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_04ac(P_c: float = None, n_i: float = None, n_nc: float = None, p: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ac__P_c(n_i: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ac__n_i(P_c: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04ac__n_nc(P_c: float, n_i: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04ac__p(P_c: float, n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ac__p_i(P_c: float, n_i: float, n_nc: float, p: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result

    @kwasak_static
    def eqn_7_05(N_i: float = None, N_nc: float = None, P: float = None, P_c: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_05__N_i(N_nc: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
        return result

    @staticmethod
    def eqn_7_05__N_nc(N_i: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return result

    @staticmethod
    def eqn_7_05__P(N_i: float, N_nc: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc*p_i/N_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_05__P_c(N_i: float, N_nc: float, P: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc*p_i/N_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_05__p_i(N_i: float, N_nc: float, P: float, P_c: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i*(P - P_c)/N_nc
        result.append(p_i)
        return result

    @kwasak_static
    def eqn_7_06(M: float = None, P: float = None, P_i_0: float = None, W_air: float = None, W_i: float = None, p_c: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_06__M(P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_06__P(M: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_06__P_i_0(M: float, P: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_06__W_air(M: float, P: float, P_i_0: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*x_i)
        result.append(W_air)
        return result

    @staticmethod
    def eqn_7_06__W_i(M: float, P: float, P_i_0: float, W_air: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*x_i/(29*(P - p_c))
        result.append(W_i)
        return result

    @staticmethod
    def eqn_7_06__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_06__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_07(M: float = None, P: float = None, P_i_0: float = None, W_air: float = None, W_i: float = None, epsilon_i: float = None, p_c: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_07__M(P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_07__P(M: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_07__P_i_0(M: float, P: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*epsilon_i*x_i)
        result.append(P_i_0)
        return result

    @staticmethod
    def eqn_7_07__W_air(M: float, P: float, P_i_0: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*epsilon_i*x_i)
        result.append(W_air)
        return result

    @staticmethod
    def eqn_7_07__W_i(M: float, P: float, P_i_0: float, W_air: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*epsilon_i*x_i/(29*(P - p_c))
        result.append(W_i)
        return result

    @staticmethod
    def eqn_7_07__epsilon_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_07__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @staticmethod
    def eqn_7_07__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*epsilon_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_7_08(L_c: float = None, Q_condensor_heat_duty: float = None, c_p: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_08__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_08__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_08__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_08__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_09(L_c: float = None, Q_condensor_heat_duty: float = None, c_p: float = None, del_T: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_09__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_09__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_09__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_09__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return result

    @staticmethod
    def eqn_7_09__rho(L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return result

    @kwasak_static
    def eqn_7_10(L_c_P: float = None, Q_condensor_heat_duty: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_10__L_c_P(Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty/(500*del_T)
        result.append(L_c_P)
        return result

    @staticmethod
    def eqn_7_10__Q_condensor_heat_duty(L_c_P: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500*L_c_P*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_10__del_T(L_c_P: float, Q_condensor_heat_duty: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(500*L_c_P)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_11(Q_condensor_heat_duty: float = None, U_v: float = None, V_c: float = None, del_T_LM: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_11__Q_condensor_heat_duty(U_v: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v*V_c*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_11__U_v(Q_condensor_heat_duty: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
        result.append(U_v)
        return result

    @staticmethod
    def eqn_7_11__V_c(Q_condensor_heat_duty: float, U_v: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty/(U_v*del_T_LM)
        result.append(V_c)
        return result

    @staticmethod
    def eqn_7_11__del_T_LM(Q_condensor_heat_duty: float, U_v: float, V_c: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(U_v*V_c)
        result.append(del_T_LM)
        return result

    @kwasak_static
    def eqn_7_12(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_12__A(Q_condensor_heat_duty: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty/(U*del_T)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_12__Q_condensor_heat_duty(A: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A*U*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty/(A*del_T)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty/(A*U)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_14a(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T_LM: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14a__A(Q_condensor_heat_duty: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty/(U*del_T_LM)
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14a__Q_condensor_heat_duty(A: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A*U*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14a__U(A: float, Q_condensor_heat_duty: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty/(A*del_T_LM)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14a__del_T_LM(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(A*U)
        result.append(del_T_LM)
        return result

    @kwasak_static
    def eqn_7_14b(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T_1: float = None, del_T_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14b__A(Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14b__Q_condensor_heat_duty(A: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A*U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14b__U(A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty/(A*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14b__del_T_1(A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_1)
        return result

    @staticmethod
    def eqn_7_14b__del_T_2(A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_2)
        return result

    @kwasak_static
    def eqn_7_15(U: float = None, sum_R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_15__U(sum_R: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        U = 1/sum_R
        result.append(U)
        return result

    @staticmethod
    def eqn_7_15__sum_R(U: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        sum_R = 1/U
        result.append(sum_R)
        return result

    @kwasak_static
    def eqn_7_16(D_0: float = None, D_LM: float = None, D_i: float = None, R_f_0: float = None, R_fi: float = None, U_0: float = None, h_0: float = None, h_i: float = None, k_w: float = None, x_w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_16__D_0(D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_f_0*U_0*h_0 - U_0 + h_0)/(U_0*h_0*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_16__D_LM(D_0: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_0*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(D_LM)
        return result

    @staticmethod
    def eqn_7_16__D_i(D_0: float, D_LM: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_0*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_0*x_w + D_LM*R_f_0*U_0*h_0*k_w + D_LM*U_0*k_w - D_LM*h_0*k_w))
        result.append(D_i)
        return result

    @staticmethod
    def eqn_7_16__R_f_0(D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
        result.append(R_f_0)
        return result

    @staticmethod
    def eqn_7_16__R_fi(D_0: float, D_LM: float, D_i: float, R_f_0: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_f_0/D_0 - D_i/(D_0*h_0) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result

    @staticmethod
    def eqn_7_16__U_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_0*h_i*k_w/(D_0*D_LM*R_fi*h_0*h_i*k_w + D_0*D_LM*h_0*k_w + D_0*D_i*h_0*h_i*x_w + D_LM*D_i*R_f_0*h_0*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result

    @staticmethod
    def eqn_7_16__h_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_0 = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_f_0*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_0)
        return result

    @staticmethod
    def eqn_7_16__h_i(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_0*k_w/(D_0*D_LM*R_fi*U_0*h_0*k_w + D_0*D_i*U_0*h_0*x_w + D_LM*D_i*R_f_0*U_0*h_0*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_0*k_w)
        result.append(h_i)
        return result

    @staticmethod
    def eqn_7_16__k_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_0*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(k_w)
        return result

    @staticmethod
    def eqn_7_16__x_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_f_0*k_w/D_0 - D_LM*k_w/(D_0*h_0) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result

    @kwasak_static
    def eqn_7_17(R_0: float = None, R_nc: float = None, h_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_17__R_0(R_nc: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1/h_c
        result.append(R_0)
        return result

    @staticmethod
    def eqn_7_17__R_nc(R_0: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1/h_c
        result.append(R_nc)
        return result

    @staticmethod
    def eqn_7_17__h_c(R_0: float, R_nc: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1/(R_0 - R_nc)
        result.append(h_c)
        return result

    @kwasak_static
    def eqn_7_18(D_0: float = None, D_LM: float = None, D_i: float = None, R_fi: float = None, R_fo: float = None, R_nc: float = None, U_0: float = None, h_c: float = None, h_i: float = None, k_w: float = None, x_w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_18__D_0(D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_fo*U_0*h_c - R_nc*U_0*h_c - U_0 + h_c)/(U_0*h_c*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_18__D_LM(D_0: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_c*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(D_LM)
        return result

    @staticmethod
    def eqn_7_18__D_i(D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_c*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_c*x_w + D_LM*R_fo*U_0*h_c*k_w + D_LM*R_nc*U_0*h_c*k_w + D_LM*U_0*k_w - D_LM*h_c*k_w))
        result.append(D_i)
        return result

    @staticmethod
    def eqn_7_18__R_fi(D_0: float, D_LM: float, D_i: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_fo/D_0 - D_i*R_nc/D_0 - D_i/(D_0*h_c) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result

    @staticmethod
    def eqn_7_18__R_fo(D_0: float, D_LM: float, D_i: float, R_fi: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fo = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_nc - 1/h_c + 1/U_0
        result.append(R_fo)
        return result

    @staticmethod
    def eqn_7_18__R_nc(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_nc = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_fo - 1/h_c + 1/U_0
        result.append(R_nc)
        return result

    @staticmethod
    def eqn_7_18__U_0(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_c*h_i*k_w/(D_0*D_LM*R_fi*h_c*h_i*k_w + D_0*D_LM*h_c*k_w + D_0*D_i*h_c*h_i*x_w + D_LM*D_i*R_fo*h_c*h_i*k_w + D_LM*D_i*R_nc*h_c*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result

    @staticmethod
    def eqn_7_18__h_c(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_c = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_fo*U_0*h_i*k_w + D_LM*D_i*R_nc*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_c)
        return result

    @staticmethod
    def eqn_7_18__h_i(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_c*k_w/(D_0*D_LM*R_fi*U_0*h_c*k_w + D_0*D_i*U_0*h_c*x_w + D_LM*D_i*R_fo*U_0*h_c*k_w + D_LM*D_i*R_nc*U_0*h_c*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_c*k_w)
        result.append(h_i)
        return result

    @staticmethod
    def eqn_7_18__k_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_c*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(k_w)
        return result

    @staticmethod
    def eqn_7_18__x_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_fo*k_w/D_0 - D_LM*R_nc*k_w/D_0 - D_LM*k_w/(D_0*h_c) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result


class SelectingPump:

    @kwasak_static
    def eqn_8_01(NC: float = None, NS: float = None, SCON: float = None, installation_cost: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_01__NC(NS: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return result

    @staticmethod
    def eqn_8_01__NS(NC: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return result

    @staticmethod
    def eqn_8_01__SCON(NC: float, NS: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        # [Sympy Failover]
