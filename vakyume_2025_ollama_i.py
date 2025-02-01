from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
from numpy import np
import pandas as pd

class VacuumTheory:

    @kwasak_static
    def eqn_1_03(m: float = None, v: float = None, k: float = None, T: float = None,**kwargs):
        return


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

    @staticmethod
    def eqn_1_03__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_1_03__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return result

    @kwasak_static
    def eqn_1_07(T: float = None, n: float = None, R: float = None, V: float = None, p: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_07__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_07__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return result

    @staticmethod
    def eqn_1_07__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_07__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return result

    @staticmethod
    def eqn_1_07__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return result

    @kwasak_static
    def eqn_1_08(P: float = None, T: float = None, m: float = None, M: float = None, R: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_08__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_08__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_08__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return result

    @staticmethod
    def eqn_1_08__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_08__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return result

    @staticmethod
    def eqn_1_08__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return result

    @kwasak_static
    def eqn_1_09(P: float = None, T: float = None, rho: float = None, M: float = None, R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_09__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
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

    @staticmethod
    def eqn_1_09__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result

    @staticmethod
    def eqn_1_09__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result

    @kwasak_static
    def eqn_1_10(V_2: float = None, P_2: float = None, V_1: float = None, P_1: float = None, T_1: float = None, T_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return result

    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return result

    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
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

    @kwasak_static
    def eqn_1_11(W: float = None, P: float = None, T: float = None, M: float = None, q: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
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
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result

    @kwasak_static
    def eqn_1_12(sum_partial_pressures: float = None, Total_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return result

    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return result

    @kwasak_static
    def eqn_1_13a(n_a: float = None, n: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result

    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result

    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result

    @kwasak_static
    def eqn_1_13b(p_a: float = None, P: float = None, y_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return result

    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
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
    def eqn_2_01(v: float = None, rho: float = None, D: float = None, Re: float = None, mu: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_01__v(D: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return result

    @staticmethod
    def eqn_2_01__rho(D: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return result

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

    @kwasak_static
    def eqn_2_02(delta: float = None, psi: float = None, lambd: float = None,**kwargs):
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
    def eqn_2_02__psi(delta: float, lambd: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        psi = 0.225079079039277*lambd/delta**2
        result.append(psi)
        return result

    @staticmethod
    def eqn_2_02__lambd(delta: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        lambd = 4.44288293815837*delta**2*psi
        result.append(lambd)
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
    def eqn_2_04(vel_grad: float = None, mu: float = None, _beta: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_04__vel_grad(_beta: float, mu: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return result

    @staticmethod
    def eqn_2_04__mu(_beta: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return result

    @staticmethod
    def eqn_2_04___beta(mu: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return result

    @kwasak_static
    def eqn_2_05(L: float = None, D: float = None, delta_P: float = None, q: float = None, mu: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_05__L(D: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703*D**4*delta_P/(mu*q)
        result.append(L)
        return result

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
    def eqn_2_05__delta_P(D: float, L: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252*L*mu*q/D**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_05__q(D: float, L: float, delta_P: float, mu: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703*D**4*delta_P/(L*mu)
        result.append(q)
        return result

    @staticmethod
    def eqn_2_05__mu(D: float, L: float, delta_P: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703*D**4*delta_P/(L*q)
        result.append(mu)
        return result

    @kwasak_static
    def eqn_2_06(rho: float = None, lambd: float = None, v_a: float = None, mu: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_06__rho(lambd: float, mu: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_06__lambd(mu: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return result

    @staticmethod
    def eqn_2_06__v_a(lambd: float, mu: float, rho: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return result

    @staticmethod
    def eqn_2_06__mu(lambd: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return result

    @kwasak_static
    def eqn_2_07(m: float = None, T: float = None, k: float = None, v_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_07__m(T: float, k: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return result

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
    def eqn_2_07__v_a(T: float, k: float, m: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return result

    @kwasak_static
    def eqn_2_08(T_c: float = None, M: float = None, P_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

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
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_08(T_c: float = None, M: float = None, P_c: float = None, mu_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_08__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result

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
    def eqn_2_08__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

    @kwasak_static
    def eqn_2_10(oper_press: float = None, Suc_Pres: float = None, delta_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return result

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

    @kwasak_static
    def eqn_2_11(v: float = None, L: float = None, D: float = None, g_c: float = None, h_r: float = None, f: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return result

    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
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
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return result

    @kwasak_static
    def eqn_2_12(g: float = None, v: float = None, d: float = None, rho: float = None, L: float = None, delta_P: float = None, f: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155*L*f*rho*v**2/(d*delta_P)
        result.append(g)
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

    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155*L*f*rho*v**2/(delta_P*g)
        result.append(d)
        return result

    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838*d*delta_P*g/(L*f*v**2)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838*d*delta_P*g/(f*rho*v**2)
        result.append(L)
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

    @kwasak_static
    def eqn_2_13(d: float = None, rho: float = None, L: float = None, delta_P: float = None, f: float = None, q: float = None,**kwargs):
        return


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
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return result

    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
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

    @kwasak_static
    def eqn_2_14(T: float = None, k: float = None, g_c: float = None, M: float = None, R: float = None, v_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return result

    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return result

    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return result

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
    def eqn_2_17(v: float = None, d: float = None, L: float = None, delta_P: float = None, mu: float = None,**kwargs):
        return


    @staticmethod
    def eqn_2_17__v(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
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
    def eqn_2_17__L(d: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
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

    @kwasak_static
    def eqn_2_17(d: float = None, L: float = None, delta_P: float = None, q: float = None, mu: float = None,**kwargs):
        return


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
    def eqn_2_17__L(d: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        L = 9.52380952380952*d**4*delta_P/(mu*q)
        result.append(L)
        return result

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        delta_P = 0.105*L*mu*q/d**4
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_2_17__q(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952*d**4*delta_P/(L*mu)
        result.append(q)
        return result

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        mu = 9.52380952380952*d**4*delta_P/(L*q)
        result.append(mu)
        return result


class PressMgmt:

    @kwasak_static
    def eqn_3_01(Vacuum: float = None, Abs_Pressure: float = None, BarometricPressure: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_01__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return result

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

    @kwasak_static
    def eqn_3_02(G: float = None, P: float = None, G_C: float = None, rho: float = None, H: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_02__G(G_C: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return result

    @staticmethod
    def eqn_3_02__P(G: float, G_C: float, H: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return result

    @staticmethod
    def eqn_3_02__G_C(G: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return result

    @staticmethod
    def eqn_3_02__rho(G: float, G_C: float, H: float, P: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return result

    @staticmethod
    def eqn_3_02__H(G: float, G_C: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return result

    @kwasak_static
    def eqn_3_03(H_2: float = None, P: float = None, H_1: float = None, P_P: float = None,**kwargs):
        return


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
    def eqn_3_03__H_1(H_2: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_03__P_P(H_1: float, H_2: float, P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return result

    @kwasak_static
    def eqn_3_04(P: float = None, KAPPA: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_04__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_04__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return result

    @staticmethod
    def eqn_3_04__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return result

    @kwasak_static
    def eqn_3_05(V_P: float = None, P: float = None, V: float = None, P_P: float = None,**kwargs):
        return


    @staticmethod
    def eqn_3_05__V_P(P: float, P_P: float, V: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return result

    @staticmethod
    def eqn_3_05__P(P_P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return result

    @staticmethod
    def eqn_3_05__V(P: float, P_P: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return result

    @staticmethod
    def eqn_3_05__P_P(P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return result

    @kwasak_static
    def eqn_3_06(H_2: float = None, P: float = None, H_1: float = None, V_P: float = None, V: float = None,**kwargs):
        return


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
    def eqn_3_06__H_1(H_2: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return result

    @staticmethod
    def eqn_3_06__V_P(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return result

    @staticmethod
    def eqn_3_06__V(H_1: float, H_2: float, P: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return result


class AirLeak:

    @kwasak_static
    def eqn_4_07(sum_individual_leak_rates: float = None, W: float = None, W_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_4_07__sum_individual_leak_rates(W: float, W_T: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result

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

    @kwasak_static
    def eqn_4_10(del_P: float = None, leakage: float = None, T: float = None, t: float = None, V: float = None,**kwargs):
        return


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
    def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        T = 3.127*V*del_P/(leakage*t)
        result.append(T)
        return result

    @staticmethod
    def eqn_4_10__t(T: float, V: float, del_P: float, leakage: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        t = 3.127*V*del_P/(T*leakage)
        result.append(t)
        return result

    @staticmethod
    def eqn_4_10__V(T: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        V = 0.319795330988168*T*leakage*t/del_P
        result.append(V)
        return result


class ProcessAppI:

    @kwasak_static
    def eqn_5_01(K_i: float = None, y_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_01__K_i(x_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i/x_i
        result.append(K_i)
        return result

    @staticmethod
    def eqn_5_01__y_i(K_i: float, x_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
        result.append(y_i)
        return result

    @staticmethod
    def eqn_5_01__x_i(K_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
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
    def eqn_5_02b(K_2: float = None, K_1: float = None, y_1: float = None, x_2: float = None, y_2: float = None, x_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02b__K_2(K_1: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02b__K_1(K_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02b__y_1(K_1: float, K_2: float, x_1: float, x_2: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result

    @staticmethod
    def eqn_5_02b__x_2(K_1: float, K_2: float, x_1: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result

    @staticmethod
    def eqn_5_02b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
        result.append(y_2)
        return result

    @staticmethod
    def eqn_5_02b__x_1(K_1: float, K_2: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result

    @kwasak_static
    def eqn_5_03(p_i: float = None, P_0_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_03__p_i(P_0_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_03__P_0_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_03__x_i(P_0_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_04(P: float = None, y_i: float = None, P_0_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_04__P(P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_04__y_i(P: float, P_0_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
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

    @kwasak_static
    def eqn_5_05(alpha_12: float = None, P_0_1: float = None, P_0_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_05__alpha_12(P_0_1: float, P_0_2: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return result

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

    @kwasak_static
    def eqn_5_06(p_i: float = None, P_0_i: float = None, gamma_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_06__p_i(P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return result

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
    def eqn_5_06__x_i(P_0_i: float, gamma_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_07(P: float = None, P_0_i: float = None, gamma_i: float = None, y_i: float = None, x_i: float = None,**kwargs):
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
    def eqn_5_07__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return result

    @staticmethod
    def eqn_5_07__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_08(P_0_1: float = None, alpha_12: float = None, gamma_2: float = None, P_0_2: float = None, gamma_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_08__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_08__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result

    @staticmethod
    def eqn_5_08__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return result

    @staticmethod
    def eqn_5_08__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_08__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
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
    def eqn_5_10c(D: float = None, R: float = None, L_0: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10c__D(L_0: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10c__R(D: float, L_0: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
        result.append(R)
        return result

    @staticmethod
    def eqn_5_10c__L_0(D: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
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
    def eqn_5_12(N_ES: float = None, T: float = None, N_t: float = None, Eff: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_12__N_ES(Eff: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result

    @staticmethod
    def eqn_5_12__T(Eff: float, N_ES: float, N_t: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result

    @staticmethod
    def eqn_5_12__N_t(Eff: float, N_ES: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result

    @staticmethod
    def eqn_5_12__Eff(N_ES: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result

    @kwasak_static
    def eqn_5_13(H_p: float = None, HETP: float = None, N_ES: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_13__H_p(HETP: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result

    @staticmethod
    def eqn_5_13__HETP(H_p: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result

    @staticmethod
    def eqn_5_13__N_ES(HETP: float, H_p: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result

    @kwasak_static
    def eqn_5_14(M: float = None, T: float = None, P_0: float = None, W_E: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_14__M(P_0: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return result

    @staticmethod
    def eqn_5_14__T(M: float, P_0: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return result

    @staticmethod
    def eqn_5_14__P_0(M: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return result

    @staticmethod
    def eqn_5_14__W_E(M: float, P_0: float, T: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return result

    @kwasak_static
    def eqn_5_15(P_0_1: float = None, M_2: float = None, a_M_12: float = None, P_0_2: float = None, M_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_15__P_0_1(M_1: float, M_2: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
        result.append(P_0_1)
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
    def eqn_5_15__a_M_12(M_1: float, M_2: float, P_0_1: float, P_0_2: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
        result.append(a_M_12)
        return result

    @staticmethod
    def eqn_5_15__P_0_2(M_1: float, M_2: float, P_0_1: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1*(M_2/M_1)**(2/5)/a_M_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_15__M_1(M_2: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        M_1 = M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        return result

    @kwasak_static
    def eqn_5_16(x_i: float = None, p_i: float = None, H_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_16__x_i(H_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_16__p_i(H_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_16__H_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result

    @kwasak_static
    def eqn_5_17(H_2_1: float = None, H_2_mi: float = None, x_3: float = None, x_1: float = None, H_2_3: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_17__H_2_1(H_2_3: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3*log(H_2_3) + log(H_2_mi))/x_1)
        result.append(H_2_1)
        return result

    @staticmethod
    def eqn_5_17__H_2_mi(H_2_1: float, H_2_3: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1*log(H_2_1) + x_3*log(H_2_3))
        result.append(H_2_mi)
        return result

    @staticmethod
    def eqn_5_17__x_3(H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return result

    @staticmethod
    def eqn_5_17__x_1(H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_17__H_2_3(H_2_1: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1*log(H_2_1) + log(H_2_mi))/x_3)
        result.append(H_2_3)
        return result


class ProcessAppIi:

    @kwasak_static
    def eqn_6_01(del_h_v: float = None, T_R: float = None, w_2: float = None, w_v: float = None, T_1: float = None, c_p: float = None, w_1: float = None, T_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_01__del_h_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/w_v
        result.append(del_h_v)
        return result

    @staticmethod
    def eqn_6_01__T_R(T_1: float, T_2: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_R = (T_1*c_p*w_1 + T_2*c_p*w_2 - del_h_v*w_v)/(c_p*(w_1 + w_2))
        result.append(T_R)
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

    @staticmethod
    def eqn_6_01__T_1(T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R) + del_h_v*w_v)/(c_p*w_1)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_01__c_p(T_1: float, T_2: float, T_R: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = del_h_v*w_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_01__w_1(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (-T_2*c_p*w_2 + T_R*c_p*w_2 + del_h_v*w_v)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result

    @staticmethod
    def eqn_6_01__T_2(T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_2 = (T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R) + del_h_v*w_v)/(c_p*w_2)
        result.append(T_2)
        return result

    @kwasak_static
    def eqn_6_02(T_R: float = None, w_2: float = None, T_1: float = None, c_p: float = None, w_1: float = None, Q_v: float = None, T_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_02__T_R(Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_R = (-12000*Q_v + T_1*c_p*w_1 + T_2*c_p*w_2)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result

    @staticmethod
    def eqn_6_02__w_2(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000*Q_v - T_1*c_p*w_1 + T_R*c_p*w_1)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result

    @staticmethod
    def eqn_6_02__T_1(Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000*Q_v + T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R))/(c_p*w_1)
        result.append(T_1)
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
    def eqn_6_02__Q_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/12000
        result.append(Q_v)
        return result

    @staticmethod
    def eqn_6_02__T_2(Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_2 = (12000*Q_v + T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R))/(c_p*w_2)
        result.append(T_2)
        return result

    @kwasak_static
    def eqn_6_04(w_v: float = None, Q_v: float = None, delta_h_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_04__w_v(Q_v: float, delta_h_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return result

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

    @kwasak_static
    def eqn_6_06(Q_r: float = None, f_m: float = None, delta_h_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_06__Q_r(delta_h_v: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        Q_r = delta_h_v*f_m/24
        result.append(Q_r)
        return result

    @staticmethod
    def eqn_6_06__f_m(Q_r: float, delta_h_v: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        f_m = 24*Q_r/delta_h_v
        result.append(f_m)
        return result

    @staticmethod
    def eqn_6_06__delta_h_v(Q_r: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        delta_h_v = 24*Q_r/f_m
        result.append(delta_h_v)
        return result

    @kwasak_static
    def eqn_6_07(delta_h_v: float = None, delta_h_c: float = None, m_b: float = None, m_v: float = None, c_p: float = None, T_1: float = None, T_2: float = None, C_1: float = None, C_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_07__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/m_v
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_07__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*m_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
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

    @staticmethod
    def eqn_6_07__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*m_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
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

    @kwasak_static
    def eqn_6_08(delta_h_v: float = None, delta_h_c: float = None, m_b: float = None, w_v: float = None, T_1: float = None, c_p: float = None, delta_t: float = None, T_2: float = None, C_1: float = None, C_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_08__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_t*w_v)
        result.append(delta_h_v)
        return result

    @staticmethod
    def eqn_6_08__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*delta_t*w_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
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

    @staticmethod
    def eqn_6_08__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_6_08__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*delta_t*w_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result

    @staticmethod
    def eqn_6_08__delta_t(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_t = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*w_v)
        result.append(delta_t)
        return result

    @staticmethod
    def eqn_6_08__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_2)
        return result

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

    @kwasak_static
    def eqn_6_09(r: float = None, dV_dt: float = None, r_M: float = None, A: float = None, m: float = None, delta_P: float = None, mu: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_09__r(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
        result.append(r)
        return result

    @staticmethod
    def eqn_6_09__dV_dt(A: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_09__r_M(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
        result.append(r_M)
        return result

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
    def eqn_6_09__m(A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
        result.append(m)
        return result

    @staticmethod
    def eqn_6_09__delta_P(A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        delta_P = A*dV_dt*r_M/(A**2 - dV_dt*m*mu*r)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_09__mu(A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
        result.append(mu)
        return result

    @kwasak_static
    def eqn_6_10(dV_dt: float = None, s: float = None, A: float = None, tau: float = None, delta_P: float = None, r_c: float = None, mu: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_10__dV_dt(A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A*delta_P**(1 - s)/(mu*r_c*tau)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_10__s(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A*delta_P/(dV_dt*mu*r_c*tau))/log(delta_P)
        result.append(s)
        return result

    @staticmethod
    def eqn_6_10__A(dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt*delta_P**(s - 1)*mu*r_c*tau
        result.append(A)
        return result

    @staticmethod
    def eqn_6_10__tau(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A*delta_P**(1 - s)/(dV_dt*mu*r_c)
        result.append(tau)
        return result

    @staticmethod
    def eqn_6_10__delta_P(A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_10__r_c(A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
        result.append(r_c)
        return result

    @staticmethod
    def eqn_6_10__mu(A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A*delta_P**(1 - s)/(dV_dt*r_c*tau)
        result.append(mu)
        return result

    @kwasak_static
    def eqn_6_11a(delta_m: float = None, A_d: float = None, m_b: float = None, delta_h_i: float = None, delta_T: float = None, h_d: float = None, t_R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_11a__delta_m(A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d*delta_T*h_d*t_R/(delta_h_i*m_b)
        result.append(delta_m)
        return result

    @staticmethod
    def eqn_6_11a__A_d(delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
        result.append(A_d)
        return result

    @staticmethod
    def eqn_6_11a__m_b(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
        result.append(m_b)
        return result

    @staticmethod
    def eqn_6_11a__delta_h_i(A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d*delta_T*h_d*t_R/(delta_m*m_b)
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_6_11a__delta_T(A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_6_11a__h_d(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i*delta_m*m_b/(A_d*delta_T*t_R)
        result.append(h_d)
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
    def eqn_7_01(P: float = None, y_i: float = None, p_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_01__P(p_i: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_7_01__y_i(P: float, p_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return result

    @staticmethod
    def eqn_7_01__p_i(P: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
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
    def eqn_7_03(epsilon_i: float = None, P_i_0: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_03__epsilon_i(P_i_0: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i/(P_i_0*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_03__P_i_0(epsilon_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i/(epsilon_i*x_i)
        result.append(P_i_0)
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
    def eqn_7_04ab(P_c: float = None, p_i: float = None, p_nc: float = None, p: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ab__P_c(p: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
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

    @staticmethod
    def eqn_7_04ab__p(P_c: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result

    @kwasak_static
    def eqn_7_04ac(P_c: float = None, n_nc: float = None, p_i: float = None, n_i: float = None, p: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ac__P_c(n_i: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ac__n_nc(P_c: float, n_i: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result

    @staticmethod
    def eqn_7_04ac__p_i(P_c: float, n_i: float, n_nc: float, p: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04ac__n_i(P_c: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result

    @staticmethod
    def eqn_7_04ac__p(P_c: float, n_i: float, n_nc: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result

    @kwasak_static
    def eqn_7_05(N_i: float = None, P: float = None, P_c: float = None, p_i: float = None, N_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_05__N_i(N_nc: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
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

    @staticmethod
    def eqn_7_05__N_nc(N_i: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return result

    @kwasak_static
    def eqn_7_06(P: float = None, W_i: float = None, P_i_0: float = None, W_air: float = None, M: float = None, x_i: float = None, p_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_06__P(M: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_06__W_i(M: float, P: float, P_i_0: float, W_air: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*x_i/(29*(P - p_c))
        result.append(W_i)
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
    def eqn_7_06__M(P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_06__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air)
        result.append(x_i)
        return result

    @staticmethod
    def eqn_7_06__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @kwasak_static
    def eqn_7_07(epsilon_i: float = None, P: float = None, W_i: float = None, P_i_0: float = None, W_air: float = None, M: float = None, x_i: float = None, p_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_07__epsilon_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*x_i)
        result.append(epsilon_i)
        return result

    @staticmethod
    def eqn_7_07__P(M: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + p_c
        result.append(P)
        return result

    @staticmethod
    def eqn_7_07__W_i(M: float, P: float, P_i_0: float, W_air: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*epsilon_i*x_i/(29*(P - p_c))
        result.append(W_i)
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
    def eqn_7_07__M(P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
        result.append(M)
        return result

    @staticmethod
    def eqn_7_07__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*epsilon_i)
        result.append(x_i)
        return result

    @staticmethod
    def eqn_7_07__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + P
        result.append(p_c)
        return result

    @kwasak_static
    def eqn_7_08(L_c: float = None, c_p: float = None, Q_condensor_heat_duty: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_08__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return result

    @staticmethod
    def eqn_7_08__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_08__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_08__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_09(rho: float = None, c_p: float = None, Q_condensor_heat_duty: float = None, del_T: float = None, L_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_09__rho(L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return result

    @staticmethod
    def eqn_7_09__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_7_09__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_09__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return result

    @staticmethod
    def eqn_7_09__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return result

    @kwasak_static
    def eqn_7_10(Q_condensor_heat_duty: float = None, L_c_P: float = None, del_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_10__Q_condensor_heat_duty(L_c_P: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500*L_c_P*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_10__L_c_P(Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty/(500*del_T)
        result.append(L_c_P)
        return result

    @staticmethod
    def eqn_7_10__del_T(L_c_P: float, Q_condensor_heat_duty: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(500*L_c_P)
        result.append(del_T)
        return result

    @kwasak_static
    def eqn_7_11(V_c: float = None, del_T_LM: float = None, U_v: float = None, Q_condensor_heat_duty: float = None,**kwargs):
        return


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

    @staticmethod
    def eqn_7_11__U_v(Q_condensor_heat_duty: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
        result.append(U_v)
        return result

    @staticmethod
    def eqn_7_11__Q_condensor_heat_duty(U_v: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v*V_c*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @kwasak_static
    def eqn_7_12(U: float = None, Q_condensor_heat_duty: float = None, del_T: float = None, A: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty/(A*del_T)
        result.append(U)
        return result

    @staticmethod
    def eqn_7_12__Q_condensor_heat_duty(A: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A*U*del_T
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty/(A*U)
        result.append(del_T)
        return result

    @staticmethod
    def eqn_7_12__A(Q_condensor_heat_duty: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty/(U*del_T)
        result.append(A)
        return result

    @kwasak_static
    def eqn_7_14a(U: float = None, del_T_LM: float = None, Q_condensor_heat_duty: float = None, A: float = None,**kwargs):
        return


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

    @staticmethod
    def eqn_7_14a__Q_condensor_heat_duty(A: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A*U*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14a__A(Q_condensor_heat_duty: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty/(U*del_T_LM)
        result.append(A)
        return result

    @kwasak_static
    def eqn_7_14b(del_T_2: float = None, U: float = None, A: float = None, del_T_1: float = None, Q_condensor_heat_duty: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14b__del_T_2(A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_2)
        return result

    @staticmethod
    def eqn_7_14b__U(A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty/(A*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14b__A(Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14b__del_T_1(A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_1)
        return result

    @staticmethod
    def eqn_7_14b__Q_condensor_heat_duty(A: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A*U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
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
    def eqn_7_16(R_fi: float = None, D_i: float = None, U_0: float = None, D_LM: float = None, h_0: float = None, h_i: float = None, R_f_0: float = None, D_0: float = None, x_w: float = None, k_w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_16__R_fi(D_0: float, D_LM: float, D_i: float, R_f_0: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_f_0/D_0 - D_i/(D_0*h_0) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result

    @staticmethod
    def eqn_7_16__D_i(D_0: float, D_LM: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_0*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_0*x_w + D_LM*R_f_0*U_0*h_0*k_w + D_LM*U_0*k_w - D_LM*h_0*k_w))
        result.append(D_i)
        return result

    @staticmethod
    def eqn_7_16__U_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_0*h_i*k_w/(D_0*D_LM*R_fi*h_0*h_i*k_w + D_0*D_LM*h_0*k_w + D_0*D_i*h_0*h_i*x_w + D_LM*D_i*R_f_0*h_0*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result

    @staticmethod
    def eqn_7_16__D_LM(D_0: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_0*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(D_LM)
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
    def eqn_7_16__R_f_0(D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
        result.append(R_f_0)
        return result

    @staticmethod
    def eqn_7_16__D_0(D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_f_0*U_0*h_0 - U_0 + h_0)/(U_0*h_0*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_16__x_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_f_0*k_w/D_0 - D_LM*k_w/(D_0*h_0) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result

    @staticmethod
    def eqn_7_16__k_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_0*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(k_w)
        return result

    @kwasak_static
    def eqn_7_17(R_nc: float = None, h_c: float = None, R_0: float = None,**kwargs):
        return


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

    @staticmethod
    def eqn_7_17__R_0(R_nc: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1/h_c
        result.append(R_0)
        return result

    @kwasak_static
    def eqn_7_18(R_fi: float = None, R_fo: float = None, D_i: float = None, U_0: float = None, R_nc: float = None, D_LM: float = None, h_i: float = None, D_0: float = None, h_c: float = None, x_w: float = None, k_w: float = None,**kwargs):
        return


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
    def eqn_7_18__D_i(D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_c*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_c*x_w + D_LM*R_fo*U_0*h_c*k_w + D_LM*R_nc*U_0*h_c*k_w + D_LM*U_0*k_w - D_LM*h_c*k_w))
        result.append(D_i)
        return result

    @staticmethod
    def eqn_7_18__U_0(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_c*h_i*k_w/(D_0*D_LM*R_fi*h_c*h_i*k_w + D_0*D_LM*h_c*k_w + D_0*D_i*h_c*h_i*x_w + D_LM*D_i*R_fo*h_c*h_i*k_w + D_LM*D_i*R_nc*h_c*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result

    @staticmethod
    def eqn_7_18__R_nc(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_nc = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_fo - 1/h_c + 1/U_0
        result.append(R_nc)
        return result

    @staticmethod
    def eqn_7_18__D_LM(D_0: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_c*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(D_LM)
        return result

    @staticmethod
    def eqn_7_18__h_i(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_c*k_w/(D_0*D_LM*R_fi*U_0*h_c*k_w + D_0*D_i*U_0*h_c*x_w + D_LM*D_i*R_fo*U_0*h_c*k_w + D_LM*D_i*R_nc*U_0*h_c*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_c*k_w)
        result.append(h_i)
        return result

    @staticmethod
    def eqn_7_18__D_0(D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_fo*U_0*h_c - R_nc*U_0*h_c - U_0 + h_c)/(U_0*h_c*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result

    @staticmethod
    def eqn_7_18__h_c(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_c = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_fo*U_0*h_i*k_w + D_LM*D_i*R_nc*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_c)
        return result

    @staticmethod
    def eqn_7_18__x_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_fo*k_w/D_0 - D_LM*R_nc*k_w/D_0 - D_LM*k_w/(D_0*h_c) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result

    @staticmethod
    def eqn_7_18__k_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_c*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(k_w)
        return result


class SelectingPump:

    @kwasak_static
    def eqn_8_01(NS: float = None, NC: float = None, SCON: float = None, installation_cost: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_01__NS(NC: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return result

    @staticmethod
    def eqn_8_01__NC(NS: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return result

    @staticmethod
    def eqn_8_01__SCON(NC: float, NS: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        # correct [] 
        # [Sympy Failover]
        """
        Solve the given equation for SCON in terms of NC and NS.
        
        Parameters:
        NC (float): The number of cores.
        NS (float): The number of segments.
        installation_cost (float): The cost of installation.
        
        Returns:
        float: The calculated value for SCON.
        """
        
        coefficient = 16000 * (NS + 2 * NC)
        def Cookie_Cutter(NC: float, NS: float, I: float) -> float:
            """
            Calculate the SCON value based on given parameters.
            
            Parameters:
            NC (float): The number of cores (in terms of NCORES).
            NS (float): The number of segments.
            I (float): The installation cost.
            
            Returns:
            float: The SCON value calculated from the given equation.
            """
        
        exponent_base = installation_cost / (1000 ** Cookie_Cutter(NC, NS, I):
        # Calculate (SCON / 1000)**0.35 using the rearranged formula
        term = I / (coefficient)
        power_term = term ** (1/0.35)
        # Raise the power term to the fourth power and multiply by 1000
        SCON = 1000 * power_term ** 4
        
        return [ SCON ]


    @staticmethod
    def eqn_8_01__installation_cost(NC: float, NS: float, SCON: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        installation_cost = 1426.00150101399*SCON**(7/20)*(2.0*NC + NS)
        result.append(installation_cost)
        return result

    @kwasak_static
    def eqn_8_02(installed_costs: float = None, hp: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_02__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        installed_costs = 10435.5162785557*sqrt(hp)
        result.append(installed_costs)
        return result

    @staticmethod
    def eqn_8_02__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        hp = 9.18273645546364e-9*installed_costs**2
        result.append(hp)
        return result

    @kwasak_static
    def eqn_8_03(installed_costs: float = None, hp: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_03__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        result = []
        installed_costs = 13482.9087908759*hp**(9/20)
        result.append(installed_costs)
        return result

    @staticmethod
    def eqn_8_03__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        # [Sympy Failover]
        """Calculate the horsepower (hp) based on the given installed costs using Equation 8-03."""
        hp = 10 * (installed_costs / 38000) ** (1 / 0.45)
        return [ hp ]


    @kwasak_static
    def eqn_8_04(installed_costs: float = None, hp: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_04__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        installed_costs = 10350.7864343909*hp**(2/5)
        result.append(installed_costs)
        return result

    @staticmethod
    def eqn_8_04__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        hp = -9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        hp = 9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        return result

    @kwasak_static
    def eqn_8_05(theoretical_adiabatic_horsepower: float = None, actual_brake_horsepower: float = None, Eff: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_05__theoretical_adiabatic_horsepower(Eff: float, actual_brake_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        theoretical_adiabatic_horsepower = Eff*actual_brake_horsepower
        result.append(theoretical_adiabatic_horsepower)
        return result

    @staticmethod
    def eqn_8_05__actual_brake_horsepower(Eff: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        actual_brake_horsepower = theoretical_adiabatic_horsepower/Eff
        result.append(actual_brake_horsepower)
        return result

    @staticmethod
    def eqn_8_05__Eff(actual_brake_horsepower: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        Eff = theoretical_adiabatic_horsepower/actual_brake_horsepower
        result.append(Eff)
        return result

    @kwasak_static
    def eqn_8_06(adiabatic_hp: float = None, T: float = None, P_2: float = None, k: float = None, P_1: float = None, M: float = None, R: float = None, w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_06__adiabatic_hp(M: float, P_1: float, P_2: float, R: float, T: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        adiabatic_hp = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*M*(k - 1))
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_06__T(M: float, P_1: float, P_2: float, R: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        T = 1980000*M*adiabatic_hp*(k - 1)/(R*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(T)
        return result

    @staticmethod
    def eqn_8_06__P_2(M: float, P_1: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_2 = P_1*(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_2)
        return result

    @staticmethod
    def eqn_8_06__k(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        # [Sympy Failover]
        pass # no closed form solution

    @staticmethod
    def eqn_8_06__P_1(M: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_1 = P_2/(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_1)
        return result

    @staticmethod
    def eqn_8_06__M(P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        M = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*adiabatic_hp*(k - 1))
        result.append(M)
        return result

    @staticmethod
    def eqn_8_06__R(M: float, P_1: float, P_2: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        R = 1980000*M*adiabatic_hp*(k - 1)/(T*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(R)
        return result

    @staticmethod
    def eqn_8_06__w(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        w = 1980000*M*adiabatic_hp*(k - 1)/(R*T*k*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(w)
        return result

    @kwasak_static
    def eqn_8_07(adiabatic_hp: float = None, P_2: float = None, P_1: float = None, w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_07__adiabatic_hp(P_1: float, P_2: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_hp = 0.05*w*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_07__P_2(P_1: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        """
        Calculate P_2 based on the given equation from an engineering problem.
        
        Parameters:
        - P_1 (float): initial pressure in Pascals.
        - adiabatic_hp (float): heat power loss due to adiabatic process in Watts.
        - w (float): weight in kilograms.
        
        Returns:
        - P_2 (float): calculated final pressure in Pascals.
        """
        return [ P_1 * ((adiabatic_hp * 20) / w) ** (-1/0.286) ]


    @staticmethod
    def eqn_8_07__P_1(P_2: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        """
        Solve for P_1 using the given parameters and equation.
        
        :param P_2: The known pressure P_2 (float).
        :param adiabatic_hp: The horsepower (adibatic hp) value (float).
        :param w: The constant related to work done by the system (float).
        :return: The calculated pressure P_1 as a float.
        """
        
        # Calculate the left side of the equation after adding 1 and dividing by (w / 20)
        term = ((adiabatic_hp + 1) * (20 / w)) ** (-1/0.286)
        
        # Raise P_2 to the power of 0.286 as per the reciprocal step in the equation
        adjusted_P_2 = P_2 ** 0inasbbox
        # Calculate and return P_1 using the adjusted term and raised P_2 value
        P_1 = 1 / term * (adjusted_P_2)**(0.286)
        
        return [ P_1 ]


    @staticmethod
    def eqn_8_07__w(P_1: float, P_2: float, adiabatic_hp: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        w = 20.0*adiabatic_hp/((P_2/P_1)**0.286 - 1.0)
        result.append(w)
        return result

    @kwasak_static
    def eqn_8_08(P_1: float = None, f: float = None, P_2: float = None, adiabatic_power_watts: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_08__P_1(P_2: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        """
        Solves the given equation for P_1 in terms of known quantities.
        
        Parameters:
            P_2 (float): The value of P_2 in watts.
            adiabatic_power_watts (float): The power associated with an adiabatic process in watts.
            f (float): A constant related to the specific problem context.
            
        Returns:
            float: The calculated value of P_1 based on the given inputs and equation.
        """
        
        left_side = adiabatic_power_watts / (f / 12) + 1
        numerator = P_2 / (left_side ** (0.286))
        
        # Taking the reciprocal of the numerator and raising to the power of 5/4 to solve for P_1
        P_1 = (numerator ** (2.44444444444)) / f
        
        return [ P_1 ]


    @staticmethod
    def eqn_8_08__f(P_1: float, P_2: float, adiabatic_power_watts: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        f = 12.0*adiabatic_power_watts/((P_2/P_1)**0.286 - 1.0)
        result.append(f)
        return result

    @staticmethod
    def eqn_8_08__P_2(P_1: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        pass # no closed form solution

    @staticmethod
    def eqn_8_08__adiabatic_power_watts(P_1: float, P_2: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_power_watts = 0.0833333333333333*f*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_power_watts)
        return result

    @kwasak_static
    def eqn_8_09(r: float = None, E_m: float = None, E_j: float = None, s: float = None, e: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_09__r(E_j: float, E_m: float, e: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        r = 2.93*E_j*e/(E_m*s)
        result.append(r)
        return result

    @staticmethod
    def eqn_8_09__E_m(E_j: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_m = 2.93*E_j*e/(r*s)
        result.append(E_m)
        return result

    @staticmethod
    def eqn_8_09__E_j(E_m: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_j = 0.341296928327645*E_m*r*s/e
        result.append(E_j)
        return result

    @staticmethod
    def eqn_8_09__s(E_j: float, E_m: float, e: float, r: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        s = 2.93*E_j*e/(E_m*r)
        result.append(s)
        return result

    @staticmethod
    def eqn_8_09__e(E_j: float, E_m: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        e = 0.341296928327645*E_m*r*s/E_j
        result.append(e)
        return result


class SteamJetInjectors:

    @kwasak_static
    def eqn_9_01(rho_s: float = None, w_s: float = None, v: float = None, A: float = None,**kwargs):
        return


    @staticmethod
    def eqn_9_01__rho_s(A: float, v: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        rho_s = w_s/(A*v)
        result.append(rho_s)
        return result

    @staticmethod
    def eqn_9_01__w_s(A: float, rho_s: float, v: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        w_s = A*rho_s*v
        result.append(w_s)
        return result

    @staticmethod
    def eqn_9_01__v(A: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        v = w_s/(A*rho_s)
        result.append(v)
        return result

    @staticmethod
    def eqn_9_01__A(rho_s: float, v: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        A = w_s/(rho_s*v)
        result.append(A)
        return result

    @kwasak_static
    def eqn_9_02(rho_s: float = None, P_m: float = None, w_s: float = None, d_n: float = None,**kwargs):
        return


    @staticmethod
    def eqn_9_02__rho_s(P_m: float, d_n: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        rho_s = 1.334027668054e-6*w_s**2/(P_m*d_n**4)
        result.append(rho_s)
        return result

    @staticmethod
    def eqn_9_02__P_m(d_n: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        P_m = 1.334027668054e-6*w_s**2/(d_n**4*rho_s)
        result.append(P_m)
        return result

    @staticmethod
    def eqn_9_02__w_s(P_m: float, d_n: float, rho_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        w_s = 865.8*d_n**2*sqrt(P_m*rho_s)
        result.append(w_s)
        return result

    @staticmethod
    def eqn_9_02__d_n(P_m: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        d_n = -0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        d_n = 0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        return result

    @kwasak_static
    def eqn_9_03(t_e: float = None, P_s: float = None, w_j: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_9_03__t_e(P_s: float, V: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        t_e = 0.001*V*(2300.0 - 3.0*P_s)/w_j
        result.append(t_e)
        return result

    @staticmethod
    def eqn_9_03__P_s(V: float, t_e: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        P_s = 33.3333333333333*(23.0*V - 10.0*t_e*w_j)/V
        result.append(P_s)
        return result

    @staticmethod
    def eqn_9_03__w_j(P_s: float, V: float, t_e: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        w_j = 0.001*V*(2300.0 - 3.0*P_s)/t_e
        result.append(w_j)
        return result

    @staticmethod
    def eqn_9_03__V(P_s: float, t_e: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        V = -1000.0*t_e*w_j/(3.0*P_s - 2300.0)
        result.append(V)
        return result

    @kwasak_static
    def eqn_9_04(r: float = None, w_s: float = None, AEL: float = None, SC: float = None,**kwargs):
        return


    @staticmethod
    def eqn_9_04__r(AEL: float, SC: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        r = w_s/(AEL*SC)
        result.append(r)
        return result

    @staticmethod
    def eqn_9_04__w_s(AEL: float, SC: float, r: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        w_s = AEL*SC*r
        result.append(w_s)
        return result

    @staticmethod
    def eqn_9_04__AEL(SC: float, r: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        AEL = w_s/(SC*r)
        result.append(AEL)
        return result

    @staticmethod
    def eqn_9_04__SC(AEL: float, r: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        SC = w_s/(AEL*r)
        result.append(SC)
        return result

    @kwasak_static
    def eqn_9_05(V: float = None, w_h: float = None, r_h: float = None, t_h: float = None,**kwargs):
        return


    @staticmethod
    def eqn_9_05__V(r_h: float, t_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        V = t_h*w_h/r_h
        result.append(V)
        return result

    @staticmethod
    def eqn_9_05__w_h(V: float, r_h: float, t_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        w_h = V*r_h/t_h
        result.append(w_h)
        return result

    @staticmethod
    def eqn_9_05__r_h(V: float, t_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        r_h = t_h*w_h/V
        result.append(r_h)
        return result

    @staticmethod
    def eqn_9_05__t_h(V: float, r_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        t_h = V*r_h/w_h
        result.append(t_h)
        return result


class LiquidRing:

    @kwasak_static
    def eqn_10_01(sig_R: float = None, w: float = None, D_r: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_01__sig_R(D_r: float, w: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        sig_R = 0.00436*D_r*w
        result.append(sig_R)
        return result

    @staticmethod
    def eqn_10_01__w(D_r: float, sig_R: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        w = 229.357798165138*sig_R/D_r
        result.append(w)
        return result

    @staticmethod
    def eqn_10_01__D_r(sig_R: float, w: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        D_r = 229.357798165138*sig_R/w
        result.append(D_r)
        return result

    @kwasak_static
    def eqn_10_02(PS: float = None, Q_gas: float = None, dP: float = None, V: float = None, dt: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_02__PS(Q_gas: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        PS = Q_gas - V*dP/dt
        result.append(PS)
        return result

    @staticmethod
    def eqn_10_02__Q_gas(PS: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        Q_gas = PS + V*dP/dt
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_02__dP(PS: float, Q_gas: float, V: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dP = dt*(-PS + Q_gas)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_10_02__V(PS: float, Q_gas: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        V = dt*(-PS + Q_gas)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_10_02__dt(PS: float, Q_gas: float, V: float, dP: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dt = -V*dP/(PS - Q_gas)
        result.append(dt)
        return result

    @kwasak_static
    def eqn_10_03(T: float = None, Q_gas: float = None, N_mfw: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_03__T(N_mfw: float, Q_gas: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        T = 0.108108108108108*Q_gas/N_mfw
        result.append(T)
        return result

    @staticmethod
    def eqn_10_03__Q_gas(N_mfw: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        Q_gas = 9.25*N_mfw*T
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_03__N_mfw(Q_gas: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        N_mfw = 0.108108108108108*Q_gas/T
        result.append(N_mfw)
        return result

    @kwasak_static
    def eqn_10_04(SP_1: float = None, Q_gas: float = None, S_p: float = None, SP_2: float = None, t: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_04__SP_1(Q_gas: float, SP_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_1 = Q_gas + (-Q_gas + SP_2)*exp(S_p*t/V)
        result.append(SP_1)
        return result

    @staticmethod
    def eqn_10_04__Q_gas(SP_1: float, SP_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        Q_gas = -(SP_1 - SP_2*exp(S_p*t/V))/(exp(S_p*t/V) - 1)
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_04__S_p(Q_gas: float, SP_1: float, SP_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        S_p = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_04__SP_2(Q_gas: float, SP_1: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_2 = (Q_gas*exp(S_p*t/V) - Q_gas + SP_1)*exp(-S_p*t/V)
        result.append(SP_2)
        return result

    @staticmethod
    def eqn_10_04__t(Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        t = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/S_p
        result.append(t)
        return result

    @staticmethod
    def eqn_10_04__V(Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        V = S_p*t/log((Q_gas - SP_1)/(Q_gas - SP_2))
        result.append(V)
        return result

    @kwasak_static
    def eqn_10_05(P_2: float = None, S_p: float = None, P_1: float = None, t: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_05__P_2(P_1: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_p*t/V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_05__S_p(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        S_p = V*log(P_1/P_2)/t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_05__P_1(P_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_p*t/V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_05__t(P_1: float, P_2: float, S_p: float, V: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_p
        result.append(t)
        return result

    @staticmethod
    def eqn_10_05__V(P_1: float, P_2: float, S_p: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        V = S_p*t/log(P_1/P_2)
        result.append(V)
        return result

    @kwasak_static
    def eqn_10_06(P_2: float = None, S_a: float = None, P_1: float = None, t: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_06__P_2(P_1: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_a*t/V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_06__S_a(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        S_a = V*log(P_1/P_2)/t
        result.append(S_a)
        return result

    @staticmethod
    def eqn_10_06__P_1(P_2: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_a*t/V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_06__t(P_1: float, P_2: float, S_a: float, V: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_a
        result.append(t)
        return result

    @staticmethod
    def eqn_10_06__V(P_1: float, P_2: float, S_a: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        V = S_a*t/log(P_1/P_2)
        result.append(V)
        return result

    @kwasak_static
    def eqn_10_08(bhp: float = None, delta_h_i: float = None, delta_T: float = None, rho: float = None, c_p: float = None, w_i: float = None, f_a: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_08__bhp(c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        bhp = 0.00315127701375246*c_p*delta_T*f_a*rho - 0.000392927308447937*delta_h_i*w_i
        result.append(bhp)
        return result

    @staticmethod
    def eqn_10_08__delta_h_i(bhp: float, c_p: float, delta_T: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_h_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/w_i
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_10_08__delta_T(bhp: float, c_p: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_T = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*f_a*rho)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_10_08__rho(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        rho = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*f_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_10_08__c_p(bhp: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        c_p = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(delta_T*f_a*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_10_08__w_i(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        w_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/delta_h_i
        result.append(w_i)
        return result

    @staticmethod
    def eqn_10_08__f_a(bhp: float, c_p: float, delta_T: float, delta_h_i: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        f_a = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*rho)
        result.append(f_a)
        return result

    @kwasak_static
    def eqn_10_09(T_c: float = None, delta_T: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_09__T_c(T_s: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_c = T_s + delta_T
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_09__delta_T(T_c: float, T_s: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        delta_T = T_c - T_s
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_10_09__T_s(T_c: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_s = T_c - delta_T
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_10(rho: float = None, bhp_0: float = None, mu: float = None, bhp: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_10__rho(bhp: float, bhp_0: float, mu: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        # [Sympy Failover]
        numerator = (bhp - bhp_0) / bhp_0
        denominator = 0.0155 * mu ** 0.16
        rho_powered = numerator / denominator
        
        # Raise the result to the power of 1/0.84 to get rho
        rho = rho_powered ** (1/0.84)
        
        return [ rho ]


    @staticmethod
    def eqn_10_10__bhp_0(bhp: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp_0 = 2000.0*bhp/(31.0*mu**0.16*rho**0.84 + 1000.0)
        result.append(bhp_0)
        return result

    @staticmethod
    def eqn_10_10__mu(bhp: float, bhp_0: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        mu = -204374584201.104*I*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = 204374584201.104*I*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = -204374584201.104*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = 204374584201.104*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        return result

    @staticmethod
    def eqn_10_10__bhp(bhp_0: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp = 0.0005*bhp_0*(31.0*mu**(4/25)*rho**(21/25) + 1000.0)
        result.append(bhp)
        return result

    @kwasak_static
    def eqn_10_11(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_11__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_c = T_s + 10
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_11__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_s = T_c - 10
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_12(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_12__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_c = T_s + 5
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_12__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_s = T_c - 5
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_13(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_13__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_c = T_s + 25
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_13__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_s = T_c - 25
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_14(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_14__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_c = T_s + 12
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_14__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_s = T_c - 12
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_15(P: float = None, S_p: float = None, p_s: float = None, S_Th: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_15__P(S_Th: float, S_p: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        P = S_Th*p_s/(S_Th - S_p)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_15__S_p(P: float, S_Th: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_p = S_Th*(P - p_s)/P
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_15__p_s(P: float, S_Th: float, S_p: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        p_s = P*(S_Th - S_p)/S_Th
        result.append(p_s)
        return result

    @staticmethod
    def eqn_10_15__S_Th(P: float, S_p: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_Th = P*S_p/(P - p_s)
        result.append(S_Th)
        return result

    @kwasak_static
    def eqn_10_16(S_0: float = None, P: float = None, p_0: float = None, S_Th: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_16__S_0(P: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/(P/(P - p_0))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_16__P(S_0: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        P = p_0*(S_Th/S_0)**(5/3)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_16__p_0(P: float, S_0: float, S_Th: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        p_0 = P - P/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = P - 2.05280095711867*P/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = P - 2.05280095711867*P/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result

    @staticmethod
    def eqn_10_16__S_Th(P: float, S_0: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*(P/(P - p_0))**(3/5)
        result.append(S_Th)
        return result

    @kwasak_static
    def eqn_10_17(S_Th: float = None, P: float = None, S_0: float = None, p_0: float = None, p_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_17__S_Th(P: float, S_0: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*((P - p_s)/(P - p_0))**(3/5)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_17__P(S_0: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        P = (p_0*(S_Th/S_0)**(5/3) - p_s)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = (0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = (0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_17__S_0(P: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/((P - p_s)/(P - p_0))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_17__p_0(P: float, S_0: float, S_Th: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_0 = (P*(S_Th/S_0)**(5/3) - P + p_s)/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = 2.05280095711867*(0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = 2.05280095711867*(0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result

    @staticmethod
    def eqn_10_17__p_s(P: float, S_0: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_s = -P*(S_Th/S_0)**(5/3) + P + p_0*(S_Th/S_0)**(5/3)
        result.append(p_s)
        p_s = -0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 + P + 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        p_s = -0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 + P + 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        return result

    @kwasak_static
    def eqn_10_18(T_i: float = None, S_Th: float = None, P: float = None, S_p: float = None, p_s: float = None, p_c: float = None, T_e: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_18__T_i(P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_i = (-460*P*S_Th + P*S_p*T_e + 460*P*S_p + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*(P - p_s))
        result.append(T_i)
        return result

    @staticmethod
    def eqn_10_18__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_Th = S_p*(P*T_e + 460*P - T_e*p_c - 460*p_c)/(P*T_i + 460*P - T_i*p_s - 460*p_s)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_18__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        P = (S_Th*T_i*p_s + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*T_i + 460*S_Th - S_p*T_e - 460*S_p)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_18__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_p = S_Th*(P*T_i + 460*P - T_i*p_s - 460*p_s)/(P*T_e + 460*P - T_e*p_c - 460*p_c)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_18__p_s(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_s = (P*S_Th*T_i + 460*P*S_Th - P*S_p*T_e - 460*P*S_p + S_p*T_e*p_c + 460*S_p*p_c)/(S_Th*(T_i + 460))
        result.append(p_s)
        return result

    @staticmethod
    def eqn_10_18__p_c(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_c = (-P*S_Th*T_i - 460*P*S_Th + P*S_p*T_e + 460*P*S_p + S_Th*T_i*p_s + 460*S_Th*p_s)/(S_p*(T_e + 460))
        result.append(p_c)
        return result

    @staticmethod
    def eqn_10_18__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_e = (P*S_Th*T_i + 460*P*S_Th - 460*P*S_p - S_Th*T_i*p_s - 460*S_Th*p_s + 460*S_p*p_c)/(S_p*(P - p_c))
        result.append(T_e)
        return result

    @kwasak_static
    def eqn_10_19(T_i: float = None, S_Th: float = None, P: float = None, S_p: float = None, p_s: float = None, p_c: float = None, T_e: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_19__T_i(P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_i = (P*T_e*(S_p/S_Th)**(5/3) + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P - T_e*p_c*(S_p/S_Th)**(5/3) - 460.0*p_c*(S_p/S_Th)**(5/3) + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        return result

    @staticmethod
    def eqn_10_19__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_Th = S_p/((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_19__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        # [Sympy Failover]
        """
        Solves the given equation for P using the provided parameters.
        
        Parameters:
        - S_Th (float): The thermal efficiency of the engine.
        - S_p (float): The power output from the engine.
        - T_e (float): The temperature in Kelvin at which heat is rejected.
        - T_i (float): The temperature in Kelvin at which heat is added.
        - p_c (float): The cost per unit of electricity produced by a cogeneration plant.
        - p_s (float): The selling price per unit of electricity from the engine.
        
        Returns:
        float: The value of P for which the equation holds true.
        """
        # Step 1: Multiply both sides by [(P - p_c)(460 + T_e)] to get rid of denominators
        numerator = (S_p / ((T_i + S_Th * 460)**0ranks.keys()) for ranks in self.hands):
                return hand_rank(sorted_ranks, suit=False)
            else:
                raise ValueError('Invalid input format!')
        
        def rank_hands(self):
            """Return the highest ranked hand."""
            if not self.has_cards():
                return 0
            
            ranks = [rank for card in self.cards for rank, _ in enumerate(['23456789TJQKA', 'shdc'])]
            hands = []
            for ranks_, suits in zip(*map(list, itertools.product(*[ranks]*len(self.hands)))):
                cards_str = ['%s%s' % (rank, suit) for rank, suit in zip(ranks_, suits)]
                hands.append([PokerCard(card) for card in cards_str])
            return max(hands, key=lambda x: hand_rank(sorted(map(get_suit, map(get_rank, x))), ranks=x))
        
        def __repr__(self):
            if self.has_cards():
                cards = ', '.join([str(card) for card in self.cards])
                return '%s(%s)' % (type(self).__name__, cards)
            else:
                return type(self).__name__
            
        def __hash__(self):
            if not self.has_cards():
                raise ValueError('PokerHand must have at least one card to be hashable')
            return hash(tuple([str(card) for card in self.cards]))
        
        @staticmethod
        def parse_hands(*args, **kwargs):
            """Return a list of PokerHand given arguments."""
            if len(args) == 1 and isinstance(args[0], Iterable):
                cards = [card for card in args[0]]
            else:
                cards = []
                for arg in args:
                    if not (isinstance(arg, PokerCard) or isinstance(arg, str)):
                        raise ValueError('Given argument(s) are invalid')
                    elif isinstance(arg, str):
                        rank_suit = arg.strip()
                        if len(rank_suit) > 3:
                            raise ValueError('Input string too long!')
                        
                        rank = get_rank(rank_suit[:-1])
                        suit = get_suit(rank_suit[-1].lower())
                    else:
                        rank, suit = arg.get_card()
                    cards.append(PokerCard(rank, suit))
            return PokerHand(*cards, **kwargs)
        
        """Class for a standard game of poker."""
        
        def __init__(self, players=None):
            self.players = [] if not players else [Player(i+1, 'Player %d' % (i+1)) for i in range(len(players))]
            self._deck = Deck()
            
        @property
        def deck_size(self):
            """Return the number of cards left in the game."""
            return len(self._deck) - 2*self.num_players
        
        @property
            if (occurs % 10 == 0) or (occurrences <= 5):
                raise ValueError('Too few examples of the class "{}"'.format(label))
            
            # If not yet filtered, return all labels.
            if label_idx is None:
                return np.array([l for l in classes_list])
                
            return np.array([classes_list[label_idx]])
        
        def sample(self):
            """Draw a sample from the dataset."""
            
            # Drawing a sample consists of two steps: first, decide how many samples to draw;
            # second, choose which examples will be drawn. The choice is made by sampling uniformly
            # over all examples and then deciding whether each example should be chosen or not based
            # on its class label (if this is a multiclass problem). Note that the number of negative
            # samples for any given positive sample can differ from other classes since we randomly choose
            # which examples to draw in order to create balanced batches.
            
            n_samples = np.random.randint(1, self.num_examples + 1)
            
            label_idx = None
            
            if len(self.classes_) == 2:
                # If it is a binary classification problem, we only care about the class of the first sample
                label_idx = np.random.randint(0, 2)
                
            # We draw 'n_samples' examples uniformly from the dataset and then decide which one to keep based on their labels.
            idx = np.random.choice(self.num_examples, n_samples, replace=False)
            
            samples, targets = [], []
            
            for i in idx:
                sample, target = self[i]
                
                # Check if the example should be kept or discarded depending on its label and problem type
                keep_sample = True
                if len(self.classes_) == 2:
                    if not np.random.choice([True, False], p=[target==label_idx, 1-target==label_idx]):
                        keep_sample = False
                
                # If the sample should be kept we store it in our list of samples and targets
                if keep_sample:
                    samples += [sample]
                    targets += [target]
            
            return np CV, random_state=self.random_state)
        
        def split(self):
            """Return train/test indices to split data in train/test sets."""
            if self._splits is None:
                # We first calculate the number of splits that we can make with our samples dataset, i.e., the maximum possible number 
                # of subsets. Note that each subset should contain at least one sample from every class (otherwise there would be no way to create a multiclass
                # classification problem) and because these classes may have different sizes this constraint needs to be taken into account:
                
                n_samples = self.num_examples
                min_per_class = max([self.num_classes_[i] for i in range(len(self.classes_))])
                
                # Then we compute the maximum number of subsets that can be made by calculating how many samples we have from each class and dividing it with 
                # the minimum per class (since this is a constraint, this gives us an approximation instead of an exact answer). In practice there will likely not be enough samples
                # to split into exactly those numbers of subsets.
                
                max_n_splits = np.floor(np.array([self.num_examples // min_per_class])).astype('int')
                
                while True:
                    # We then calculate how many samples we have from each class and store them as a list of arrays, since this is the output format required by scikit-learn for CV splitting functions.
                    
                    n_samples_per_class = [self.num_examples // min_per_class] * self.num_classes_
                    
                    # We check whether each array in this list has exactly one element (which would mean that the minimum number of samples per class could be respected) 
                    # and if they do we store it, otherwise we increment our maximum n_splits variable by one to indicate that there is an additional sample that needs to be distributed between subsets.
                    
                    if all([len(n_samples_per_class[i]) == 1 for i in range(self.num_classes_)]) and self.num_classes_ == len(np.unique(self.targets)):
                        break
                    else:
                        max_n_splits += 1
                
                # At this point we have computed the maximum number of splits that could be made with our dataset while respecting each class size constraint, and stored it in 'max_n_splits'. We then convert this value to a list so it matches the format required by scikit-learn.
                self._splits = list(range(1, max_n_splits+1))
                
            return super().split()
        
        def _validate_index(self, key):
            """Validate index."""
            
            if isinstance(key, str) and (not self._sample_type == 'numpy'):
                raise NotImplementedError('String indexing currently only supported for numpy arrays.')
            
            return super()._validate_index(key)
        
        """Wrapper around a pandas dataframe."""
        
        def __init__(self, df, random_state=None, data_type='pandas', verbose=False):
            if not isinstance(df, pd.DataFrame):
                raise TypeError('The DataFramePandas constructor requires that the first argument be of type \'pandas.DataFrame\'')
            
            super().__init__(data_type, df.columns.to_numpy(), df.index)
            
            self._random_state = utils.check_random_state(random_state)
            self._df = deepcopy(df)
            if verbose: print('DataFramePandas data loaded')
        
        def __len__(self):
            """Return number of samples in the dataset."""
            
            return len(self.index)
        
        @property
        def num_classes_(self):
            """Return number of classes if target column was provided at initialization, else None."""
            
            # We first check whether a target variable (column) exists in the DataFrame and then count how many unique values are present. If there is only one class we return 1 otherwise we return the total number of unique labels found. Note that if no classes were provided when initializing this object it will default to None as explained in '__init__'.
            try:
                target_column = self._target_columns[0] # The first column is taken to be the target variable since we assume that only one was passed during initialization.
                
                return len(np.unique(self._df[target_column]))
            
            except KeyError:
                return None
        
        def sample(self):
            """Draw a random sample from dataframe."""
            
            i = self._random_state.choice(len(self)) # We draw an index uniformly at random over all rows in the dataset and use it to get its data.
            
            return self._df.iloc[i], None
        
        def apply(self, func):
            """Apply a function element-wise."""
            
            for i in range(len(self)):
                self._df.iloc[i] = func(self._df.iloc[i]) # For each row we call the passed 'func' and assign its output to that same position. This will return a new DataFrame since pandas does not allow changing elements by index in place. Note that if the user wanted to change values directly they could use self._df.iat instead of iloc.
            
            return None
        
        def _validate_index(self, key):
            """Validate index."""
            
            if isinstance(key, str) and (not self._sample_type == 'numpy'):
                raise NotImplementedError('String indexing currently only supported for numpy arrays.')
            
            return super()._validate_index(key)
        
        """Wrapper around a numpy array."""
        
        def __init__(self, data, target=None, dtype='float', index=None, columns=None, verbose=False):
            if not isinstance(data, np.ndarray) and (not is_listy(data)):
                raise TypeError('The Numpy constructor requires that the first argument be of type \'numpy.ndarray\' or a list-like object that can convert to an array automatically')
            
            data = np.asarray(data, dtype=dtype) # We cast our input as a numpy array and store it in self._data attribute since we might need to change its values later on (e.g. after applying transformations). Also note that if the user passed a list-like object without explicitly setting the 'dtype' argument during construction, np.asarray will use this by default which is fine as long as it makes sense for your data format.
            
            self._data_shape = (len(data),) # We store in _data_shape what the shape of our array should be after flattening any multi-dimensional numpy arrays we may have received from user input. If there are no dimensions to remove, this will just give us a one-element tuple with len(self._data) as its only value.
            
            super().__init__(dtype, columns, index) # We also call the initializer of our parent class in order to setup everything needed for indexing (e.g. setting up column names). Note that this will set self._sample_type='numpy' if no target vector is passed during initialization.
            
            self._data = data # And now we store a copy of the original input array as our internal data variable so that it can be changed later on when needed (e.g. after applying transformations).
            
            if verbose: print('Numpy data loaded')
        
        def __len__(self):
            """Return number of samples in dataset."""
            
            return len(self._data) # We just need to find the length of our self._data array.
        
        @property
        def num_classes_(self):
            """Return number of classes if target column was provided at initialization, else None."""
            
            try:
                return len(np.unique(self._df[self._target_columns])) # If there is a target variable we use it to find the unique values in that vector and store their count as self._num_classes_. Otherwise this will just raise an error because there are no classes and thus np.unique will fail since our DataFrame does not have a 'target' column (by default).
            
            except AttributeError: # If we get here it means there is no target variable in the dataset so self._num_classes_ must be None.
                return None
        
        def sample(self):
            """Draw random samples from array."""
            
            i = self._random_state.choice(len(self)) # We draw an index uniformly at random over all rows in the dataset and use it to get its data.
            
            return (self._data[i], None)
        
        def apply(self, func):
            """Apply a function element-wise."""
            
            for i in range(len(self)): # For each row we call the passed 'func' and assign its output to that same position. This will return a new Numpy object since numpy does not allow changing elements by index in place (we could do this with self._data.flat[i] = func(self._data.flat[i]) but it would be much less readable).
                self._data[i] = func(self._data[i]) # If the passed function is applied along one dimension we will only affect that dimension and all other dimensions in our Numpy array will remain unchanged (e.g. a reshape operation will not modify the shape of those additional axes).
            
            return None
        
        def _validate_index(self, key):
            """Validate index."""
            
            if isinstance(key, str) and (not self._sample_type == 'numpy'):
                raise NotImplementedError('String indexing currently only supported for numpy arrays.')
            
            return super()._validate_index(key)
        
        """Wrapper around a pandas DataFrame or Series."""
        
        def __init__(self, data, target=None, columns=None, index=None, dtype='float', verbose=False):
            if isinstance(data, (pd.DataFrame, pd.Series)): # We will initialize our dataset as a Pandas DataFrame or Series object depending on whether we passed one of these types to the constructor or not
                self._df = data.copy() # We store a copy so that it does not change whenever transformations are applied later.
            elif isinstance(data, pd.Series):
                self._df = pd.DataFrame(data) # If we were only passed a Series object then create a DataFrame with this as its sole column using the values from our input data. Note that if you pass in a pandas Series to Pandas constructor it will first convert into numpy array which would make sense since internally all arrays are stored in memory as np.ndarrays, but we don't want to do any such conversion here because sometimes having a DataFrame object might be more convenient (e.g. for handling missing data).
            else: # In the general case where no other types have been passed that would not fit well with our constructor design, just treat it as an array and create a new DataFrame from those values using pd.DataFrame(data) which will also convert any list-like objects (e.g. arrays or tuples) into DataFrames by default (note however that the resulting object might be very inefficient if you try to index it since internally pandas stores data as lists of lists).
                self._df = pd.DataFrame(is_listy(data), columns=columns, index=index) # We create a new empty dataframe with our list-like input passed as values and optionally set the column names and row indices from the arguments if available. This will also convert any array-like objects into DataFrames (including lists of tuples or arrays since these are not considered valid pandas Series inputs).
            
            self._num_samples = len(self._df) # We store in _num_samples what is the length of our dataset (i.e. its number of rows which corresponds to the size of a 2D numpy array after flattening any multi-dimensional structure it may have).
            
            self._dtype = dtype if dtype != 'float' else np.dtype(dtype) # We store in _dtype what type our data will be stored as internally (this is set by default to float64 but can be overridden with the optional input argument here if needed). Note that this only matters when creating new objects from existing ones since internally pandas stores all data as np.ndarrays and then cast them back into the proper dtype when accessing individual entries using integer indices or slices (this does not happen when indexing via label-based methods such as loc[] because those labels are kept intact).
            
            super().__init__(self._df, target, columns, index) # We store in self._target_columns what is the name of our class's target variable if it exists (or None otherwise), and we call the initializer of parent PandasDataset to setup everything needed for indexing. This will also set self._sample_type='pandas' so that users can know whether they should use integer-based or label-based indexing methods with this object.
            
            if verbose: print('Pandas data loaded') # Let the user know their dataset has been created and initialized successfully.
        
        def __len__(self):
            """Return number of samples in dataset."""
            
            return self._num_samples
        
        @property
        def num_classes_(self):
            """Return number of classes if target column was provided at initialization, else None."""
            
            try: # If there is a target variable we can get the unique values from it and count them.
                return len(np.unique(self._df[self._target_column]))
            except KeyError:
                return None
        
        def _validate_index(self, key):
            """Validate index."""
            
            if isinstance(key, str) and (not self._sample_type == 'pandas'): # In case of a string index we want to verify that the given label exists in our dataset.
                try:
            else:
                raise ValueError("Unsupported method for computing pairwise distance")
        def _get_distance_matrix(self, a):
            """Compute the matrix of pairwise distances between samples in `a`."""
            # compute using scipy.spatial.distance.pdist
            X = a
            if not self.pbc:
                # Use euclidean distance (scipy.spatial.distance) or specialized function for fastest computation based on data type and dimensionality
                D = scipy.pdist(X, metric=self._metric)
            else:  # with PBC we can compute more efficiently using the vectorial method in spglib
                n_atoms = X.shape[1]
                symmops = SpacegroupAnalyzer(ase.Atoms(a, pbc=True)).get_symmetry_operations()
                # get matrix of all combinations between symmetry operations
                ops = npoccsd.OCCSRSimilaritySupercellGrouping().get_operation_matrix(symmops)  # TODO: this may be optimized using the group_finder parameter in SpacegroupAnalyzer instead (see https://github.com/schdimian/spglib/issues/27)
                d = ops.dot(X).dot(ops.T)
                n_pairs = n_atoms * (n_atoms - 1) / 2  # for efficient computation we do not use the full pairwise distance matrix, only a vector with all unique distances among atoms in each unit cell
                D = npcsd.pairdistancevec(d[:n_pairs], d[n_pairs:], metric=self._metric)
            return squareform(D), ops
        def _get_nn_matrix(self, a):
            """Compute the matrix of neighbor indices for samples in `a`."""
            if not self.pbc:
                # use scipy.spatial.distance.cdist to efficiently compute all pairwise distance and then find the nearest neighbors (this is very fast)
                D = pdist(a, metric=self._metric)  # note that this is a flat array of distances between each atom in `a` with every other one
                I = npccsd.pynndescent_tree_nearestneighbors().get_indices(D, self.n_neighbors)
            else:  # for periodic systems we can exploit the translational symmetry to compute distances and indices more efficiently using spglib
                ops = self._ops
                D = ops.dot(a).dot(ops.T)
                n_pairs = a.shape[0] * (a.shape[0] - 1) // 2
                I = npccsd.spgsimd_getnnzindices(D[:n_pairs], D[n_pairs:], self._metric, self.n_neighbors).reshape(-1, self.n_neighbors) # TODO: this may be optimized using the group_finder parameter in SpacegroupAnalyzer instead (see https://github.com/schdimian/spglib/issues/27)
            return I
        def fit(self, X):
            """Compute pairwise distances and neighbor indices for samples in `X`.
            
            Parameters
            ----------
            X : array_like of shape [n_samples, n_atoms] or [n_samples, n_frames, n_atoms]
                Array of PBC simulated atomic positions. It can either be a NxM (where M is the number of atoms) 2D array containing only one structure, or an NxMx(P) 3D array containing P frames for each sample. The input format does not change after fitting, it depends on the user's needs before fitting.
            """
            
            if self.n_jobs < 1:
                self.n_jobs = cpu_count() # Use all available cores to compute pairwise distances and neighbor indices at once
            
            X = np.asarray(X, order="C", dtype=float)
            n_samples, *n_dims = X.shape
            if self.pbc:
                self._ops = SpacegroupAnalyzer(ase.Atoms(X)).get_symmetry_operations()  # get operations matrix that map the input data into its periodically repeated version (see https://github. Written by Michiel Coenraads, adapted from PyTorch's implementation
        """A fast 2D Box Filter implemented using convolution."""
        assert data.dim() == 3
        kern = (np.ones([r*2+1, r*2+1], dtype=np.float32) + 1)/(r**2)
        ret = np.zeros_like(data)
        paddings = ((r, r), (r, r), (0, 0))
        data = np.pad(data, paddings, mode='constant')
        data = Variable(torch.FloatTensor(data).unsqueeze(0)).cuda()
        kern = torch.FloatTensor(kern).unsqueeze(0)
        return (Variable(data) * kern).sum().cpu().data[0]  # This is the result, but as a Variable object!
    	""" A fast implementation of a Box Filter for images with batches.
    	Inputs:
    	    x (torch FloatTensor) size(B,C,H,W)
    	Outputs:
    	    y (torch FloatTensor) size(B,C,H,W)"""
    	assert len(x.size()) == 4 # BCHW format
    	r = x.size(3)/2
    	if isinstance(r, int): r=float(r)
    	return boxfilter(x[:,:,:,:],int(np.ceil(r)))
        """Generate a Gaussian kernel with standard deviation sigma and size 2*sigma+1"""
        x = torch.arange(-sigma, sigma+1).to(device)
        grid = torchquadgauss(x, x).to(device)
        return grid/grid.sum()[None]
        """Convolve data with a Gaussian kernel and returns result"""
        assert len(data.shape) == 3 # BCHW format
        g = torch.from_numpy(gaussian_kernel(sigma).reshape(-1)).to(data)
        if isinstance(data, Variable):
            data = data.cpu()
        result = conv2d(data[:,:,:], g[:,None,None]).cuda()
        return torch.clamp(result, 0, 1e6).view_as(data)
        """Computes the output of a Gaussian kernel with input x."""
        # Compute (2 pi)**{-1/2} * e**{-(xx+yy)/2}, size x*y.
        return [ torch.exp(-0.5*(x**2 + y**2)) ]


    @staticmethod
    def eqn_10_19__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_p = S_Th*((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_19__p_s(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_s = (-P*T_e*(S_p/S_Th)**(5/3) + P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P + T_e*p_c*(S_p/S_Th)**(5/3) + 460.0*p_c*(S_p/S_Th)**(5/3))/(T_i + 460.0)
        result.append(p_s)
        p_s = (-0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        p_s = (-0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        return result

    @staticmethod
    def eqn_10_19__p_c(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_c = (P*T_e*(S_p/S_Th)**(5/3) - P*T_i + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P + T_i*p_s + 460.0*p_s)/((S_p/S_Th)**(5/3)*(T_e + 460.0))
        result.append(p_c)
        p_c = 2.05280095711867*(0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        p_c = 2.05280095711867*(0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        return result

    @staticmethod
    def eqn_10_19__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_e = (P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P - T_i*p_s + 460.0*p_c*(S_p/S_Th)**(5/3) - 460.0*p_s)/((S_p/S_Th)**(5/3)*(P - p_c))
        result.append(T_e)
        T_e = 2.05280095711867*(P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        T_e = 2.05280095711867*(P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        return result

    @kwasak_static
    def eqn_10_20(T_i: float = None, P: float = None, S_p: float = None, S_0: float = None, p_0: float = None, p_s: float = None, p_c: float = None, T_e: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_20__T_i(P: float, S_0: float, S_p: float, T_e: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # no closed form solution

    @staticmethod
    def eqn_10_20__P(S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # no closed form solution

    @staticmethod
    def eqn_10_20__S_p(P: float, S_0: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_p = S_0/((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_20__S_0(P: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_0 = S_p*((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_20__p_0(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # no closed form solution

    @staticmethod
    def eqn_10_20__p_s(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        numerator = ((P - p_0)*(460 + T_i)) * (P - p_c)
        denominator = S_0**(1/0.6) * S_p * (460 + T_e)
        
        return [ P - (numerator / denominator) ]


    @staticmethod
    def eqn_10_20__p_c(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        numerator = (S_0 / S_p) ** (2/3)
        denominator = (P * (P - p_s)) * (460 + T_e)
        p_c = P - numerator * denominator
        return [ p_c ]


    @staticmethod
    def eqn_10_20__T_e(P: float, S_0: float, S_p: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # no closed form solution

    @kwasak_static
    def eqn_10_21(P: float = None, P_prime: float = None, P_d: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_21__P(P_d: float, P_prime: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P = P_d*P_prime/760
        result.append(P)
        return result

    @staticmethod
    def eqn_10_21__P_prime(P: float, P_d: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_prime = 760*P/P_d
        result.append(P_prime)
        return result

    @staticmethod
    def eqn_10_21__P_d(P: float, P_prime: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_d = 760*P/P_prime
        result.append(P_d)
        return result


class RotaryPistonVane:

    @kwasak_static
    def eqn_11_01(PS: float = None, Q_0: float = None, Q_external_gas_throughput: float = None, dT: float = None, dP: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_01__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_01__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_01__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_01__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

    @staticmethod
    def eqn_11_01__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_01__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @kwasak_static
    def eqn_11_02(Q_0: float = None, Q_external_gas_throughput: float = None, S_vol_pump_speed: float = None, SP_2: float = None, t: float = None, SP_1: float = None, Q: float = None, V: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_02__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        exp_term = np.exp((t / (V / S_vol_pump_speed)))
        right_side = SP_1 - Q_external_gas_throughput - (SP_2 - Q) * exp_term
        
        Q0 = SP_1 - Q_external_gas_throughput - right_side
        
        return [ Q0 ]


    @staticmethod
    def eqn_11_02__Q_external_gas_throughput(Q: float, Q_0: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        """
        Calculate the external gas throughput (Q_ext_gas_thrput) using an algebraic equation.
        
        Parameters:
        Q (float): The current value of a variable to solve for.
        Q_0 (float): A constant term in the equation related to previous gas throuput.
        SP_1 (float): Partial pressure at state point 1.
        SP_2 (float): Partial pressure at state point 2.
        S_vol_pump_speed (float): Pump speed volume per unit time.
        V (float): Total system volume.
        t (float): Time variable in the equation.
        
        Returns:
        float: The calculated value of external gas throughput (Q_ext_gas_thrput).
        """
        numerator = SP_1 - Q_0
        denominator = SP_2 - (Q + Q_0)
        
        # Exponentiate the term t * S_vol_pump_speed / V to remove the natural logarithm.
        exponentiation_term = exp(t * S_vol_pump_speed / V)
        
        # Calculate Q_external_gas_throughput using the derived equation.
        Q_ext_gas_thrput = numerator * exponentiation_term - denominator
        
        return [ Q_ext_gas_thrput ]


    @staticmethod
    def eqn_11_02__S_vol_pump_speed(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        """
        Calculate the volumetric pump speed given an equation.
        
        Parameters:
        - Q (float): The gas throughput in terms of a variable quantity measure.
        - Q_0 (float): Base quantity measurement for external gas throughput.
        - Q_external_gas_throughput (float): External gas flow rate.
        - SP_1 (float): Point pressure 1 (constant).
        - SP_2 (float): Point pressure 2 (constant).
        - V (float): Total volume in terms of a constant measurement.
        - t (float): Time variable for the equation.
        
        Returns:
        - S_vol_pump_speed (float): The isolated volumetric pump speed.
        """
        
        # Calculate the natural logarithm term using given pressures and throughputs
        ln_term = log((SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0)))
        
        # Isolate S_vol_pump_speed by rearranging the equation
        S_vol_pump_speed = V / (t * ln_term)
        
        return [ S_vol_pump_speed ]


    @staticmethod
    def eqn_11_02__SP_2(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        
        import math
        
        SP_2 = (exp(t * (S_vol_pump_speed / V))) * (SP_1 - (Q_external_gas_throughput + Q_0) + exp(t * (S_vol_pump_speed / V)) * (Q + Q_0))
        
        return [ SP_2 ]


    @staticmethod
    def eqn_11_02__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        """
        Solve the equation for t using given parameters in terms of Q.
        
        Parameters:
        Q (float): The throughput of gas flowing into the system (in appropriate units).
        Q_0 (float): The base flow rate at equilibrium or initial condition.
        Q_external_gas_throughput (float): The external gas throughput in the system.
        SP_1 (float): The first standard pressure point in the system.
        SP_2 (float): The second standard pressure point in the system.
        S_vol_pump_speed (float): The speed of the volume pump (in appropriate units).
        V (float): The total volume at which gas is being processed or analyzed.
        
        Returns:
        float: The calculated time 't' based on the given equation parameters.
        """
        # Calculate the natural logarithm part of the equation, ensuring that denominator is not zero to avoid math domain error
        ln_ratio = (SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0))
        
        if SP_2 - (Q + Q_0) == 0:
            raise ValueError("The denominator in the logarithm cannot be zero.")
            
        t = V / S_vol_pump_speed * log(ln_ratio)
        
        return [ t ]


    @staticmethod
    def eqn_11_02__SP_1(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        
        SP_1 = exp((t * S_vol_pump_speed) / V) * (SP_2 - (Q + Q_0)) + (Q_external_gas_throughput + Q_0)
        
        return [ SP_1 ]


    @staticmethod
    def eqn_11_02__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        result = SP_2 - SP_1 - Q_external_gas_throughput / (S_vol_pump_speed/V * exp(t)) + 1
        return [ result ]


    @staticmethod
    def eqn_11_02__V(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        numerator = 1 / (S_vol_pump_speed * log((SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0))))
        return [ numerator * exp(-t) ]


    @kwasak_static
    def eqn_11_03(F_s: float = None, t: float = None, t_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_03__F_s(t: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        F_s = t/t_c
        result.append(F_s)
        return result

    @staticmethod
    def eqn_11_03__t(F_s: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t = F_s*t_c
        result.append(t)
        return result

    @staticmethod
    def eqn_11_03__t_c(F_s: float, t: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t_c = t/F_s
        result.append(t_c)
        return result

    @kwasak_static
    def eqn_11_04(p_g: float = None, p_s: float = None, p_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_04__p_g(p_s: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_g = p_s - p_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_04__p_s(p_g: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_s = p_g + p_v
        result.append(p_s)
        return result

    @staticmethod
    def eqn_11_04__p_v(p_g: float, p_s: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_v = 0
        result.append(p_v)
        p_v = -p_g + p_s
        result.append(p_v)
        return result

    @kwasak_static
    def eqn_11_05(P_0_v: float = None, P_D: float = None, p_v_max: float = None, p_g: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_05__P_0_v(P_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D*p_v_max/(p_g + p_v_max)
        result.append(P_0_v)
        return result

    @staticmethod
    def eqn_11_05__P_D(P_0_v: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v*(p_g + p_v_max)/p_v_max
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_05__p_v_max(P_0_v: float, P_D: float, p_g: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v*p_g/(P_0_v - P_D)
        result.append(p_v_max)
        return result

    @staticmethod
    def eqn_11_05__p_g(P_0_v: float, P_D: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max*(-P_0_v + P_D)/P_0_v
        result.append(p_g)
        return result

    @kwasak_static
    def eqn_11_06(S_D: float = None, p_b: float = None, S_B: float = None, P_0_V: float = None, p_v_max: float = None, P_D: float = None, p_g: float = None, P_v_0: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_06__S_D(P_0_V: float, P_D: float, P_v_0: float, S_B: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = P_D*S_B*(P_0_V - p_b)/(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)
        result.append(S_D)
        return result

    @staticmethod
    def eqn_11_06__p_b(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (P_0_V*P_D*S_B - P_D*S_D*p_v_max + P_v_0*S_D*p_g + P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(p_b)
        return result

    @staticmethod
    def eqn_11_06__S_B(P_0_V: float, P_D: float, P_v_0: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = S_D*(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)/(P_D*(P_0_V - p_b))
        result.append(S_B)
        return result

    @staticmethod
    def eqn_11_06__P_0_V(P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_g - P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(P_0_V)
        return result

    @staticmethod
    def eqn_11_06__p_v_max(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V*P_D*S_B - P_D*S_B*p_b + P_v_0*S_D*p_g)/(S_D*(P_D - P_v_0))
        result.append(p_v_max)
        return result

    @staticmethod
    def eqn_11_06__P_D(P_0_V: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0*S_D*(p_g + p_v_max)/(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_06__p_g(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (-P_0_V*P_D*S_B + P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_v_max)/(P_v_0*S_D)
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_06__P_v_0(P_0_V: float, P_D: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = P_D*(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)/(S_D*(p_g + p_v_max))
        result.append(P_v_0)
        return result
""
import tru
y = {}
for u,o in enumerate(filter(lambda o:str(o)[0].isalpha() and str(o)[0].capitalize()==str(o)[0] and str(o) not in map(lambda a:a.strip(),'I, Piecewise, LambertW, Eq, symbols'.split(',')),dir())):
    print(f'@@@{u+1}.',o, type(o))
    # try:
    truth = False
    for tempt in range(budget:=5):
        try:
            truth = truth or tru.Verify(vars()[o]).verify() 
        except ValueError as ve:
            if (m:="math domain error") in str(ve):pass
            # elif(m:=)
            print("[ERROR]"+":"*99,m)
            # print(str(ve));1/0
    print("+"*8*8,*((truth,) if (b:=isinstance(truth,bool)) else (truth.items())),sep=('\n\t'if not b else ''))
    y[o] = truth
print(*[yo for yo in y.items()],sep=('\n\t'))

    