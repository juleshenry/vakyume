from math import log, sqrt, exp
from sympy import I, Piecewise, LambertW, Eq


class VacuumTheory:

    @staticmethod
    def eqn_1_3__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T

    @staticmethod
    def eqn_1_3__v(T: float, k: float, m: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T

    @staticmethod
    def eqn_1_3__m(T: float, k: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T

    @staticmethod
    def eqn_1_3__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T

    @staticmethod
    def eqn_1_7__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T

    @staticmethod
    def eqn_1_7__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T

    @staticmethod
    def eqn_1_7__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T

    @staticmethod
    def eqn_1_7__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T

    @staticmethod
    def eqn_1_7__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T

    @staticmethod
    def eqn_1_8__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T

    @staticmethod
    def eqn_1_8__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T

    @staticmethod
    def eqn_1_8__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T

    @staticmethod
    def eqn_1_8__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T

    @staticmethod
    def eqn_1_8__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T

    @staticmethod
    def eqn_1_8__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T

    @staticmethod
    def eqn_1_9__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)

    @staticmethod
    def eqn_1_9__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)

    @staticmethod
    def eqn_1_9__T(M: float, P: float, R: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)

    @staticmethod
    def eqn_1_9__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)

    @staticmethod
    def eqn_1_9__rho(M: float, P: float, R: float, T: float):
        # [.pyeqn] rho = P * M / (R * T)

    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2

    @staticmethod
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2

    @staticmethod
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2

    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2

    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2

    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)

    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)

    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)

    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)

    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)

    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
        # [.pyeqn] Total_P = sum_partial_pressures

    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
        # [.pyeqn] Total_P = sum_partial_pressures

    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
        # [.pyeqn] y_a = n_a / n

    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
        # [.pyeqn] y_a = n_a / n

    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
        # [.pyeqn] y_a = n_a / n

    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P

    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
        # [.pyeqn] y_a = p_a / P

    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P


class FluidFlowVacuumLines:

    @staticmethod
    def eqn_2_1__Re(D: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu

    @staticmethod
    def eqn_2_1__D(Re: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu

    @staticmethod
    def eqn_2_1__v(D: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = rho * D * v / mu

    @staticmethod
    def eqn_2_1__mu(D: float, Re: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu

    @staticmethod
    def eqn_2_1__rho(D: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu

    @staticmethod
    def eqn_2_2__psi(delta: float, lambd: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5

    @staticmethod
    def eqn_2_2__lambd(delta: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5

    @staticmethod
    def eqn_2_2__delta(lambd: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5

    @staticmethod
    def eqn_2_3__kn(D: float, lambd: float):
        # [.pyeqn] kn = lambd / D

    @staticmethod
    def eqn_2_3__lambd(D: float, kn: float):
        # [.pyeqn] kn = lambd / D

    @staticmethod
    def eqn_2_3__D(kn: float, lambd: float):
        # [.pyeqn] kn = lambd / D

    @staticmethod
    def eqn_2_4__vel_grad(_beta: float, mu: float):
        # [.pyeqn] _beta = mu * vel_grad

    @staticmethod
    def eqn_2_4___beta(mu: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad

    @staticmethod
    def eqn_2_4__mu(_beta: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad

    @staticmethod
    def eqn_2_5__q(D: float, L: float, delta_P: float, mu: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)

    @staticmethod
    def eqn_2_5__delta_P(D: float, L: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)

    @staticmethod
    def eqn_2_5__D(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)

    @staticmethod
    def eqn_2_5__mu(D: float, L: float, delta_P: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)

    @staticmethod
    def eqn_2_5__L(D: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)

    @staticmethod
    def eqn_2_6__v_a(lambd: float, mu: float, rho: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a

    @staticmethod
    def eqn_2_6__lambd(mu: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a

    @staticmethod
    def eqn_2_6__mu(lambd: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a

    @staticmethod
    def eqn_2_6__rho(lambd: float, mu: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a

    @staticmethod
    def eqn_2_7__k(T: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5

    @staticmethod
    def eqn_2_7__v_a(T: float, k: float, m: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5

    @staticmethod
    def eqn_2_7__m(T: float, k: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5

    @staticmethod
    def eqn_2_7__T(k: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5

    @staticmethod
    def eqn_2_8__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_8__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_8__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_8__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_8__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_8__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_8__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_8__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)

    @staticmethod
    def eqn_2_10__delta_P(Suc_Pres: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P

    @staticmethod
    def eqn_2_10__Suc_Pres(delta_P: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P

    @staticmethod
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P

    @staticmethod
    def eqn_2_11__g_c(D: float, L: float, f: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)

    @staticmethod
    def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)

    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)

    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)

    @staticmethod
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)

    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)

    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)

    @staticmethod
    def eqn_2_12__delta_P(L: float, d: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)

    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)

    @staticmethod
    def eqn_2_12__v(L: float, d: float, delta_P: float, f: float, g: float, rho: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)

    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)

    @staticmethod
    def eqn_2_12__f(L: float, d: float, delta_P: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)

    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)

    @staticmethod
    def eqn_2_13__q(L: float, d: float, delta_P: float, f: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)

    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)

    @staticmethod
    def eqn_2_13__delta_P(L: float, d: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)

    @staticmethod
    def eqn_2_13__d(L: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)

    @staticmethod
    def eqn_2_13__f(L: float, d: float, delta_P: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)

    @staticmethod
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)

    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5

    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5

    @staticmethod
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5

    @staticmethod
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5

    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5

    @staticmethod
    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5

    @staticmethod
    def eqn_2_15__f(Re: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)

    @staticmethod
    def eqn_2_15__Re(f: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2

    @staticmethod
    def eqn_2_17__v(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2

    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2

    @staticmethod
    def eqn_2_17__q(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4

    @staticmethod
    def eqn_2_17__delta_P(L: float, d: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4

    @staticmethod
    def eqn_2_17__d(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4

    @staticmethod
    def eqn_2_17__mu(L: float, d: float, delta_P: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4

    @staticmethod
    def eqn_2_17__L(d: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4


class PressMgmt:

    @staticmethod
    def eqn_3_1__BarometricPressure(Abs_Pressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum

    @staticmethod
    def eqn_3_1__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum

    @staticmethod
    def eqn_3_1__Abs_Pressure(BarometricPressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum

    @staticmethod
    def eqn_3_2__G(G_C: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)

    @staticmethod
    def eqn_3_2__G_C(G: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)

    @staticmethod
    def eqn_3_2__H(G: float, G_C: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)

    @staticmethod
    def eqn_3_2__P(G: float, G_C: float, H: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)

    @staticmethod
    def eqn_3_2__rho(G: float, G_C: float, H: float, P: float):
        # [.pyeqn] P = G / (G_C * rho * H)

    @staticmethod
    def eqn_3_3__P_P(H_1: float, H_2: float, P: float):
        # [.pyeqn] P_P - P = H_2 - H_1

    @staticmethod
    def eqn_3_3__H_2(H_1: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1

    @staticmethod
    def eqn_3_3__P(H_1: float, H_2: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1

    @staticmethod
    def eqn_3_3__H_1(H_2: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1

    @staticmethod
    def eqn_3_4__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA

    @staticmethod
    def eqn_3_4__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA

    @staticmethod
    def eqn_3_4__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA

    @staticmethod
    def eqn_3_5__P_P(P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)

    @staticmethod
    def eqn_3_5__V_P(P: float, P_P: float, V: float):
        # [.pyeqn] P_P = P * (V / V_P)

    @staticmethod
    def eqn_3_5__P(P_P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)

    @staticmethod
    def eqn_3_5__V(P: float, P_P: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)

    @staticmethod
    def eqn_3_6__V_P(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)

    @staticmethod
    def eqn_3_6__P(H_1: float, H_2: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)

    @staticmethod
    def eqn_3_6__H_1(H_2: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)

    @staticmethod
    def eqn_3_6__H_2(H_1: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)

    @staticmethod
    def eqn_3_6__V(H_1: float, H_2: float, P: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
