from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak_static
from vakyume.config import UnsolvedException
import numpy as np

class AirLeak:
    @kwasak_static
    def eqn_4_10(T=None, V=None, del_P=None, leakage=None, t=None, **kwargs):
        return

    def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float, **kwargs):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        T = 3.127*V*del_P/(leakage*t)
        result.append(T)
        return result
    def eqn_4_10__V(T: float, del_P: float, leakage: float, t: float, **kwargs):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        V = 0.319795330988168*T*leakage*t/del_P
        result.append(V)
        return result
    def eqn_4_10__del_P(T: float, V: float, leakage: float, t: float, **kwargs):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        del_P = 0.319795330988168*T*leakage*t/V
        result.append(del_P)
        return result
    def eqn_4_10__leakage(T: float, V: float, del_P: float, t: float, **kwargs):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        leakage = 3.127*V*del_P/(T*t)
        result.append(leakage)
        return result
    def eqn_4_10__t(T: float, V: float, del_P: float, leakage: float, **kwargs):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        t = 3.127*V*del_P/(T*leakage)
        result.append(t)
        return result
    @kwasak_static
    def eqn_4_7(W=None, W_T=None, sum_individual_leak_rates=None, **kwargs):
        return

    def eqn_4_7__W_T(W: float, sum_individual_leak_rates: float, **kwargs):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return result
    def eqn_4_7__W(W_T: float, sum_individual_leak_rates: float, **kwargs):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return result
    def eqn_4_7__sum_individual_leak_rates(W: float, W_T: float, **kwargs):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result
class FluidFlowVacuumLines:
    @kwasak_static
    def eqn_2_1(D=None, Re=None, mu=None, rho=None, v=None, **kwargs):
        return

    def eqn_2_1__D(Re: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return result
    def eqn_2_1__Re(D: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return result
    def eqn_2_1__mu(D: float, Re: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return result
    def eqn_2_1__rho(D: float, Re: float, mu: float, v: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return result
    def eqn_2_1__v(D: float, Re: float, mu: float, rho: float, **kwargs):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return result
    @kwasak_static
    def eqn_2_10(Suc_Pres=None, delta_P=None, oper_press=None, **kwargs):
        return

    def eqn_2_10__Suc_Pres(delta_P: float, oper_press: float, **kwargs):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        Suc_Pres = -delta_P + oper_press
        result.append(Suc_Pres)
        return result
    def eqn_2_10__delta_P(Suc_Pres: float, oper_press: float, **kwargs):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        delta_P = -Suc_Pres + oper_press
        result.append(delta_P)
        return result
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float, **kwargs):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return result
    @kwasak_static
    def eqn_2_14(M=None, R=None, T=None, g_c=None, k=None, v_s=None, **kwargs):
        return

    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return result
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return result
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return result
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return result
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return result
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float, **kwargs):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return result
    @kwasak_static
    def eqn_2_15(Re=None, f=None, **kwargs):
        return

    def eqn_2_15__Re(f: float, **kwargs):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736/f**4
        result.append(Re)
        return result
    def eqn_2_15__f(Re: float, **kwargs):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316/Re**(1/4)
        result.append(f)
        return result
    @kwasak_static
    def eqn_2_16(Re=None, f=None, **kwargs):
        return

    def eqn_2_16__Re(f: float, **kwargs):
        # [.pyeqn] f = 64 / Re
        result = []
        Re = 64/f
        result.append(Re)
        return result
    def eqn_2_16__f(Re: float, **kwargs):
        # [.pyeqn] f = 64 / Re
        result = []
        f = 64/Re
        result.append(f)
        return result
    @kwasak_static
    def eqn_2_18a(D_eq=None, R_ll=None, **kwargs):
        return

    def eqn_2_18a__D_eq(R_ll: float, **kwargs):
        # [.pyeqn] D_eq = 4 * R_ll
        result = []
        D_eq = 4*R_ll
        result.append(D_eq)
        return result
    def eqn_2_18a__R_ll(D_eq: float, **kwargs):
        # [.pyeqn] D_eq = 4 * R_ll
        result = []
        R_ll = D_eq/4
        result.append(R_ll)
        return result
    @kwasak_static
    def eqn_2_18b(R_ll=None, h=None, w=None, **kwargs):
        return

    def eqn_2_18b__R_ll(h: float, w: float, **kwargs):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        R_ll = h*w/(2*(h + w))
        result.append(R_ll)
        return result
    def eqn_2_18b__h(R_ll: float, w: float, **kwargs):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        h = 2*R_ll*w/(-2*R_ll + w)
        result.append(h)
        return result
    def eqn_2_18b__w(R_ll: float, h: float, **kwargs):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        w = 2*R_ll*h/(-2*R_ll + h)
        result.append(w)
        return result
    @kwasak_static
    def eqn_2_19a(R_ll=None, Re=None, mu=None, rho=None, v=None, **kwargs):
        return

    def eqn_2_19a__R_ll(Re: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        R_ll = Re*mu/(4*rho*v)
        result.append(R_ll)
        return result
    def eqn_2_19a__Re(R_ll: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        Re = 4*R_ll*rho*v/mu
        result.append(Re)
        return result
    def eqn_2_19a__mu(R_ll: float, Re: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        mu = 4*R_ll*rho*v/Re
        result.append(mu)
        return result
    def eqn_2_19a__rho(R_ll: float, Re: float, mu: float, v: float, **kwargs):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        rho = Re*mu/(4*R_ll*v)
        result.append(rho)
        return result
    def eqn_2_19a__v(R_ll: float, Re: float, mu: float, rho: float, **kwargs):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        v = Re*mu/(4*R_ll*rho)
        result.append(v)
        return result
    @kwasak_static
    def eqn_2_19b(Re=None, h=None, mu=None, rho=None, v=None, w=None, **kwargs):
        return

    def eqn_2_19b__Re(h: float, mu: float, rho: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        Re = 2*h*rho*v*w/(mu*(h + w))
        result.append(Re)
        return result
    def eqn_2_19b__h(Re: float, mu: float, rho: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        h = Re*mu*w/(-Re*mu + 2*rho*v*w)
        result.append(h)
        return result
    def eqn_2_19b__mu(Re: float, h: float, rho: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        mu = 2*h*rho*v*w/(Re*(h + w))
        result.append(mu)
        return result
    def eqn_2_19b__rho(Re: float, h: float, mu: float, v: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        rho = Re*mu*(h + w)/(2*h*v*w)
        result.append(rho)
        return result
    def eqn_2_19b__v(Re: float, h: float, mu: float, rho: float, w: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        v = Re*mu*(h + w)/(2*h*rho*w)
        result.append(v)
        return result
    def eqn_2_19b__w(Re: float, h: float, mu: float, rho: float, v: float, **kwargs):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        w = Re*h*mu/(-Re*mu + 2*h*rho*v)
        result.append(w)
        return result
    @kwasak_static
    def eqn_2_20(L=None, sum_equivalent_length=None, sum_pipe=None, **kwargs):
        return

    def eqn_2_20__L(sum_equivalent_length: float, sum_pipe: float, **kwargs):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        L = sum_equivalent_length + sum_pipe
        result.append(L)
        return result
    def eqn_2_20__sum_equivalent_length(L: float, sum_pipe: float, **kwargs):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        sum_equivalent_length = L - sum_pipe
        result.append(sum_equivalent_length)
        return result
    def eqn_2_20__sum_pipe(L: float, sum_equivalent_length: float, **kwargs):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        sum_pipe = L - sum_equivalent_length
        result.append(sum_pipe)
        return result
    @kwasak_static
    def eqn_2_22(P_s=None, Q_throughput=None, S_p=None, **kwargs):
        return

    def eqn_2_22__P_s(Q_throughput: float, S_p: float, **kwargs):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        P_s = Q_throughput/S_p
        result.append(P_s)
        return result
    def eqn_2_22__Q_throughput(P_s: float, S_p: float, **kwargs):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        Q_throughput = P_s*S_p
        result.append(Q_throughput)
        return result
    def eqn_2_22__S_p(P_s: float, Q_throughput: float, **kwargs):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        S_p = Q_throughput/P_s
        result.append(S_p)
        return result
    @kwasak_static
    def eqn_2_25(C=None, P_1=None, P_2=None, Q_throughput=None, **kwargs):
        return

    def eqn_2_25__C(P_1: float, P_2: float, Q_throughput: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        C = Q_throughput/(P_1 - P_2)
        result.append(C)
        return result
    def eqn_2_25__P_1(C: float, P_2: float, Q_throughput: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        P_1 = P_2 + Q_throughput/C
        result.append(P_1)
        return result
    def eqn_2_25__P_2(C: float, P_1: float, Q_throughput: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        P_2 = P_1 - Q_throughput/C
        result.append(P_2)
        return result
    def eqn_2_25__Q_throughput(C: float, P_1: float, P_2: float, **kwargs):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        Q_throughput = C*(P_1 - P_2)
        result.append(Q_throughput)
        return result
    @kwasak_static
    def eqn_2_29(C=None, S_1=None, S_2=None, **kwargs):
        return

    def eqn_2_29__C(S_1: float, S_2: float, **kwargs):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        C = -S_1*S_2/(S_1 - S_2)
        result.append(C)
        return result
    def eqn_2_29__S_1(C: float, S_2: float, **kwargs):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_1 = C*S_2/(C + S_2)
        result.append(S_1)
        return result
    def eqn_2_29__S_2(C: float, S_1: float, **kwargs):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_2 = C*S_1/(C - S_1)
        result.append(S_2)
        return result
    @kwasak_static
    def eqn_2_3(D=None, kn=None, lambd=None, **kwargs):
        return

    def eqn_2_3__D(kn: float, lambd: float, **kwargs):
        # [.pyeqn] kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return result
    def eqn_2_3__kn(D: float, lambd: float, **kwargs):
        # [.pyeqn] kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return result
    def eqn_2_3__lambd(D: float, kn: float, **kwargs):
        # [.pyeqn] kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return result
    @kwasak_static
    def eqn_2_31(C=None, S_p=None, S_pump_speed=None, **kwargs):
        return

    def eqn_2_31__C(S_p: float, S_pump_speed: float, **kwargs):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        C = S_p*S_pump_speed/(S_p - S_pump_speed)
        result.append(C)
        return result
    def eqn_2_31__S_p(C: float, S_pump_speed: float, **kwargs):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_p = C*S_pump_speed/(C - S_pump_speed)
        result.append(S_p)
        return result
    def eqn_2_31__S_pump_speed(C: float, S_p: float, **kwargs):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_pump_speed = C*S_p/(C + S_p)
        result.append(S_pump_speed)
        return result
    @kwasak_static
    def eqn_2_32(C_series=None, geometric_sum_C=None, **kwargs):
        return

    def eqn_2_32__C_series(geometric_sum_C: float, **kwargs):
        # [.pyeqn] 1 / C_series = geometric_sum_C
        result = []
        C_series = 1/geometric_sum_C
        result.append(C_series)
        return result
    def eqn_2_32__geometric_sum_C(C_series: float, **kwargs):
        # [.pyeqn] 1 / C_series = geometric_sum_C
        result = []
        geometric_sum_C = 1/C_series
        result.append(geometric_sum_C)
        return result
    @kwasak_static
    def eqn_2_33(C_paralell=None, arithmetic_sum_C=None, **kwargs):
        return

    def eqn_2_33__C_paralell(arithmetic_sum_C: float, **kwargs):
        # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
        result = []
        C_paralell = 1/arithmetic_sum_C
        result.append(C_paralell)
        return result
    def eqn_2_33__arithmetic_sum_C(C_paralell: float, **kwargs):
        # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
        result = []
        arithmetic_sum_C = 1/C_paralell
        result.append(arithmetic_sum_C)
        return result
    @kwasak_static
    def eqn_2_36(C=None, C_0=None, F_t=None, **kwargs):
        return

    def eqn_2_36__C_0(C: float, F_t: float, **kwargs):
        # [.pyeqn] C = C_0 * F_t
        result = []
        C_0 = C/F_t
        result.append(C_0)
        return result
    def eqn_2_36__C(C_0: float, F_t: float, **kwargs):
        # [.pyeqn] C = C_0 * F_t
        result = []
        C = C_0*F_t
        result.append(C)
        return result
    def eqn_2_36__F_t(C: float, C_0: float, **kwargs):
        # [.pyeqn] C = C_0 * F_t
        result = []
        F_t = C/C_0
        result.append(F_t)
        return result
    @kwasak_static
    def eqn_2_37(A=None, C=None, F_t=None, M=None, T=None, **kwargs):
        return

    def eqn_2_37__A(C: float, F_t: float, M: float, T: float, **kwargs):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        A = 0.000681714375311032*C**2*M/(F_t*T)
        result.append(A)
        return result
    def eqn_2_37__C(A: float, F_t: float, M: float, T: float, **kwargs):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        C = 38.3*sqrt(A*F_t*T/M)
        result.append(C)
        return result
    def eqn_2_37__F_t(A: float, C: float, M: float, T: float, **kwargs):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        F_t = 0.000681714375311032*C**2*M/(A*T)
        result.append(F_t)
        return result
    def eqn_2_37__M(A: float, C: float, F_t: float, T: float, **kwargs):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        M = 1466.89*A*F_t*T/C**2
        result.append(M)
        return result
    def eqn_2_37__T(A: float, C: float, F_t: float, M: float, **kwargs):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        T = 0.000681714375311032*C**2*M/(A*F_t)
        result.append(T)
        return result
    @kwasak_static
    def eqn_2_4(_beta=None, mu=None, vel_grad=None, **kwargs):
        return

    def eqn_2_4___beta(mu: float, vel_grad: float, **kwargs):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return result
    def eqn_2_4__mu(_beta: float, vel_grad: float, **kwargs):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return result
    def eqn_2_4__vel_grad(_beta: float, mu: float, **kwargs):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return result
    @kwasak_static
    def eqn_2_6(lambd=None, mu=None, rho=None, v_a=None, **kwargs):
        return

    def eqn_2_6__lambd(mu: float, rho: float, v_a: float, **kwargs):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return result
    def eqn_2_6__mu(lambd: float, rho: float, v_a: float, **kwargs):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return result
    def eqn_2_6__rho(lambd: float, mu: float, v_a: float, **kwargs):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return result
    def eqn_2_6__v_a(lambd: float, mu: float, rho: float, **kwargs):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return result
    @kwasak_static
    def eqn_2_7(T=None, k=None, m=None, v_a=None, **kwargs):
        return

    def eqn_2_7__T(k: float, m: float, v_a: float, **kwargs):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        T = 0.392699081698724*m*v_a**2/k
        result.append(T)
        return result
    def eqn_2_7__k(T: float, m: float, v_a: float, **kwargs):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        k = 0.392699081698724*m*v_a**2/T
        result.append(k)
        return result
    def eqn_2_7__m(T: float, k: float, v_a: float, **kwargs):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return result
    def eqn_2_7__v_a(T: float, k: float, m: float, **kwargs):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return result
class LiquidRing:
    @kwasak_static
    def eqn_10_1(D_r=None, sig_R=None, w=None, **kwargs):
        return

    def eqn_10_1__D_r(sig_R: float, w: float, **kwargs):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        D_r = 229.357798165138*sig_R/w
        result.append(D_r)
        return result
    def eqn_10_1__sig_R(D_r: float, w: float, **kwargs):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        sig_R = 0.00436*D_r*w
        result.append(sig_R)
        return result
    def eqn_10_1__w(D_r: float, sig_R: float, **kwargs):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        w = 229.357798165138*sig_R/D_r
        result.append(w)
        return result
    @kwasak_static
    def eqn_10_11(T_c=None, T_s=None, **kwargs):
        return

    def eqn_10_11__T_c(T_s: float, **kwargs):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_c = T_s + 10
        result.append(T_c)
        return result
    def eqn_10_11__T_s(T_c: float, **kwargs):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_s = T_c - 10
        result.append(T_s)
        return result
    @kwasak_static
    def eqn_10_12(T_c=None, T_s=None, **kwargs):
        return

    def eqn_10_12__T_c(T_s: float, **kwargs):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_c = T_s + 5
        result.append(T_c)
        return result
    def eqn_10_12__T_s(T_c: float, **kwargs):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_s = T_c - 5
        result.append(T_s)
        return result
    @kwasak_static
    def eqn_10_13(T_c=None, T_s=None, **kwargs):
        return

    def eqn_10_13__T_c(T_s: float, **kwargs):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_c = T_s + 25
        result.append(T_c)
        return result
    def eqn_10_13__T_s(T_c: float, **kwargs):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_s = T_c - 25
        result.append(T_s)
        return result
    @kwasak_static
    def eqn_10_14(T_c=None, T_s=None, **kwargs):
        return

    def eqn_10_14__T_c(T_s: float, **kwargs):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_c = T_s + 12
        result.append(T_c)
        return result
    def eqn_10_14__T_s(T_c: float, **kwargs):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_s = T_c - 12
        result.append(T_s)
        return result
    @kwasak_static
    def eqn_10_15(P=None, S_Th=None, S_p=None, p_s=None, **kwargs):
        return

    def eqn_10_15__P(S_Th: float, S_p: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        P = S_Th*p_s/(S_Th - S_p)
        result.append(P)
        return result
    def eqn_10_15__S_Th(P: float, S_p: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_Th = P*S_p/(P - p_s)
        result.append(S_Th)
        return result
    def eqn_10_15__S_p(P: float, S_Th: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_p = S_Th*(P - p_s)/P
        result.append(S_p)
        return result
    def eqn_10_15__p_s(P: float, S_Th: float, S_p: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        p_s = P*(S_Th - S_p)/S_Th
        result.append(p_s)
        return result
    @kwasak_static
    def eqn_10_17(P=None, S_0=None, S_Th=None, p_0=None, p_s=None, **kwargs):
        return

    def eqn_10_17__P(S_0: float, S_Th: float, p_0: float, p_s: float, **kwargs):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        P = (p_0*(S_Th/S_0)**(5/3) - p_s)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = (p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/((-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = (p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/((-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result
    def eqn_10_17__S_0(P: float, S_Th: float, p_0: float, p_s: float, **kwargs):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/((P - p_s)/(P - p_0))**(3/5)
        result.append(S_0)
        return result
    def eqn_10_17__S_Th(P: float, S_0: float, p_0: float, p_s: float, **kwargs):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*((P - p_s)/(P - p_0))**(3/5)
        result.append(S_Th)
        return result
    def eqn_10_17__p_0(P: float, S_0: float, S_Th: float, p_s: float, **kwargs):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_0 = (P*(S_Th/S_0)**(5/3) - P + p_s)/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = (P*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = (P*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result
    def eqn_10_17__p_s(P: float, S_0: float, S_Th: float, p_0: float, **kwargs):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_s = -P*(S_Th/S_0)**(5/3) + P + p_0*(S_Th/S_0)**(5/3)
        result.append(p_s)
        p_s = -P*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 + P + p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        p_s = -P*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 + P + p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        return result
    @kwasak_static
    def eqn_10_18(P=None, S_Th=None, S_p=None, T_e=None, T_i=None, p_c=None, p_s=None, **kwargs):
        return

    def eqn_10_18__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        P = (S_Th*T_i*p_s + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*T_i + 460*S_Th - S_p*T_e - 460*S_p)
        result.append(P)
        return result
    def eqn_10_18__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_Th = S_p*(P*T_e + 460*P - T_e*p_c - 460*p_c)/(P*T_i + 460*P - T_i*p_s - 460*p_s)
        result.append(S_Th)
        return result
    def eqn_10_18__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_p = S_Th*(P*T_i + 460*P - T_i*p_s - 460*p_s)/(P*T_e + 460*P - T_e*p_c - 460*p_c)
        result.append(S_p)
        return result
    def eqn_10_18__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_e = (P*S_Th*T_i + 460*P*S_Th - 460*P*S_p - S_Th*T_i*p_s - 460*S_Th*p_s + 460*S_p*p_c)/(S_p*(P - p_c))
        result.append(T_e)
        return result
    def eqn_10_18__T_i(P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_i = (-460*P*S_Th + P*S_p*T_e + 460*P*S_p + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*(P - p_s))
        result.append(T_i)
        return result
    def eqn_10_18__p_c(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_c = (-P*S_Th*T_i - 460*P*S_Th + P*S_p*T_e + 460*P*S_p + S_Th*T_i*p_s + 460*S_Th*p_s)/(S_p*(T_e + 460))
        result.append(p_c)
        return result
    def eqn_10_18__p_s(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, **kwargs):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_s = (P*S_Th*T_i + 460*P*S_Th - P*S_p*T_e - 460*P*S_p + S_p*T_e*p_c + 460*S_p*p_c)/(S_Th*(T_i + 460))
        result.append(p_s)
        return result
    @kwasak_static
    def eqn_10_2(PS=None, Q_gas=None, V=None, dP=None, dt=None, **kwargs):
        return

    def eqn_10_2__PS(Q_gas: float, V: float, dP: float, dt: float, **kwargs):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        PS = Q_gas - V*dP/dt
        result.append(PS)
        return result
    def eqn_10_2__Q_gas(PS: float, V: float, dP: float, dt: float, **kwargs):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        Q_gas = PS + V*dP/dt
        result.append(Q_gas)
        return result
    def eqn_10_2__V(PS: float, Q_gas: float, dP: float, dt: float, **kwargs):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        V = dt*(-PS + Q_gas)/dP
        result.append(V)
        return result
    def eqn_10_2__dP(PS: float, Q_gas: float, V: float, dt: float, **kwargs):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dP = dt*(-PS + Q_gas)/V
        result.append(dP)
        return result
    def eqn_10_2__dt(PS: float, Q_gas: float, V: float, dP: float, **kwargs):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dt = -V*dP/(PS - Q_gas)
        result.append(dt)
        return result
    @kwasak_static
    def eqn_10_21(P=None, P_d=None, P_prime=None, **kwargs):
        return

    def eqn_10_21__P(P_d: float, P_prime: float, **kwargs):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P = P_d*P_prime/760
        result.append(P)
        return result
    def eqn_10_21__P_d(P: float, P_prime: float, **kwargs):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_d = 760*P/P_prime
        result.append(P_d)
        return result
    def eqn_10_21__P_prime(P: float, P_d: float, **kwargs):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_prime = 760*P/P_d
        result.append(P_prime)
        return result
    @kwasak_static
    def eqn_10_3(N_mfw=None, Q_gas=None, T=None, **kwargs):
        return

    def eqn_10_3__N_mfw(Q_gas: float, T: float, **kwargs):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        N_mfw = 0.108108108108108*Q_gas/T
        result.append(N_mfw)
        return result
    def eqn_10_3__Q_gas(N_mfw: float, T: float, **kwargs):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        Q_gas = 9.25*N_mfw*T
        result.append(Q_gas)
        return result
    def eqn_10_3__T(N_mfw: float, Q_gas: float, **kwargs):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        T = 0.108108108108108*Q_gas/N_mfw
        result.append(T)
        return result
    @kwasak_static
    def eqn_10_4(Q_gas=None, SP_1=None, SP_2=None, S_p=None, V=None, t=None, **kwargs):
        return

    def eqn_10_4__Q_gas(SP_1: float, SP_2: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        Q_gas = -(SP_1 - SP_2*exp(S_p*t/V))/(exp(S_p*t/V) - 1)
        result.append(Q_gas)
        return result
    def eqn_10_4__SP_1(Q_gas: float, SP_2: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_1 = Q_gas + (-Q_gas + SP_2)*exp(S_p*t/V)
        result.append(SP_1)
        return result
    def eqn_10_4__SP_2(Q_gas: float, SP_1: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_2 = (Q_gas*exp(S_p*t/V) - Q_gas + SP_1)*exp(-S_p*t/V)
        result.append(SP_2)
        return result
    def eqn_10_4__S_p(Q_gas: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        S_p = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/t
        result.append(S_p)
        return result
    def eqn_10_4__V(Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        V = S_p*t/log((Q_gas - SP_1)/(Q_gas - SP_2))
        result.append(V)
        return result
    def eqn_10_4__t(Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float, **kwargs):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        t = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/S_p
        result.append(t)
        return result
    @kwasak_static
    def eqn_10_5(P_1=None, P_2=None, S_p=None, V=None, t=None, **kwargs):
        return

    def eqn_10_5__P_1(P_2: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_p*t/V)
        result.append(P_1)
        return result
    def eqn_10_5__P_2(P_1: float, S_p: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_p*t/V)
        result.append(P_2)
        return result
    def eqn_10_5__S_p(P_1: float, P_2: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        S_p = V*log(P_1/P_2)/t
        result.append(S_p)
        return result
    def eqn_10_5__V(P_1: float, P_2: float, S_p: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        V = S_p*t/log(P_1/P_2)
        result.append(V)
        return result
    def eqn_10_5__t(P_1: float, P_2: float, S_p: float, V: float, **kwargs):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_p
        result.append(t)
        return result
    @kwasak_static
    def eqn_10_6(P_1=None, P_2=None, S_a=None, V=None, t=None, **kwargs):
        return

    def eqn_10_6__P_1(P_2: float, S_a: float, V: float, t: float, **kwargs):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_a*t/V)
        result.append(P_1)
        return result
    def eqn_10_6__P_2(P_1: float, S_a: float, V: float, t: float, **kwargs):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_a*t/V)
        result.append(P_2)
        return result
    def eqn_10_6__S_a(P_1: float, P_2: float, V: float, t: float, **kwargs):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        S_a = V*log(P_1/P_2)/t
        result.append(S_a)
        return result
    def eqn_10_6__V(P_1: float, P_2: float, S_a: float, t: float, **kwargs):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        V = S_a*t/log(P_1/P_2)
        result.append(V)
        return result
    def eqn_10_6__t(P_1: float, P_2: float, S_a: float, V: float, **kwargs):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_a
        result.append(t)
        return result
    @kwasak_static
    def eqn_10_8(bhp=None, c_p=None, delta_T=None, delta_h_i=None, f_a=None, rho=None, w_i=None, **kwargs):
        return

    def eqn_10_8__bhp(c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        bhp = 0.00315127701375246*c_p*delta_T*f_a*rho - 0.000392927308447937*delta_h_i*w_i
        result.append(bhp)
        return result
    def eqn_10_8__c_p(bhp: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        c_p = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(delta_T*f_a*rho)
        result.append(c_p)
        return result
    def eqn_10_8__delta_T(bhp: float, c_p: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_T = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*f_a*rho)
        result.append(delta_T)
        return result
    def eqn_10_8__delta_h_i(bhp: float, c_p: float, delta_T: float, f_a: float, rho: float, w_i: float, **kwargs):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_h_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/w_i
        result.append(delta_h_i)
        return result
    def eqn_10_8__f_a(bhp: float, c_p: float, delta_T: float, delta_h_i: float, rho: float, w_i: float, **kwargs):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        f_a = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*rho)
        result.append(f_a)
        return result
    def eqn_10_8__rho(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, w_i: float, **kwargs):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        rho = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*f_a)
        result.append(rho)
        return result
    def eqn_10_8__w_i(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, **kwargs):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        w_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/delta_h_i
        result.append(w_i)
        return result
    @kwasak_static
    def eqn_10_9(T_c=None, T_s=None, delta_T=None, **kwargs):
        return

    def eqn_10_9__T_c(T_s: float, delta_T: float, **kwargs):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_c = T_s + delta_T
        result.append(T_c)
        return result
    def eqn_10_9__T_s(T_c: float, delta_T: float, **kwargs):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_s = T_c - delta_T
        result.append(T_s)
        return result
    def eqn_10_9__delta_T(T_c: float, T_s: float, **kwargs):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        delta_T = T_c - T_s
        result.append(delta_T)
        return result
class Precondensors:
    @kwasak_static
    def eqn_7_1(P=None, p_i=None, y_i=None, **kwargs):
        return

    def eqn_7_1__P(p_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return result
    def eqn_7_1__p_i(P: float, y_i: float, **kwargs):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
        return result
    def eqn_7_1__y_i(P: float, p_i: float, **kwargs):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return result
    @kwasak_static
    def eqn_7_10(L_c_P=None, Q_condensor_heat_duty=None, del_T=None, **kwargs):
        return

    def eqn_7_10__L_c_P(Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty/(500*del_T)
        result.append(L_c_P)
        return result
    def eqn_7_10__Q_condensor_heat_duty(L_c_P: float, del_T: float, **kwargs):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500*L_c_P*del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_10__del_T(L_c_P: float, Q_condensor_heat_duty: float, **kwargs):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(500*L_c_P)
        result.append(del_T)
        return result
    @kwasak_static
    def eqn_7_11(Q_condensor_heat_duty=None, U_v=None, V_c=None, del_T_LM=None, **kwargs):
        return

    def eqn_7_11__Q_condensor_heat_duty(U_v: float, V_c: float, del_T_LM: float, **kwargs):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v*V_c*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_11__U_v(Q_condensor_heat_duty: float, V_c: float, del_T_LM: float, **kwargs):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
        result.append(U_v)
        return result
    def eqn_7_11__V_c(Q_condensor_heat_duty: float, U_v: float, del_T_LM: float, **kwargs):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty/(U_v*del_T_LM)
        result.append(V_c)
        return result
    def eqn_7_11__del_T_LM(Q_condensor_heat_duty: float, U_v: float, V_c: float, **kwargs):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(U_v*V_c)
        result.append(del_T_LM)
        return result
    @kwasak_static
    def eqn_7_12(A=None, Q_condensor_heat_duty=None, U=None, del_T=None, **kwargs):
        return

    def eqn_7_12__A(Q_condensor_heat_duty: float, U: float, del_T: float, **kwargs):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty/(U*del_T)
        result.append(A)
        return result
    def eqn_7_12__Q_condensor_heat_duty(A: float, U: float, del_T: float, **kwargs):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A*U*del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty/(A*del_T)
        result.append(U)
        return result
    def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float, **kwargs):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty/(A*U)
        result.append(del_T)
        return result
    @kwasak_static
    def eqn_7_14a(A=None, Q_condensor_heat_duty=None, U=None, del_T_LM=None, **kwargs):
        return

    def eqn_7_14a__A(Q_condensor_heat_duty: float, U: float, del_T_LM: float, **kwargs):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty/(U*del_T_LM)
        result.append(A)
        return result
    def eqn_7_14a__Q_condensor_heat_duty(A: float, U: float, del_T_LM: float, **kwargs):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A*U*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_14a__U(A: float, Q_condensor_heat_duty: float, del_T_LM: float, **kwargs):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty/(A*del_T_LM)
        result.append(U)
        return result
    def eqn_7_14a__del_T_LM(A: float, Q_condensor_heat_duty: float, U: float, **kwargs):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(A*U)
        result.append(del_T_LM)
        return result
    @kwasak_static
    def eqn_7_15(U=None, sum_R=None, **kwargs):
        return

    def eqn_7_15__U(sum_R: float, **kwargs):
        # [.pyeqn] 1 / U = sum_R
        result = []
        U = 1/sum_R
        result.append(U)
        return result
    def eqn_7_15__sum_R(U: float, **kwargs):
        # [.pyeqn] 1 / U = sum_R
        result = []
        sum_R = 1/U
        result.append(sum_R)
        return result
    @kwasak_static
    def eqn_7_16(D_0=None, D_LM=None, D_i=None, R_f_0=None, R_fi=None, U_0=None, h_0=None, h_i=None, k_w=None, x_w=None, **kwargs):
        return

    def eqn_7_16__D_0(D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_f_0*U_0*h_0 - U_0 + h_0)/(U_0*h_0*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result
    def eqn_7_16__D_LM(D_0: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_0*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(D_LM)
        return result
    def eqn_7_16__D_i(D_0: float, D_LM: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_0*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_0*x_w + D_LM*R_f_0*U_0*h_0*k_w + D_LM*U_0*k_w - D_LM*h_0*k_w))
        result.append(D_i)
        return result
    def eqn_7_16__R_f_0(D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
        result.append(R_f_0)
        return result
    def eqn_7_16__R_fi(D_0: float, D_LM: float, D_i: float, R_f_0: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_f_0/D_0 - D_i/(D_0*h_0) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result
    def eqn_7_16__U_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_0*h_i*k_w/(D_0*D_LM*R_fi*h_0*h_i*k_w + D_0*D_LM*h_0*k_w + D_0*D_i*h_0*h_i*x_w + D_LM*D_i*R_f_0*h_0*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result
    def eqn_7_16__h_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_0 = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_f_0*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_0)
        return result
    def eqn_7_16__h_i(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_0*k_w/(D_0*D_LM*R_fi*U_0*h_0*k_w + D_0*D_i*U_0*h_0*x_w + D_LM*D_i*R_f_0*U_0*h_0*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_0*k_w)
        result.append(h_i)
        return result
    def eqn_7_16__k_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_0*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(k_w)
        return result
    def eqn_7_16__x_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_f_0*k_w/D_0 - D_LM*k_w/(D_0*h_0) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result
    @kwasak_static
    def eqn_7_17(R_0=None, R_nc=None, h_c=None, **kwargs):
        return

    def eqn_7_17__R_0(R_nc: float, h_c: float, **kwargs):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1/h_c
        result.append(R_0)
        return result
    def eqn_7_17__R_nc(R_0: float, h_c: float, **kwargs):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1/h_c
        result.append(R_nc)
        return result
    def eqn_7_17__h_c(R_0: float, R_nc: float, **kwargs):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1/(R_0 - R_nc)
        result.append(h_c)
        return result
    @kwasak_static
    def eqn_7_18(D_0=None, D_LM=None, D_i=None, R_fi=None, R_fo=None, R_nc=None, U_0=None, h_c=None, h_i=None, k_w=None, x_w=None, **kwargs):
        return

    def eqn_7_18__D_0(D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_fo*U_0*h_c - R_nc*U_0*h_c - U_0 + h_c)/(U_0*h_c*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result
    def eqn_7_18__D_LM(D_0: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_c*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(D_LM)
        return result
    def eqn_7_18__D_i(D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_c*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_c*x_w + D_LM*R_fo*U_0*h_c*k_w + D_LM*R_nc*U_0*h_c*k_w + D_LM*U_0*k_w - D_LM*h_c*k_w))
        result.append(D_i)
        return result
    def eqn_7_18__R_fi(D_0: float, D_LM: float, D_i: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_fo/D_0 - D_i*R_nc/D_0 - D_i/(D_0*h_c) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result
    def eqn_7_18__R_fo(D_0: float, D_LM: float, D_i: float, R_fi: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fo = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_nc - 1/h_c + 1/U_0
        result.append(R_fo)
        return result
    def eqn_7_18__R_nc(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_nc = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_fo - 1/h_c + 1/U_0
        result.append(R_nc)
        return result
    def eqn_7_18__U_0(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_c*h_i*k_w/(D_0*D_LM*R_fi*h_c*h_i*k_w + D_0*D_LM*h_c*k_w + D_0*D_i*h_c*h_i*x_w + D_LM*D_i*R_fo*h_c*h_i*k_w + D_LM*D_i*R_nc*h_c*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result
    def eqn_7_18__h_c(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_c = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_fo*U_0*h_i*k_w + D_LM*D_i*R_nc*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_c)
        return result
    def eqn_7_18__h_i(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, k_w: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_c*k_w/(D_0*D_LM*R_fi*U_0*h_c*k_w + D_0*D_i*U_0*h_c*x_w + D_LM*D_i*R_fo*U_0*h_c*k_w + D_LM*D_i*R_nc*U_0*h_c*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_c*k_w)
        result.append(h_i)
        return result
    def eqn_7_18__k_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, x_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_c*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(k_w)
        return result
    def eqn_7_18__x_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, **kwargs):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_fo*k_w/D_0 - D_LM*R_nc*k_w/D_0 - D_LM*k_w/(D_0*h_c) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result
    @kwasak_static
    def eqn_7_2(P_i_0=None, p_i=None, x_i=None, **kwargs):
        return

    def eqn_7_2__P_i_0(p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i/x_i
        result.append(P_i_0)
        return result
    def eqn_7_2__p_i(P_i_0: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        p_i = P_i_0*x_i
        result.append(p_i)
        return result
    def eqn_7_2__x_i(P_i_0: float, p_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        x_i = p_i/P_i_0
        result.append(x_i)
        return result
    @kwasak_static
    def eqn_7_3(P_i_0=None, epsilon_i=None, p_i=None, x_i=None, **kwargs):
        return

    def eqn_7_3__P_i_0(epsilon_i: float, p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i/(epsilon_i*x_i)
        result.append(P_i_0)
        return result
    def eqn_7_3__epsilon_i(P_i_0: float, p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i/(P_i_0*x_i)
        result.append(epsilon_i)
        return result
    def eqn_7_3__p_i(P_i_0: float, epsilon_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0*epsilon_i*x_i
        result.append(p_i)
        return result
    def eqn_7_3__x_i(P_i_0: float, epsilon_i: float, p_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i/(P_i_0*epsilon_i)
        result.append(x_i)
        return result
    @kwasak_static
    def eqn_7_4a(P=None, p_c=None, p_nc=None, **kwargs):
        return

    def eqn_7_4a__P(p_c: float, p_nc: float, **kwargs):
        # [.pyeqn] p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return result
    def eqn_7_4a__p_c(P: float, p_nc: float, **kwargs):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return result
    def eqn_7_4a__p_nc(P: float, p_c: float, **kwargs):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return result
    @kwasak_static
    def eqn_7_4aa(n_i=None, n_nc=None, p_i=None, p_nc=None, **kwargs):
        return

    def eqn_7_4aa__n_i(n_nc: float, p_i: float, p_nc: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_i = n_nc*p_i/p_nc
        result.append(n_i)
        return result
    def eqn_7_4aa__n_nc(n_i: float, p_i: float, p_nc: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        n_nc = n_i*p_nc/p_i
        result.append(n_nc)
        return result
    def eqn_7_4aa__p_i(n_i: float, n_nc: float, p_nc: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_i = n_i*p_nc/n_nc
        result.append(p_i)
        return result
    def eqn_7_4aa__p_nc(n_i: float, n_nc: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / p_nc
        result = []
        p_nc = n_nc*p_i/n_i
        result.append(p_nc)
        return result
    @kwasak_static
    def eqn_7_4ac(P_c=None, n_i=None, n_nc=None, p=None, p_i=None, **kwargs):
        return

    def eqn_7_4ac__P_c(n_i: float, n_nc: float, p: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result
    def eqn_7_4ac__n_i(P_c: float, n_nc: float, p: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result
    def eqn_7_4ac__n_nc(P_c: float, n_i: float, p: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result
    def eqn_7_4ac__p(P_c: float, n_i: float, n_nc: float, p_i: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result
    def eqn_7_4ac__p_i(P_c: float, n_i: float, n_nc: float, p: float, **kwargs):
        # [.pyeqn] n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result
    @kwasak_static
    def eqn_7_5(N_i=None, N_nc=None, P=None, P_c=None, p_i=None, **kwargs):
        return

    def eqn_7_5__N_i(N_nc: float, P: float, P_c: float, p_i: float, **kwargs):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
        return result
    def eqn_7_5__N_nc(N_i: float, P: float, P_c: float, p_i: float, **kwargs):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return result
    def eqn_7_5__P_c(N_i: float, N_nc: float, P: float, p_i: float, **kwargs):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc*p_i/N_i
        result.append(P_c)
        return result
    def eqn_7_5__P(N_i: float, N_nc: float, P_c: float, p_i: float, **kwargs):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc*p_i/N_i
        result.append(P)
        return result
    def eqn_7_5__p_i(N_i: float, N_nc: float, P: float, P_c: float, **kwargs):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i*(P - P_c)/N_nc
        result.append(p_i)
        return result
    @kwasak_static
    def eqn_7_6(M=None, P=None, P_i_0=None, W_air=None, W_i=None, p_c=None, x_i=None, **kwargs):
        return

    def eqn_7_6__M(P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*x_i)
        result.append(M)
        return result
    def eqn_7_6__P(M: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*x_i/(29*W_i) + p_c
        result.append(P)
        return result
    def eqn_7_6__P_i_0(M: float, P: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*x_i)
        result.append(P_i_0)
        return result
    def eqn_7_6__W_air(M: float, P: float, P_i_0: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*x_i)
        result.append(W_air)
        return result
    def eqn_7_6__W_i(M: float, P: float, P_i_0: float, W_air: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*x_i/(29*(P - p_c))
        result.append(W_i)
        return result
    def eqn_7_6__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*x_i/(29*W_i) + P
        result.append(p_c)
        return result
    def eqn_7_6__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air)
        result.append(x_i)
        return result
    @kwasak_static
    def eqn_7_7(M=None, P=None, P_i_0=None, W_air=None, W_i=None, epsilon_i=None, p_c=None, x_i=None, **kwargs):
        return

    def eqn_7_7__M(P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
        result.append(M)
        return result
    def eqn_7_7__P(M: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + p_c
        result.append(P)
        return result
    def eqn_7_7__P_i_0(M: float, P: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*epsilon_i*x_i)
        result.append(P_i_0)
        return result
    def eqn_7_7__W_air(M: float, P: float, P_i_0: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*epsilon_i*x_i)
        result.append(W_air)
        return result
    def eqn_7_7__W_i(M: float, P: float, P_i_0: float, W_air: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*epsilon_i*x_i/(29*(P - p_c))
        result.append(W_i)
        return result
    def eqn_7_7__epsilon_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*x_i)
        result.append(epsilon_i)
        return result
    def eqn_7_7__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, x_i: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + P
        result.append(p_c)
        return result
    def eqn_7_7__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, **kwargs):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*epsilon_i)
        result.append(x_i)
        return result
    @kwasak_static
    def eqn_7_8(L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, **kwargs):
        return

    def eqn_7_8__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return result
    def eqn_7_8__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_8__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return result
    def eqn_7_8__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return result
    @kwasak_static
    def eqn_7_9(L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, rho=None, **kwargs):
        return

    def eqn_7_9__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return result
    def eqn_7_9__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, rho: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_9__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return result
    def eqn_7_9__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return result
    def eqn_7_9__rho(L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float, **kwargs):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return result
class PressMgmt:
    @kwasak_static
    def eqn_3_1(Abs_Pressure=None, BarometricPressure=None, Vacuum=None, **kwargs):
        return

    def eqn_3_1__Abs_Pressure(BarometricPressure: float, Vacuum: float, **kwargs):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return result
    def eqn_3_1__BarometricPressure(Abs_Pressure: float, Vacuum: float, **kwargs):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return result
    def eqn_3_1__Vacuum(Abs_Pressure: float, BarometricPressure: float, **kwargs):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return result
    @kwasak_static
    def eqn_3_13(H_1=None, H_2=None, KAPPA_2=None, P=None, **kwargs):
        return

    def eqn_3_13__H_1(H_2: float, KAPPA_2: float, P: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_1 = H_2 - P/KAPPA_2
        result.append(H_1)
        return result
    def eqn_3_13__H_2(H_1: float, KAPPA_2: float, P: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_2 = H_1 + P/KAPPA_2
        result.append(H_2)
        return result
    def eqn_3_13__KAPPA_2(H_1: float, H_2: float, P: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        KAPPA_2 = -P/(H_1 - H_2)
        result.append(KAPPA_2)
        return result
    def eqn_3_13__P(H_1: float, H_2: float, KAPPA_2: float, **kwargs):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        P = KAPPA_2*(-H_1 + H_2)
        result.append(P)
        return result
    @kwasak_static
    def eqn_3_2(G=None, G_C=None, H=None, P=None, rho=None, **kwargs):
        return

    def eqn_3_2__G_C(G: float, H: float, P: float, rho: float, **kwargs):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return result
    def eqn_3_2__G(G_C: float, H: float, P: float, rho: float, **kwargs):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return result
    def eqn_3_2__H(G: float, G_C: float, P: float, rho: float, **kwargs):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return result
    def eqn_3_2__P(G: float, G_C: float, H: float, rho: float, **kwargs):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return result
    def eqn_3_2__rho(G: float, G_C: float, H: float, P: float, **kwargs):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return result
    @kwasak_static
    def eqn_3_3(H_1=None, H_2=None, P=None, P_P=None, **kwargs):
        return

    def eqn_3_3__H_1(H_2: float, P: float, P_P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return result
    def eqn_3_3__H_2(H_1: float, P: float, P_P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 - P + P_P
        result.append(H_2)
        return result
    def eqn_3_3__P_P(H_1: float, H_2: float, P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return result
    def eqn_3_3__P(H_1: float, H_2: float, P_P: float, **kwargs):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P = H_1 - H_2 + P_P
        result.append(P)
        return result
    @kwasak_static
    def eqn_3_4(KAPPA=None, P=None, V=None, **kwargs):
        return

    def eqn_3_4__KAPPA(P: float, V: float, **kwargs):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return result
    def eqn_3_4__P(KAPPA: float, V: float, **kwargs):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return result
    def eqn_3_4__V(KAPPA: float, P: float, **kwargs):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return result
    @kwasak_static
    def eqn_3_5(P=None, P_P=None, V=None, V_P=None, **kwargs):
        return

    def eqn_3_5__P_P(P: float, V: float, V_P: float, **kwargs):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return result
    def eqn_3_5__P(P_P: float, V: float, V_P: float, **kwargs):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return result
    def eqn_3_5__V_P(P: float, P_P: float, V: float, **kwargs):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return result
    def eqn_3_5__V(P: float, P_P: float, V_P: float, **kwargs):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return result
    @kwasak_static
    def eqn_3_6(H_1=None, H_2=None, P=None, V=None, V_P=None, **kwargs):
        return

    def eqn_3_6__H_1(H_2: float, P: float, V: float, V_P: float, **kwargs):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return result
    def eqn_3_6__H_2(H_1: float, P: float, V: float, V_P: float, **kwargs):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P*V/V_P - P
        result.append(H_2)
        return result
    def eqn_3_6__P(H_1: float, H_2: float, V: float, V_P: float, **kwargs):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P*(-H_1 + H_2)/(V - V_P)
        result.append(P)
        return result
    def eqn_3_6__V_P(H_1: float, H_2: float, P: float, V: float, **kwargs):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return result
    def eqn_3_6__V(H_1: float, H_2: float, P: float, V_P: float, **kwargs):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return result
    @kwasak_static
    def eqn_3_8(A_C=None, H_2=None, V_P=None, **kwargs):
        return

    def eqn_3_8__A_C(H_2: float, V_P: float, **kwargs):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        A_C = V_P/H_2
        result.append(A_C)
        return result
    def eqn_3_8__H_2(A_C: float, V_P: float, **kwargs):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        H_2 = V_P/A_C
        result.append(H_2)
        return result
    def eqn_3_8__V_P(A_C: float, H_2: float, **kwargs):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        V_P = A_C*H_2
        result.append(V_P)
        return result
class ProcessApp1:
    @kwasak_static
    def eqn_5_1(K_i=None, x_i=None, y_i=None, **kwargs):
        return

    def eqn_5_1__K_i(x_i: float, y_i: float, **kwargs):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i/x_i
        result.append(K_i)
        return result
    def eqn_5_1__x_i(K_i: float, y_i: float, **kwargs):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
        return result
    def eqn_5_1__y_i(K_i: float, x_i: float, **kwargs):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
        result.append(y_i)
        return result
    @kwasak_static
    def eqn_5_10b(L_0=None, R=None, V_1=None, **kwargs):
        return

    def eqn_5_10b__L_0(R: float, V_1: float, **kwargs):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return result
    def eqn_5_10b__R(L_0: float, V_1: float, **kwargs):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return result
    def eqn_5_10b__V_1(L_0: float, R: float, **kwargs):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return result
    @kwasak_static
    def eqn_5_10c(D=None, L_0=None, R=None, **kwargs):
        return

    def eqn_5_10c__D(L_0: float, R: float, **kwargs):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result
    def eqn_5_10c__L_0(D: float, R: float, **kwargs):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
        return result
    def eqn_5_10c__R(D: float, L_0: float, **kwargs):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
        result.append(R)
        return result
    @kwasak_static
    def eqn_5_11(B=None, L_N=None, V_0=None, **kwargs):
        return

    def eqn_5_11__B(L_N: float, V_0: float, **kwargs):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        B = L_N - V_0
        result.append(B)
        return result
    def eqn_5_11__L_N(B: float, V_0: float, **kwargs):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        L_N = B + V_0
        result.append(L_N)
        return result
    def eqn_5_11__V_0(B: float, L_N: float, **kwargs):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        V_0 = -B + L_N
        result.append(V_0)
        return result
    @kwasak_static
    def eqn_5_12(Eff=None, N_ES=None, N_t=None, T=None, **kwargs):
        return

    def eqn_5_12__Eff(N_ES: float, N_t: float, T: float, **kwargs):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result
    def eqn_5_12__N_ES(Eff: float, N_t: float, T: float, **kwargs):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result
    def eqn_5_12__N_t(Eff: float, N_ES: float, T: float, **kwargs):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result
    def eqn_5_12__T(Eff: float, N_ES: float, N_t: float, **kwargs):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result
    @kwasak_static
    def eqn_5_13(HETP=None, H_p=None, N_ES=None, **kwargs):
        return

    def eqn_5_13__HETP(H_p: float, N_ES: float, **kwargs):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result
    def eqn_5_13__H_p(HETP: float, N_ES: float, **kwargs):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result
    def eqn_5_13__N_ES(HETP: float, H_p: float, **kwargs):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result
    @kwasak_static
    def eqn_5_14(M=None, P_0=None, T=None, W_E=None, **kwargs):
        return

    def eqn_5_14__M(P_0: float, T: float, W_E: float, **kwargs):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return result
    def eqn_5_14__P_0(M: float, T: float, W_E: float, **kwargs):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return result
    def eqn_5_14__T(M: float, P_0: float, W_E: float, **kwargs):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return result
    def eqn_5_14__W_E(M: float, P_0: float, T: float, **kwargs):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return result
    @kwasak_static
    def eqn_5_16(H_i=None, p_i=None, x_i=None, **kwargs):
        return

    def eqn_5_16__H_i(p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result
    def eqn_5_16__p_i(H_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result
    def eqn_5_16__x_i(H_i: float, p_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result
    @kwasak_static
    def eqn_5_17(H_2_1=None, H_2_3=None, H_2_mi=None, x_1=None, x_3=None, **kwargs):
        return

    def eqn_5_17__H_2_1(H_2_3: float, H_2_mi: float, x_1: float, x_3: float, **kwargs):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3*log(H_2_3) + log(H_2_mi))/x_1)
        result.append(H_2_1)
        return result
    def eqn_5_17__H_2_3(H_2_1: float, H_2_mi: float, x_1: float, x_3: float, **kwargs):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1*log(H_2_1) + log(H_2_mi))/x_3)
        result.append(H_2_3)
        return result
    def eqn_5_17__H_2_mi(H_2_1: float, H_2_3: float, x_1: float, x_3: float, **kwargs):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1*log(H_2_1) + x_3*log(H_2_3))
        result.append(H_2_mi)
        return result
    def eqn_5_17__x_1(H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float, **kwargs):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return result
    def eqn_5_17__x_3(H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float, **kwargs):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return result
    @kwasak_static
    def eqn_5_2a(K_1=None, K_2=None, alpha_1_2=None, **kwargs):
        return

    def eqn_5_2a__K_1(K_2: float, alpha_1_2: float, **kwargs):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2*alpha_1_2
        result.append(K_1)
        return result
    def eqn_5_2a__K_2(K_1: float, alpha_1_2: float, **kwargs):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1/alpha_1_2
        result.append(K_2)
        return result
    def eqn_5_2a__alpha_1_2(K_1: float, K_2: float, **kwargs):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1/K_2
        result.append(alpha_1_2)
        return result
    @kwasak_static
    def eqn_5_2b(K_1=None, K_2=None, x_1=None, x_2=None, y_1=None, y_2=None, **kwargs):
        return

    def eqn_5_2b__K_1(K_2: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result
    def eqn_5_2b__K_2(K_1: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result
    def eqn_5_2b__x_1(K_1: float, K_2: float, x_2: float, y_1: float, y_2: float, **kwargs):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result
    def eqn_5_2b__x_2(K_1: float, K_2: float, x_1: float, y_1: float, y_2: float, **kwargs):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result
    def eqn_5_2b__y_1(K_1: float, K_2: float, x_1: float, x_2: float, y_2: float, **kwargs):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result
    def eqn_5_2b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float, **kwargs):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
        result.append(y_2)
        return result
    @kwasak_static
    def eqn_5_3(P_0_i=None, p_i=None, x_i=None, **kwargs):
        return

    def eqn_5_3__P_0_i(p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return result
    def eqn_5_3__p_i(P_0_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return result
    def eqn_5_3__x_i(P_0_i: float, p_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return result
    @kwasak_static
    def eqn_5_4(P=None, P_0_i=None, x_i=None, y_i=None, **kwargs):
        return

    def eqn_5_4__P_0_i(P: float, x_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P_0_i = P*y_i/x_i
        result.append(P_0_i)
        return result
    def eqn_5_4__P(P_0_i: float, x_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return result
    def eqn_5_4__x_i(P: float, P_0_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        x_i = P*y_i/P_0_i
        result.append(x_i)
        return result
    def eqn_5_4__y_i(P: float, P_0_i: float, x_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
        return result
    @kwasak_static
    def eqn_5_5(P_0_1=None, P_0_2=None, alpha_12=None, **kwargs):
        return

    def eqn_5_5__P_0_1(P_0_2: float, alpha_12: float, **kwargs):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2*alpha_12
        result.append(P_0_1)
        return result
    def eqn_5_5__P_0_2(P_0_1: float, alpha_12: float, **kwargs):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1/alpha_12
        result.append(P_0_2)
        return result
    def eqn_5_5__alpha_12(P_0_1: float, P_0_2: float, **kwargs):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return result
    @kwasak_static
    def eqn_5_6(P_0_i=None, gamma_i=None, p_i=None, x_i=None, **kwargs):
        return

    def eqn_5_6__P_0_i(gamma_i: float, p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result
    def eqn_5_6__gamma_i(P_0_i: float, p_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result
    def eqn_5_6__p_i(P_0_i: float, gamma_i: float, x_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return result
    def eqn_5_6__x_i(P_0_i: float, gamma_i: float, p_i: float, **kwargs):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result
    @kwasak_static
    def eqn_5_7(P=None, P_0_i=None, gamma_i=None, x_i=None, y_i=None, **kwargs):
        return

    def eqn_5_7__P_0_i(P: float, gamma_i: float, x_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P*y_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result
    def eqn_5_7__P(P_0_i: float, gamma_i: float, x_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i*gamma_i*x_i/y_i
        result.append(P)
        return result
    def eqn_5_7__gamma_i(P: float, P_0_i: float, x_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P*y_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result
    def eqn_5_7__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result
    def eqn_5_7__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float, **kwargs):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return result
    @kwasak_static
    def eqn_5_8(P_0_1=None, P_0_2=None, alpha_12=None, gamma_1=None, gamma_2=None, **kwargs):
        return

    def eqn_5_8__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float, **kwargs):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result
    def eqn_5_8__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float, **kwargs):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result
    def eqn_5_8__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float, **kwargs):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result
    def eqn_5_8__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float, **kwargs):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return result
    def eqn_5_8__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float, **kwargs):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return result
class ProcessApp2:
    @kwasak_static
    def eqn_6_1(T_1=None, T_2=None, T_R=None, c_p=None, del_h_v=None, w_1=None, w_2=None, w_v=None, **kwargs):
        return

    def eqn_6_1__T_1(T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R) + del_h_v*w_v)/(c_p*w_1)
        result.append(T_1)
        return result
    def eqn_6_1__T_2(T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_2 = (T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R) + del_h_v*w_v)/(c_p*w_2)
        result.append(T_2)
        return result
    def eqn_6_1__T_R(T_1: float, T_2: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_R = (T_1*c_p*w_1 + T_2*c_p*w_2 - del_h_v*w_v)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result
    def eqn_6_1__c_p(T_1: float, T_2: float, T_R: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = del_h_v*w_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result
    def eqn_6_1__del_h_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/w_v
        result.append(del_h_v)
        return result
    def eqn_6_1__w_1(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_2: float, w_v: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (-T_2*c_p*w_2 + T_R*c_p*w_2 + del_h_v*w_v)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result
    def eqn_6_1__w_2(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_v: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_2 = (-T_1*c_p*w_1 + T_R*c_p*w_1 + del_h_v*w_v)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result
    def eqn_6_1__w_v(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/del_h_v
        result.append(w_v)
        return result
    @kwasak_static
    def eqn_6_10(A=None, dV_dt=None, delta_P=None, mu=None, r_c=None, s=None, tau=None, **kwargs):
        return

    def eqn_6_10__A(dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float, **kwargs):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt*delta_P**(s - 1)*mu*r_c*tau
        result.append(A)
        return result
    def eqn_6_10__dV_dt(A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float, **kwargs):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A*delta_P**(1 - s)/(mu*r_c*tau)
        result.append(dV_dt)
        return result
    def eqn_6_10__delta_P(A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float, **kwargs):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
        result.append(delta_P)
        return result
    def eqn_6_10__mu(A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float, **kwargs):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A*delta_P**(1 - s)/(dV_dt*r_c*tau)
        result.append(mu)
        return result
    def eqn_6_10__r_c(A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float, **kwargs):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
        result.append(r_c)
        return result
    def eqn_6_10__s(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float, **kwargs):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A*delta_P/(dV_dt*mu*r_c*tau))/log(delta_P)
        result.append(s)
        return result
    def eqn_6_10__tau(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, **kwargs):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A*delta_P**(1 - s)/(dV_dt*mu*r_c)
        result.append(tau)
        return result
    @kwasak_static
    def eqn_6_11a(A_d=None, delta_T=None, delta_h_i=None, delta_m=None, h_d=None, m_b=None, t_R=None, **kwargs):
        return

    def eqn_6_11a__A_d(delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
        result.append(A_d)
        return result
    def eqn_6_11a__delta_T(A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
        result.append(delta_T)
        return result
    def eqn_6_11a__delta_h_i(A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d*delta_T*h_d*t_R/(delta_m*m_b)
        result.append(delta_h_i)
        return result
    def eqn_6_11a__delta_m(A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d*delta_T*h_d*t_R/(delta_h_i*m_b)
        result.append(delta_m)
        return result
    def eqn_6_11a__h_d(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, m_b: float, t_R: float, **kwargs):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i*delta_m*m_b/(A_d*delta_T*t_R)
        result.append(h_d)
        return result
    def eqn_6_11a__m_b(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float, **kwargs):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
        result.append(m_b)
        return result
    def eqn_6_11a__t_R(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, **kwargs):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        t_R = delta_h_i*delta_m*m_b/(A_d*delta_T*h_d)
        result.append(t_R)
        return result
    @kwasak_static
    def eqn_6_2(Q_v=None, T_1=None, T_2=None, T_R=None, c_p=None, w_1=None, w_2=None, **kwargs):
        return

    def eqn_6_2__Q_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/12000
        result.append(Q_v)
        return result
    def eqn_6_2__T_1(Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000*Q_v + T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R))/(c_p*w_1)
        result.append(T_1)
        return result
    def eqn_6_2__T_2(Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_2 = (12000*Q_v + T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R))/(c_p*w_2)
        result.append(T_2)
        return result
    def eqn_6_2__T_R(Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_R = (-12000*Q_v + T_1*c_p*w_1 + T_2*c_p*w_2)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result
    def eqn_6_2__c_p(Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        c_p = 12000*Q_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result
    def eqn_6_2__w_1(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_2: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_1 = (12000*Q_v - T_2*c_p*w_2 + T_R*c_p*w_2)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result
    def eqn_6_2__w_2(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, **kwargs):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000*Q_v - T_1*c_p*w_1 + T_R*c_p*w_1)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result
    @kwasak_static
    def eqn_6_4(Q_v=None, delta_h_v=None, w_v=None, **kwargs):
        return

    def eqn_6_4__Q_v(delta_h_v: float, w_v: float, **kwargs):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        Q_v = delta_h_v*w_v/12000
        result.append(Q_v)
        return result
    def eqn_6_4__delta_h_v(Q_v: float, w_v: float, **kwargs):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        delta_h_v = 12000*Q_v/w_v
        result.append(delta_h_v)
        return result
    def eqn_6_4__w_v(Q_v: float, delta_h_v: float, **kwargs):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return result
    @kwasak_static
    def eqn_6_7(C_1=None, C_2=None, T_1=None, T_2=None, c_p=None, delta_h_c=None, delta_h_v=None, m_b=None, m_v=None, **kwargs):
        return

    def eqn_6_7__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result
    def eqn_6_7__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result
    def eqn_6_7__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*m_v)/(c_p*m_b)
        result.append(T_1)
        return result
    def eqn_6_7__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*m_v)/(c_p*m_b)
        result.append(T_2)
        return result
    def eqn_6_7__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*m_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result
    def eqn_6_7__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*m_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result
    def eqn_6_7__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/m_v
        result.append(delta_h_v)
        return result
    def eqn_6_7__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_v: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_b = delta_h_v*m_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result
    def eqn_6_7__m_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, **kwargs):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/delta_h_v
        result.append(m_v)
        return result
    @kwasak_static
    def eqn_6_8(C_1=None, C_2=None, T_1=None, T_2=None, c_p=None, delta_h_c=None, delta_h_v=None, delta_t=None, m_b=None, w_v=None, **kwargs):
        return

    def eqn_6_8__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result
    def eqn_6_8__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result
    def eqn_6_8__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_1)
        return result
    def eqn_6_8__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_2)
        return result
    def eqn_6_8__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*delta_t*w_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result
    def eqn_6_8__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*delta_t*w_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result
    def eqn_6_8__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_t*w_v)
        result.append(delta_h_v)
        return result
    def eqn_6_8__delta_t(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_t = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*w_v)
        result.append(delta_t)
        return result
    def eqn_6_8__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, w_v: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        m_b = delta_h_v*delta_t*w_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result
    def eqn_6_8__w_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, **kwargs):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        w_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*delta_t)
        result.append(w_v)
        return result
class RotaryPistonVane:
    @kwasak_static
    def eqn_11_1(PS=None, Q_0=None, Q_external_gas_throughput=None, V=None, dP=None, dT=None, **kwargs):
        return

    def eqn_11_1__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float, **kwargs):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result
    def eqn_11_1__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float, **kwargs):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result
    def eqn_11_1__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float, **kwargs):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result
    def eqn_11_1__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float, **kwargs):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result
    def eqn_11_1__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float, **kwargs):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result
    def eqn_11_1__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, **kwargs):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result
    @kwasak_static
    def eqn_11_2(Q=None, Q_0=None, Q_external_gas_throughput=None, SP_1=None, SP_2=None, S_vol_pump_speed=None, V=None, t=None, **kwargs):
        return

    def eqn_11_2__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q_0 = -(Q*exp(S_vol_pump_speed*t/V) - Q_external_gas_throughput + SP_1 - SP_2*exp(S_vol_pump_speed*t/V))/(exp(S_vol_pump_speed*t/V) - 1)
        result.append(Q_0)
        return result
    def eqn_11_2__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float, **kwargs):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q = (Q_0 + Q_external_gas_throughput - SP_1 + (-Q_0 + SP_2)*exp(S_vol_pump_speed*t/V))*exp(-S_vol_pump_speed*t/V)
        result.append(Q)
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
    @kwasak_static
    def eqn_11_3(F_s=None, t=None, t_c=None, **kwargs):
        return

    def eqn_11_3__F_s(t: float, t_c: float, **kwargs):
        # [.pyeqn] t = t_c * F_s
        result = []
        F_s = t/t_c
        result.append(F_s)
        return result
    def eqn_11_3__t(F_s: float, t_c: float, **kwargs):
        # [.pyeqn] t = t_c * F_s
        result = []
        t = F_s*t_c
        result.append(t)
        return result
    def eqn_11_3__t_c(F_s: float, t: float, **kwargs):
        # [.pyeqn] t = t_c * F_s
        result = []
        t_c = t/F_s
        result.append(t_c)
        return result
    @kwasak_static
    def eqn_11_5(P_0_v=None, P_D=None, p_g=None, p_v_max=None, **kwargs):
        return

    def eqn_11_5__P_0_v(P_D: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D*p_v_max/(p_g + p_v_max)
        result.append(P_0_v)
        return result
    def eqn_11_5__P_D(P_0_v: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v*(p_g + p_v_max)/p_v_max
        result.append(P_D)
        return result
    def eqn_11_5__p_g(P_0_v: float, P_D: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max*(-P_0_v + P_D)/P_0_v
        result.append(p_g)
        return result
    def eqn_11_5__p_v_max(P_0_v: float, P_D: float, p_g: float, **kwargs):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v*p_g/(P_0_v - P_D)
        result.append(p_v_max)
        return result
    @kwasak_static
    def eqn_11_6(P_0_V=None, P_D=None, P_v_0=None, S_B=None, S_D=None, p_b=None, p_g=None, p_v_max=None, **kwargs):
        return

    def eqn_11_6__P_0_V(P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_g - P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(P_0_V)
        return result
    def eqn_11_6__P_D(P_0_V: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0*S_D*(p_g + p_v_max)/(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)
        result.append(P_D)
        return result
    def eqn_11_6__P_v_0(P_0_V: float, P_D: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = P_D*(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)/(S_D*(p_g + p_v_max))
        result.append(P_v_0)
        return result
    def eqn_11_6__S_B(P_0_V: float, P_D: float, P_v_0: float, S_D: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = S_D*(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)/(P_D*(P_0_V - p_b))
        result.append(S_B)
        return result
    def eqn_11_6__S_D(P_0_V: float, P_D: float, P_v_0: float, S_B: float, p_b: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = P_D*S_B*(P_0_V - p_b)/(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)
        result.append(S_D)
        return result
    def eqn_11_6__p_b(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_g: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (P_0_V*P_D*S_B - P_D*S_D*p_v_max + P_v_0*S_D*p_g + P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(p_b)
        return result
    def eqn_11_6__p_g(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_v_max: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (-P_0_V*P_D*S_B + P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_v_max)/(P_v_0*S_D)
        result.append(p_g)
        return result
    def eqn_11_6__p_v_max(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, **kwargs):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V*P_D*S_B - P_D*S_B*p_b + P_v_0*S_D*p_g)/(S_D*(P_D - P_v_0))
        result.append(p_v_max)
        return result
class SelectingPump:
    @kwasak_static
    def eqn_8_2(hp=None, installed_costs=None, **kwargs):
        return

    def eqn_8_2__hp(installed_costs: float, **kwargs):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        hp = 9.18273645546364e-9*installed_costs**2
        result.append(hp)
        return result
    def eqn_8_2__installed_costs(hp: float, **kwargs):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        installed_costs = 10435.5162785557*sqrt(hp)
        result.append(installed_costs)
        return result
    @kwasak_static
    def eqn_8_5(Eff=None, actual_brake_horsepower=None, theoretical_adiabatic_horsepower=None, **kwargs):
        return

    def eqn_8_5__Eff(actual_brake_horsepower: float, theoretical_adiabatic_horsepower: float, **kwargs):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        Eff = theoretical_adiabatic_horsepower/actual_brake_horsepower
        result.append(Eff)
        return result
    def eqn_8_5__actual_brake_horsepower(Eff: float, theoretical_adiabatic_horsepower: float, **kwargs):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        actual_brake_horsepower = theoretical_adiabatic_horsepower/Eff
        result.append(actual_brake_horsepower)
        return result
    def eqn_8_5__theoretical_adiabatic_horsepower(Eff: float, actual_brake_horsepower: float, **kwargs):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        theoretical_adiabatic_horsepower = Eff*actual_brake_horsepower
        result.append(theoretical_adiabatic_horsepower)
        return result
class SteamJetInjectors:
    @kwasak_static
    def eqn_9_1(A=None, rho_s=None, v=None, w_s=None, **kwargs):
        return

    def eqn_9_1__A(rho_s: float, v: float, w_s: float, **kwargs):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        A = w_s/(rho_s*v)
        result.append(A)
        return result
    def eqn_9_1__rho_s(A: float, v: float, w_s: float, **kwargs):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        rho_s = w_s/(A*v)
        result.append(rho_s)
        return result
    def eqn_9_1__v(A: float, rho_s: float, w_s: float, **kwargs):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        v = w_s/(A*rho_s)
        result.append(v)
        return result
    def eqn_9_1__w_s(A: float, rho_s: float, v: float, **kwargs):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        w_s = A*rho_s*v
        result.append(w_s)
        return result
    @kwasak_static
    def eqn_9_3(P_s=None, V=None, t_e=None, w_j=None, **kwargs):
        return

    def eqn_9_3__P_s(V: float, t_e: float, w_j: float, **kwargs):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        P_s = 33.3333333333333*(23.0*V - 10.0*t_e*w_j)/V
        result.append(P_s)
        return result
    def eqn_9_3__V(P_s: float, t_e: float, w_j: float, **kwargs):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        V = -1000.0*t_e*w_j/(3.0*P_s - 2300.0)
        result.append(V)
        return result
    def eqn_9_3__t_e(P_s: float, V: float, w_j: float, **kwargs):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        t_e = 0.001*V*(2300.0 - 3.0*P_s)/w_j
        result.append(t_e)
        return result
    def eqn_9_3__w_j(P_s: float, V: float, t_e: float, **kwargs):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        w_j = 0.001*V*(2300.0 - 3.0*P_s)/t_e
        result.append(w_j)
        return result
    @kwasak_static
    def eqn_9_4(AEL=None, SC=None, r=None, w_s=None, **kwargs):
        return

    def eqn_9_4__AEL(SC: float, r: float, w_s: float, **kwargs):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        AEL = w_s/(SC*r)
        result.append(AEL)
        return result
    def eqn_9_4__SC(AEL: float, r: float, w_s: float, **kwargs):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        SC = w_s/(AEL*r)
        result.append(SC)
        return result
    def eqn_9_4__r(AEL: float, SC: float, w_s: float, **kwargs):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        r = w_s/(AEL*SC)
        result.append(r)
        return result
    def eqn_9_4__w_s(AEL: float, SC: float, r: float, **kwargs):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        w_s = AEL*SC*r
        result.append(w_s)
        return result
    @kwasak_static
    def eqn_9_5(V=None, r_h=None, t_h=None, w_h=None, **kwargs):
        return

    def eqn_9_5__V(r_h: float, t_h: float, w_h: float, **kwargs):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        V = t_h*w_h/r_h
        result.append(V)
        return result
    def eqn_9_5__r_h(V: float, t_h: float, w_h: float, **kwargs):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        r_h = t_h*w_h/V
        result.append(r_h)
        return result
    def eqn_9_5__t_h(V: float, r_h: float, w_h: float, **kwargs):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        t_h = V*r_h/w_h
        result.append(t_h)
        return result
    def eqn_9_5__w_h(V: float, r_h: float, t_h: float, **kwargs):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        w_h = V*r_h/t_h
        result.append(w_h)
        return result
class VacuumTheory:
    @kwasak_static
    def eqn_1_10(P_1=None, P_2=None, T_1=None, T_2=None, V_1=None, V_2=None, **kwargs):
        return

    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float, **kwargs):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return result
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float, **kwargs):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return result
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float, **kwargs):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return result
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float, **kwargs):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return result
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float, **kwargs):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return result
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float, **kwargs):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return result
    @kwasak_static
    def eqn_1_11(M=None, P=None, T=None, W=None, q=None, **kwargs):
        return

    def eqn_1_11__M(P: float, T: float, W: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result
    def eqn_1_11__P(M: float, T: float, W: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return result
    def eqn_1_11__T(M: float, P: float, W: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return result
    def eqn_1_11__W(M: float, P: float, T: float, q: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return result
    def eqn_1_11__q(M: float, P: float, T: float, W: float, **kwargs):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result
    @kwasak_static
    def eqn_1_12(Total_P=None, sum_partial_pressures=None, **kwargs):
        return

    def eqn_1_12__Total_P(sum_partial_pressures: float, **kwargs):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return result
    def eqn_1_12__sum_partial_pressures(Total_P: float, **kwargs):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return result
    @kwasak_static
    def eqn_1_13a(n=None, n_a=None, y_a=None, **kwargs):
        return

    def eqn_1_13a__n(n_a: float, y_a: float, **kwargs):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result
    def eqn_1_13a__n_a(n: float, y_a: float, **kwargs):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result
    def eqn_1_13a__y_a(n: float, n_a: float, **kwargs):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result
    @kwasak_static
    def eqn_1_13b(P=None, p_a=None, y_a=None, **kwargs):
        return

    def eqn_1_13b__P(p_a: float, y_a: float, **kwargs):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return result
    def eqn_1_13b__p_a(P: float, y_a: float, **kwargs):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return result
    def eqn_1_13b__y_a(P: float, p_a: float, **kwargs):
        # [.pyeqn] y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return result
    @kwasak_static
    def eqn_1_7(R=None, T=None, V=None, n=None, p=None, **kwargs):
        return

    def eqn_1_7__R(T: float, V: float, n: float, p: float, **kwargs):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return result
    def eqn_1_7__T(R: float, V: float, n: float, p: float, **kwargs):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return result
    def eqn_1_7__V(R: float, T: float, n: float, p: float, **kwargs):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return result
    def eqn_1_7__n(R: float, T: float, V: float, p: float, **kwargs):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return result
    def eqn_1_7__p(R: float, T: float, V: float, n: float, **kwargs):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return result
    @kwasak_static
    def eqn_1_8(M=None, P=None, R=None, T=None, V=None, m=None, **kwargs):
        return

    def eqn_1_8__M(P: float, R: float, T: float, V: float, m: float, **kwargs):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return result
    def eqn_1_8__P(M: float, R: float, T: float, V: float, m: float, **kwargs):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return result
    def eqn_1_8__R(M: float, P: float, T: float, V: float, m: float, **kwargs):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return result
    def eqn_1_8__T(M: float, P: float, R: float, V: float, m: float, **kwargs):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return result
    def eqn_1_8__V(M: float, P: float, R: float, T: float, m: float, **kwargs):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return result
    def eqn_1_8__m(M: float, P: float, R: float, T: float, V: float, **kwargs):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return result
    @kwasak_static
    def eqn_1_9(M=None, P=None, R=None, T=None, rho=None, **kwargs):
        return

    def eqn_1_9__M(P: float, R: float, T: float, rho: float, **kwargs):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result
    def eqn_1_9__P(M: float, R: float, T: float, rho: float, **kwargs):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return result
    def eqn_1_9__R(M: float, P: float, T: float, rho: float, **kwargs):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result
    def eqn_1_9__T(M: float, P: float, R: float, rho: float, **kwargs):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return result
    def eqn_1_9__rho(M: float, P: float, R: float, T: float, **kwargs):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return result
