from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class AirLeak:
    @kwasak
    def eqn_4_10(self, T=None, V=None, del_P=None, leakage=None, t=None):
        return
    def eqn_4_10__T(self, V: float, del_P: float, leakage: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        T = 3.127*V*del_P/(leakage*t)
        result.append(T)
        return result
    def eqn_4_10__V(self, T: float, del_P: float, leakage: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        V = 0.319795330988168*T*leakage*t/del_P
        result.append(V)
        return result
    def eqn_4_10__del_P(self, T: float, V: float, leakage: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        del_P = 0.319795330988168*T*leakage*t/V
        result.append(del_P)
        return result
    def eqn_4_10__leakage(self, T: float, V: float, del_P: float, t: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        leakage = 3.127*V*del_P/(T*t)
        result.append(leakage)
        return result
    def eqn_4_10__t(self, T: float, V: float, del_P: float, leakage: float, **kwargs):
        # leakage = 0.0059 * V * del_P / t * 530 / T  # lb/hr
        result = []
        t = 3.127*V*del_P/(T*leakage)
        result.append(t)
        return result
    @kwasak
    def eqn_4_7(self, W=None, W_T=None, sum_individual_leak_rates=None):
        return
    def eqn_4_7__W(self, W_T: float, sum_individual_leak_rates: float, **kwargs):
        # W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return result
    def eqn_4_7__W_T(self, W: float, sum_individual_leak_rates: float, **kwargs):
        # W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return result
    def eqn_4_7__sum_individual_leak_rates(self, W: float, W_T: float, **kwargs):
        # W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return result

class FluidFlowVacuumLines:
    @kwasak
    def eqn_2_1(self, D=None, Re=None, mu=None, rho=None, v=None):
        """
        rho := density, lb/ft^3
        D := pipe inside diam, ft
        v := vel. ft/s
        mu := viscosity, lb/ft*s
        """
        return
    def eqn_2_1__D(self, Re: float, mu: float, rho: float, v: float, **kwargs):
        # Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return result
    def eqn_2_1__Re(self, D: float, mu: float, rho: float, v: float, **kwargs):
        # Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return result
    def eqn_2_1__mu(self, D: float, Re: float, rho: float, v: float, **kwargs):
        # Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return result
    def eqn_2_1__rho(self, D: float, Re: float, mu: float, v: float, **kwargs):
        # Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return result
    def eqn_2_1__v(self, D: float, Re: float, mu: float, rho: float, **kwargs):
        # Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return result
    @kwasak
    def eqn_2_10(self, Suc_Pres=None, delta_P=None, oper_press=None):
        """
        delta_P := pressure loss
        """
        return
    def eqn_2_10__Suc_Pres(self, delta_P: float, oper_press: float, **kwargs):
        # Suc_Pres = oper_press - delta_P
        result = []
        Suc_Pres = -delta_P + oper_press
        result.append(Suc_Pres)
        return result
    def eqn_2_10__delta_P(self, Suc_Pres: float, oper_press: float, **kwargs):
        # Suc_Pres = oper_press - delta_P
        result = []
        delta_P = -Suc_Pres + oper_press
        result.append(delta_P)
        return result
    def eqn_2_10__oper_press(self, Suc_Pres: float, delta_P: float, **kwargs):
        # Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return result
    @kwasak
    def eqn_2_11(self, D=None, L=None, f=None, g_c=None, h_r=None, v=None):
        """
        f:= Moody friction
        L:=length_pipe, ft
        v:= velocity, ft/s
        D:= inside diameter, ft
        g_c:= dimensional constant, 32.2 lb * ft / lb * s
        """
        return
    def eqn_2_11__D(self, L: float, f: float, g_c: float, h_r: float, v: float, **kwargs):
        # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
        return result
    def eqn_2_11__L(self, D: float, f: float, g_c: float, h_r: float, v: float, **kwargs):
        # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return result
    def eqn_2_11__f(self, D: float, L: float, g_c: float, h_r: float, v: float, **kwargs):
        # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return result
    def eqn_2_11__g_c(self, D: float, L: float, f: float, h_r: float, v: float, **kwargs):
        # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        g_c = L*f*v**2/(2*D*h_r)
        result.append(g_c)
        return result
    def eqn_2_11__h_r(self, D: float, L: float, f: float, g_c: float, v: float, **kwargs):
        # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        h_r = L*f*v**2/(2*D*g_c)
        result.append(h_r)
        return result
    def eqn_2_11__v(self, D: float, L: float, f: float, g_c: float, h_r: float, **kwargs):
        # h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return result
    @kwasak
    def eqn_2_12(self, L=None, d=None, delta_P=None, f=None, g=None, rho=None, v=None):
        """
        rho:= density, lb/ft^3
        d:= pipe inside diameter, in
        q:= vol. flow rate, ft^3/min
        """
        return
    def eqn_2_12__L(self, d: float, delta_P: float, f: float, g: float, rho: float, v: float, **kwargs):
        # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838*d*delta_P*g/(f*rho*v**2)
        result.append(L)
        return result
    def eqn_2_12__d(self, L: float, delta_P: float, f: float, g: float, rho: float, v: float, **kwargs):
        # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155*L*f*rho*v**2/(delta_P*g)
        result.append(d)
        return result
    def eqn_2_12__delta_P(self, L: float, d: float, f: float, g: float, rho: float, v: float, **kwargs):
        # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        delta_P = 2.155*L*f*rho*v**2/(d*g)
        result.append(delta_P)
        return result
    def eqn_2_12__f(self, L: float, d: float, delta_P: float, g: float, rho: float, v: float, **kwargs):
        # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        f = 0.464037122969838*d*delta_P*g/(L*rho*v**2)
        result.append(f)
        return result
    def eqn_2_12__g(self, L: float, d: float, delta_P: float, f: float, rho: float, v: float, **kwargs):
        # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155*L*f*rho*v**2/(d*delta_P)
        result.append(g)
        return result
    def eqn_2_12__rho(self, L: float, d: float, delta_P: float, f: float, g: float, v: float, **kwargs):
        # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838*d*delta_P*g/(L*f*v**2)
        result.append(rho)
        return result
    def eqn_2_12__v(self, L: float, d: float, delta_P: float, f: float, g: float, rho: float, **kwargs):
        # delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        v = -0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        v = 0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        return result
    @kwasak
    def eqn_2_13(self, L=None, d=None, delta_P=None, f=None, q=None, rho=None):
        """
        rho:= density, lb/ft^3
        d:= pipe inside diameter, in
        q:= vol. flow rate, ft^3/min
        """
        return
    def eqn_2_13__L(self, d: float, delta_P: float, f: float, q: float, rho: float, **kwargs):
        # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
        return result
    def eqn_2_13__d(self, L: float, delta_P: float, f: float, q: float, rho: float, **kwargs):
        # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
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
    def eqn_2_13__delta_P(self, L: float, d: float, f: float, q: float, rho: float, **kwargs):
        # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        delta_P = 2.15*L*f*q**2*rho/d**5
        result.append(delta_P)
        return result
    def eqn_2_13__f(self, L: float, d: float, delta_P: float, q: float, rho: float, **kwargs):
        # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        f = 0.465116279069767*d**5*delta_P/(L*q**2*rho)
        result.append(f)
        return result
    def eqn_2_13__q(self, L: float, d: float, delta_P: float, f: float, rho: float, **kwargs):
        # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        return result
    def eqn_2_13__rho(self, L: float, d: float, delta_P: float, f: float, q: float, **kwargs):
        # delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return result
    @kwasak
    def eqn_2_14(self, M=None, R=None, T=None, g_c=None, k=None, v_s=None):
        """
        v_s := sonic_velocity
        k:=ratio of specific heat at constant temp to the specific heat at constant volume
        """
        return
    def eqn_2_14__M(self, R: float, T: float, g_c: float, k: float, v_s: float, **kwargs):
        # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return result
    def eqn_2_14__R(self, M: float, T: float, g_c: float, k: float, v_s: float, **kwargs):
        # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return result
    def eqn_2_14__T(self, M: float, R: float, g_c: float, k: float, v_s: float, **kwargs):
        # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return result
    def eqn_2_14__g_c(self, M: float, R: float, T: float, k: float, v_s: float, **kwargs):
        # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return result
    def eqn_2_14__k(self, M: float, R: float, T: float, g_c: float, v_s: float, **kwargs):
        # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return result
    def eqn_2_14__v_s(self, M: float, R: float, T: float, g_c: float, k: float, **kwargs):
        # v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return result
    @kwasak
    def eqn_2_15(self, Re=None, f=None):
        return
    def eqn_2_15__Re(self, f: float, **kwargs):
        # f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736/f**4
        result.append(Re)
        return result
    def eqn_2_15__f(self, Re: float, **kwargs):
        # f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316/Re**(1/4)
        result.append(f)
        return result
    @kwasak
    def eqn_2_16(self, Re=None, f=None):
        return
    def eqn_2_16__Re(self, f: float, **kwargs):
        # f = 64 / Re
        result = []
        Re = 64/f
        result.append(Re)
        return result
    def eqn_2_16__f(self, Re: float, **kwargs):
        # f = 64 / Re
        result = []
        f = 64/Re
        result.append(f)
        return result
    @kwasak
    def eqn_2_17(self, L=None, d=None, delta_P=None, mu=None, v=None):
        return
    def eqn_2_17__L(self, d: float, delta_P: float, mu: float, v: float, **kwargs):
        # delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
        return result
    def eqn_2_17__d(self, L: float, delta_P: float, mu: float, v: float, **kwargs):
        # delta_P = 0.0345* mu * L * v / d**2
        result = []
        d = -0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        d = 0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        return result
    def eqn_2_17__delta_P(self, L: float, d: float, mu: float, v: float, **kwargs):
        # delta_P = 0.0345* mu * L * v / d**2
        result = []
        delta_P = 0.0345*L*mu*v/d**2
        result.append(delta_P)
        return result
    def eqn_2_17__mu(self, L: float, d: float, delta_P: float, v: float, **kwargs):
        # delta_P = 0.0345* mu * L * v / d**2
        result = []
        mu = 28.9855072463768*d**2*delta_P/(L*v)
        result.append(mu)
        return result
    def eqn_2_17__q(self, L: float, d: float, delta_P: float, mu: float, **kwargs):
        # delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952*d**4*delta_P/(L*mu)
        result.append(q)
        return result
    def eqn_2_17__v(self, L: float, d: float, delta_P: float, mu: float, **kwargs):
        # delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
        return result
    @kwasak
    def eqn_2_18a(self, D_eq=None, R_ll=None):
        return
    def eqn_2_18a__D_eq(self, R_ll: float, **kwargs):
        # D_eq = 4 * R_ll
        result = []
        D_eq = 4*R_ll
        result.append(D_eq)
        return result
    def eqn_2_18a__R_ll(self, D_eq: float, **kwargs):
        # D_eq = 4 * R_ll
        result = []
        R_ll = D_eq/4
        result.append(R_ll)
        return result
    @kwasak
    def eqn_2_18b(self, R_ll=None, h=None, w=None):
        return
    def eqn_2_18b__R_ll(self, h: float, w: float, **kwargs):
        # R_ll = w * h / (2 * (w + h))
        result = []
        R_ll = h*w/(2*(h + w))
        result.append(R_ll)
        return result
    def eqn_2_18b__h(self, R_ll: float, w: float, **kwargs):
        # R_ll = w * h / (2 * (w + h))
        result = []
        h = 2*R_ll*w/(-2*R_ll + w)
        result.append(h)
        return result
    def eqn_2_18b__w(self, R_ll: float, h: float, **kwargs):
        # R_ll = w * h / (2 * (w + h))
        result = []
        w = 2*R_ll*h/(-2*R_ll + h)
        result.append(w)
        return result
    @kwasak
    def eqn_2_19a(self, R_ll=None, Re=None, mu=None, rho=None, v=None):
        return
    def eqn_2_19a__R_ll(self, Re: float, mu: float, rho: float, v: float, **kwargs):
        # Re = 4 * R_ll * rho * v / mu
        result = []
        R_ll = Re*mu/(4*rho*v)
        result.append(R_ll)
        return result
    def eqn_2_19a__Re(self, R_ll: float, mu: float, rho: float, v: float, **kwargs):
        # Re = 4 * R_ll * rho * v / mu
        result = []
        Re = 4*R_ll*rho*v/mu
        result.append(Re)
        return result
    def eqn_2_19a__mu(self, R_ll: float, Re: float, rho: float, v: float, **kwargs):
        # Re = 4 * R_ll * rho * v / mu
        result = []
        mu = 4*R_ll*rho*v/Re
        result.append(mu)
        return result
    def eqn_2_19a__rho(self, R_ll: float, Re: float, mu: float, v: float, **kwargs):
        # Re = 4 * R_ll * rho * v / mu
        result = []
        rho = Re*mu/(4*R_ll*v)
        result.append(rho)
        return result
    def eqn_2_19a__v(self, R_ll: float, Re: float, mu: float, rho: float, **kwargs):
        # Re = 4 * R_ll * rho * v / mu
        result = []
        v = Re*mu/(4*R_ll*rho)
        result.append(v)
        return result
    @kwasak
    def eqn_2_19b(self, Re=None, h=None, mu=None, rho=None, v=None, w=None):
        return
    def eqn_2_19b__Re(self, h: float, mu: float, rho: float, v: float, w: float, **kwargs):
        # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        Re = 2*h*rho*v*w/(mu*(h + w))
        result.append(Re)
        return result
    def eqn_2_19b__h(self, Re: float, mu: float, rho: float, v: float, w: float, **kwargs):
        # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        h = Re*mu*w/(-Re*mu + 2*rho*v*w)
        result.append(h)
        return result
    def eqn_2_19b__mu(self, Re: float, h: float, rho: float, v: float, w: float, **kwargs):
        # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        mu = 2*h*rho*v*w/(Re*(h + w))
        result.append(mu)
        return result
    def eqn_2_19b__rho(self, Re: float, h: float, mu: float, v: float, w: float, **kwargs):
        # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        rho = Re*mu*(h + w)/(2*h*v*w)
        result.append(rho)
        return result
    def eqn_2_19b__v(self, Re: float, h: float, mu: float, rho: float, w: float, **kwargs):
        # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        v = Re*mu*(h + w)/(2*h*rho*w)
        result.append(v)
        return result
    def eqn_2_19b__w(self, Re: float, h: float, mu: float, rho: float, v: float, **kwargs):
        # Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        w = Re*h*mu/(-Re*mu + 2*h*rho*v)
        result.append(w)
        return result
    @kwasak
    def eqn_2_2(self, delta=None, lambd=None, psi=None):
        """
        lambd := average mean free path , in
        delta := mol. diam , in
        psi:= mol. density molecules/in^3
        """
        return
    def eqn_2_2__delta(self, lambd: float, psi: float, **kwargs):
        # lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        delta = -0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        delta = 0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        return result
    def eqn_2_2__lambd(self, delta: float, psi: float, **kwargs):
        # lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        lambd = 4.44288293815837*delta**2*psi
        result.append(lambd)
        return result
    def eqn_2_2__psi(self, delta: float, lambd: float, **kwargs):
        # lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        psi = 0.225079079039277*lambd/delta**2
        result.append(psi)
        return result
    @kwasak
    def eqn_2_20(self, L=None, sum_equivalent_length=None, sum_pipe=None):
        """
        L:= laminar flow
        """
        return
    def eqn_2_20__L(self, sum_equivalent_length: float, sum_pipe: float, **kwargs):
        # L = sum_pipe + sum_equivalent_length
        result = []
        L = sum_equivalent_length + sum_pipe
        result.append(L)
        return result
    def eqn_2_20__sum_equivalent_length(self, L: float, sum_pipe: float, **kwargs):
        # L = sum_pipe + sum_equivalent_length
        result = []
        sum_equivalent_length = L - sum_pipe
        result.append(sum_equivalent_length)
        return result
    def eqn_2_20__sum_pipe(self, L: float, sum_equivalent_length: float, **kwargs):
        # L = sum_pipe + sum_equivalent_length
        result = []
        sum_pipe = L - sum_equivalent_length
        result.append(sum_pipe)
        return result
    @kwasak
    def eqn_2_22(self, P_s=None, Q_throughput=None, S_p=None):
        """
        Q:= through_put, sucking pressure P
        S_p:= dV / Dt
        """
        return
    def eqn_2_22__P_s(self, Q_throughput: float, S_p: float, **kwargs):
        # Q_throughput = S_p * P_s
        result = []
        P_s = Q_throughput/S_p
        result.append(P_s)
        return result
    def eqn_2_22__Q_throughput(self, P_s: float, S_p: float, **kwargs):
        # Q_throughput = S_p * P_s
        result = []
        Q_throughput = P_s*S_p
        result.append(Q_throughput)
        return result
    def eqn_2_22__S_p(self, P_s: float, Q_throughput: float, **kwargs):
        # Q_throughput = S_p * P_s
        result = []
        S_p = Q_throughput/P_s
        result.append(S_p)
        return result
    @kwasak
    def eqn_2_25(self, C=None, P_1=None, P_2=None, Q_throughput=None):
        return
    def eqn_2_25__C(self, P_1: float, P_2: float, Q_throughput: float, **kwargs):
        # C = Q_throughput / (P_1 - P_2)
        result = []
        C = Q_throughput/(P_1 - P_2)
        result.append(C)
        return result
    def eqn_2_25__P_1(self, C: float, P_2: float, Q_throughput: float, **kwargs):
        # C = Q_throughput / (P_1 - P_2)
        result = []
        P_1 = P_2 + Q_throughput/C
        result.append(P_1)
        return result
    def eqn_2_25__P_2(self, C: float, P_1: float, Q_throughput: float, **kwargs):
        # C = Q_throughput / (P_1 - P_2)
        result = []
        P_2 = P_1 - Q_throughput/C
        result.append(P_2)
        return result
    def eqn_2_25__Q_throughput(self, C: float, P_1: float, P_2: float, **kwargs):
        # C = Q_throughput / (P_1 - P_2)
        result = []
        Q_throughput = C*(P_1 - P_2)
        result.append(Q_throughput)
        return result
    @kwasak
    def eqn_2_26(self, D=None, L=None, P_downstream=None, P_p=None, P_upstream=None, mu=None, q=None):
        return
    def eqn_2_26__D(self, L: float, P_downstream: float, P_p: float, P_upstream: float, mu: float, q: float, **kwargs):
        # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        D = -14953.4878122122*I*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = 14953.4878122122*I*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = -14953.4878122122*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = 14953.4878122122*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        return result
    def eqn_2_26__L(self, D: float, P_downstream: float, P_p: float, P_upstream: float, mu: float, q: float, **kwargs):
        # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        L = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(mu*q)
        result.append(L)
        return result
    def eqn_2_26__P_downstream(self, D: float, L: float, P_p: float, P_upstream: float, mu: float, q: float, **kwargs):
        # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_downstream = P_upstream - 40.7436654315252*L*mu*q/D**4
        result.append(P_downstream)
        return result
    def eqn_2_26__P_p(self, D: float, L: float, P_downstream: float, P_upstream: float, mu: float, q: float, **kwargs):
        # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_p = 0.0
        result.append(P_p)
        return result
    def eqn_2_26__P_upstream(self, D: float, L: float, P_downstream: float, P_p: float, mu: float, q: float, **kwargs):
        # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_upstream = P_downstream + 40.7436654315252*L*mu*q/D**4
        result.append(P_upstream)
        return result
    def eqn_2_26__mu(self, D: float, L: float, P_downstream: float, P_p: float, P_upstream: float, q: float, **kwargs):
        # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        mu = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(L*q)
        result.append(mu)
        return result
    def eqn_2_26__q(self, D: float, L: float, P_downstream: float, P_p: float, P_upstream: float, mu: float, **kwargs):
        # q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        q = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(L*mu)
        result.append(q)
        return result
    @kwasak
    def eqn_2_28(self, C=None, D=None, L=None, P_p=None, mu=None):
        return
    def eqn_2_28__C(self, D: float, L: float, P_p: float, mu: float, **kwargs):
        # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        C = 0.0245436926061703*D**4*P_p/(L*mu)
        result.append(C)
        return result
    def eqn_2_28__D(self, C: float, L: float, P_p: float, mu: float, **kwargs):
        # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        D = -2.52647511098426*I*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = -2.52647511098426*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = 2.52647511098426*(C*L*mu/P_p)**(1/4)
        result.append(D)
        return result
    def eqn_2_28__L(self, C: float, D: float, P_p: float, mu: float, **kwargs):
        # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        L = 0.0245436926061703*D**4*P_p/(C*mu)
        result.append(L)
        return result
    def eqn_2_28__P_p(self, C: float, D: float, L: float, mu: float, **kwargs):
        # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        P_p = 40.7436654315252*C*L*mu/D**4
        result.append(P_p)
        return result
    def eqn_2_28__mu(self, C: float, D: float, L: float, P_p: float, **kwargs):
        # C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        mu = 0.0245436926061703*D**4*P_p/(C*L)
        result.append(mu)
        return result
    @kwasak
    def eqn_2_29(self, C=None, S_1=None, S_2=None):
        return
    def eqn_2_29__C(self, S_1: float, S_2: float, **kwargs):
        # S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        C = -S_1*S_2/(S_1 - S_2)
        result.append(C)
        return result
    def eqn_2_29__S_1(self, C: float, S_2: float, **kwargs):
        # S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_1 = C*S_2/(C + S_2)
        result.append(S_1)
        return result
    def eqn_2_29__S_2(self, C: float, S_1: float, **kwargs):
        # S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_2 = C*S_1/(C - S_1)
        result.append(S_2)
        return result
    @kwasak
    def eqn_2_3(self, D=None, kn=None, lambd=None):
        """
        D:= inside diameter, in
        lambd:=avg. mean free path, in
        """
        return
    def eqn_2_3__D(self, kn: float, lambd: float, **kwargs):
        # kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return result
    def eqn_2_3__kn(self, D: float, lambd: float, **kwargs):
        # kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return result
    def eqn_2_3__lambd(self, D: float, kn: float, **kwargs):
        # kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return result
    @kwasak
    def eqn_2_31(self, C=None, S_p=None, S_pump_speed=None):
        return
    def eqn_2_31__C(self, S_p: float, S_pump_speed: float, **kwargs):
        # S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        C = S_p*S_pump_speed/(S_p - S_pump_speed)
        result.append(C)
        return result
    def eqn_2_31__S_p(self, C: float, S_pump_speed: float, **kwargs):
        # S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_p = C*S_pump_speed/(C - S_pump_speed)
        result.append(S_p)
        return result
    def eqn_2_31__S_pump_speed(self, C: float, S_p: float, **kwargs):
        # S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_pump_speed = C*S_p/(C + S_p)
        result.append(S_pump_speed)
        return result
    @kwasak
    def eqn_2_32(self, C_series=None, geometric_sum_C=None):
        return
    def eqn_2_32__C_series(self, geometric_sum_C: float, **kwargs):
        # 1 / C_series = geometric_sum_C
        result = []
        C_series = 1/geometric_sum_C
        result.append(C_series)
        return result
    def eqn_2_32__geometric_sum_C(self, C_series: float, **kwargs):
        # 1 / C_series = geometric_sum_C
        result = []
        geometric_sum_C = 1/C_series
        result.append(geometric_sum_C)
        return result
    @kwasak
    def eqn_2_33(self, C_paralell=None, arithmetic_sum_C=None):
        return
    def eqn_2_33__C_paralell(self, arithmetic_sum_C: float, **kwargs):
        # 1 / C_paralell = arithmetic_sum_C
        result = []
        C_paralell = 1/arithmetic_sum_C
        result.append(C_paralell)
        return result
    def eqn_2_33__arithmetic_sum_C(self, C_paralell: float, **kwargs):
        # 1 / C_paralell = arithmetic_sum_C
        result = []
        arithmetic_sum_C = 1/C_paralell
        result.append(arithmetic_sum_C)
        return result
    @kwasak
    def eqn_2_34(self, C=None, C_1=None, C_2=None, D=None, L=None, P_p=None, mu=None):
        return
    def eqn_2_34__C(self, C_1: float, C_2: float, D: float, L: float, P_p: float, mu: float, **kwargs):
        # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        # Solve for C by rearranging the equation
        result = []
        try:
            C = (mu * L) / ((D**4/P_p)*C_1 + D**3/L*C_2)
            result.append(float(C))  # Ensure that we return a float value, not complex numbers if any roots are imaginary
        except ZeroDivisionError:
            print("Error: Division by zero encountered")
        else:
            return [result]
    def eqn_2_34__C_1(self, C: float, C_2: float, D: float, L: float, P_p: float, mu: float, **kwargs):
        # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C_1 = mu*(C*L - C_2*D**3)/(D**4*P_p)
        result.append(C_1)
        return result
    def eqn_2_34__C_2(self, C: float, C_1: float, D: float, L: float, P_p: float, mu: float, **kwargs):
        # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C_2 = C*L/D**3 - C_1*D*P_p/mu
        result.append(C_2)
        return result
    def eqn_2_34__D(self, C: float, C_1: float, C_2: float, L: float, P_p: float, mu: float, **kwargs):
        # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        D = Piecewise((-sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (-sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((-sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (-sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        return result
    def eqn_2_34__L(self, C: float, C_1: float, C_2: float, D: float, P_p: float, mu: float, **kwargs):
        # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        L = D**3*(C_1*D*P_p + C_2*mu)/(C*mu)
        result.append(L)
        return result
    def eqn_2_34__P_p(self, C: float, C_1: float, C_2: float, D: float, L: float, mu: float, **kwargs):
        # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        P_p = mu*(C*L - C_2*D**3)/(C_1*D**4)
        result.append(P_p)
        return result
    def eqn_2_34__mu(self, C: float, C_1: float, C_2: float, D: float, L: float, P_p: float, **kwargs):
        # C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        mu = C_1*D**4*P_p/(C*L - C_2*D**3)
        result.append(mu)
        return result
    @kwasak
    def eqn_2_35(self, C_L=None, C_T=None, F_p=None):
        """
        F_P:= correction factor for Poiseuille's eqn from Figure
        """
        return
    def eqn_2_35__C_L(self, C_T: float, F_p: float, **kwargs):
        # C_T = C_L * F_p
        result = []
        C_L = C_T/F_p
        result.append(C_L)
        return result
    def eqn_2_35__C_T(self, C_L: float, F_p: float, **kwargs):
        # C_T = C_L * F_p
        result = []
        C_T = C_L*F_p
        result.append(C_T)
        return result
    def eqn_2_35__F_p(self, C_L: float, C_T: float, **kwargs):
        # C_T = C_L * F_p
        result = []
        F_p = C_T/C_L
        result.append(F_p)
        return result
    @kwasak
    def eqn_2_36(self, C=None, C_0=None, F_t=None):
        """
        C_0:=conductance thin walled aperture
        F_t:=transmission prob. for component
        """
        return
    def eqn_2_36__C(self, C_0: float, F_t: float, **kwargs):
        # C = C_0 * F_t
        result = []
        C = C_0*F_t
        result.append(C)
        return result
    def eqn_2_36__C_0(self, C: float, F_t: float, **kwargs):
        # C = C_0 * F_t
        result = []
        C_0 = C/F_t
        result.append(C_0)
        return result
    def eqn_2_36__F_t(self, C: float, C_0: float, **kwargs):
        # C = C_0 * F_t
        result = []
        F_t = C/C_0
        result.append(F_t)
        return result
    @kwasak
    def eqn_2_37(self, A=None, C=None, F_t=None, M=None, T=None):
        """
        F_t:= 1, for an aperture
        """
        return
    def eqn_2_37__A(self, C: float, F_t: float, M: float, T: float, **kwargs):
        # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        A = 0.000681714375311032*C**2*M/(F_t*T)
        result.append(A)
        return result
    def eqn_2_37__C(self, A: float, F_t: float, M: float, T: float, **kwargs):
        # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        C = 38.3*sqrt(A*F_t*T/M)
        result.append(C)
        return result
    def eqn_2_37__F_t(self, A: float, C: float, M: float, T: float, **kwargs):
        # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        F_t = 0.000681714375311032*C**2*M/(A*T)
        result.append(F_t)
        return result
    def eqn_2_37__M(self, A: float, C: float, F_t: float, T: float, **kwargs):
        # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        M = 1466.89*A*F_t*T/C**2
        result.append(M)
        return result
    def eqn_2_37__T(self, A: float, C: float, F_t: float, M: float, **kwargs):
        # C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        T = 0.000681714375311032*C**2*M/(A*F_t)
        result.append(T)
        return result
    @kwasak
    def eqn_2_4(self, _beta=None, mu=None, vel_grad=None):
        """
        mu:=coefficient of viscosity
        """
        return
    def eqn_2_4___beta(self, mu: float, vel_grad: float, **kwargs):
        # _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return result
    def eqn_2_4__mu(self, _beta: float, vel_grad: float, **kwargs):
        # _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return result
    def eqn_2_4__vel_grad(self, _beta: float, mu: float, **kwargs):
        # _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return result
    @kwasak
    def eqn_2_5(self, D=None, L=None, delta_P=None, mu=None, q=None):
        """
        q:=volumetric flow cm^3/s
        D:= pipe diam.,cm
        delta_P := upstream-downstream pressure, dyne/cm^3
        L:=length, cm
        mu:= coef. of visco., poise
        """
        return
    def eqn_2_5__D(self, L: float, delta_P: float, mu: float, q: float, **kwargs):
        # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
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
    def eqn_2_5__L(self, D: float, delta_P: float, mu: float, q: float, **kwargs):
        # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703*D**4*delta_P/(mu*q)
        result.append(L)
        return result
    def eqn_2_5__delta_P(self, D: float, L: float, mu: float, q: float, **kwargs):
        # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252*L*mu*q/D**4
        result.append(delta_P)
        return result
    def eqn_2_5__mu(self, D: float, L: float, delta_P: float, q: float, **kwargs):
        # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703*D**4*delta_P/(L*q)
        result.append(mu)
        return result
    def eqn_2_5__q(self, D: float, L: float, delta_P: float, mu: float, **kwargs):
        # q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703*D**4*delta_P/(L*mu)
        result.append(q)
        return result
    @kwasak
    def eqn_2_6(self, lambd=None, mu=None, rho=None, v_a=None):
        """
        mu :=viscosity, poise
        rho:= density, g/cm^3
        lambd:= mean free path, cm
        """
        return
    def eqn_2_6__lambd(self, mu: float, rho: float, v_a: float, **kwargs):
        # mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return result
    def eqn_2_6__mu(self, lambd: float, rho: float, v_a: float, **kwargs):
        # mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return result
    def eqn_2_6__rho(self, lambd: float, mu: float, v_a: float, **kwargs):
        # mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return result
    def eqn_2_6__v_a(self, lambd: float, mu: float, rho: float, **kwargs):
        # mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return result
    @kwasak
    def eqn_2_7(self, T=None, k=None, m=None, v_a=None):
        """
        k:=boltz
        T:= abs temp
        m:= mass of a molecule
        """
        return
    def eqn_2_7__T(self, k: float, m: float, v_a: float, **kwargs):
        # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        T = 0.392699081698724*m*v_a**2/k
        result.append(T)
        return result
    def eqn_2_7__k(self, T: float, m: float, v_a: float, **kwargs):
        # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        k = 0.392699081698724*m*v_a**2/T
        result.append(k)
        return result
    def eqn_2_7__m(self, T: float, k: float, v_a: float, **kwargs):
        # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return result
    def eqn_2_7__v_a(self, T: float, k: float, m: float, **kwargs):
        # v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return result
    @kwasak
    def eqn_2_8(self, M=None, P_c=None, T_c=None, mu_c=None):
        """
        M:= mol. weight
        T_c:= critical temp, K
        P_c:= critical pressure, atm
        """
        return
    def eqn_2_8__M(self, P_c: float, T_c: float, mu_c: float, **kwargs):
        # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return result
    def eqn_2_8__P_c(self, M: float, T_c: float, mu_c: float, **kwargs):
        # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return result
    def eqn_2_8__T_c(self, M: float, P_c: float, mu_c: float, **kwargs):
        # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return result
    def eqn_2_8__mu_c(self, M: float, P_c: float, T_c: float, **kwargs):
        # mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return result

class LiquidRing:
    @kwasak
    def eqn_10_1(self, D_r=None, sig_R=None, w=None):
        """
        sig_R := rotor tip speed ft/s
        D_r := rotor Diameter
        w := rotational speed
        """
        return
    def eqn_10_1__D_r(self, sig_R: float, w: float, **kwargs):
        # sig_R = 0.00436 * D_r * w
        result = []
        D_r = 229.357798165138*sig_R/w
        result.append(D_r)
        return result
    def eqn_10_1__sig_R(self, D_r: float, w: float, **kwargs):
        # sig_R = 0.00436 * D_r * w
        result = []
        sig_R = 0.00436*D_r*w
        result.append(sig_R)
        return result
    def eqn_10_1__w(self, D_r: float, sig_R: float, **kwargs):
        # sig_R = 0.00436 * D_r * w
        result = []
        w = 229.357798165138*sig_R/D_r
        result.append(w)
        return result
    @kwasak
    def eqn_10_10(self, bhp=None, bhp_0=None, mu=None, rho=None):
        return
    def eqn_10_10__bhp(self, bhp_0: float, mu: float, rho: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp = 0.0005*bhp_0*(31.0*mu**(4/25)*rho**(21/25) + 1000.0)
        result.append(bhp)
        return result
    def eqn_10_10__bhp_0(self, bhp: float, mu: float, rho: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp_0 = 2000.0*bhp/(31.0*mu**0.16*rho**0.84 + 1000.0)
        result.append(bhp_0)
        return result
    def eqn_10_10__mu(self, bhp: float, bhp_0: float, rho: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        mu = -4.7751763343393e-10*I*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        mu = 4.7751763343393e-10*I*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        mu = -4.7751763343393e-10*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        mu = 4.7751763343393e-10*(2000.0*bhp/(bhp_0*rho**0.84) - 1000.0/rho**0.84)**(25/4)
        result.append(mu)
        return result
    def eqn_10_10__rho(self, bhp: float, bhp_0: float, mu: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        # Solve for rho:
        # Step 1: bhp / bhp_0 = 0.5 + 0.0155 *rho ** 0.84* mu ** 0.16
        # Step 2: (bhp / bhp_0 - 0.5) = 0.0155 * rho ** 0.84 * mu ** 0.16
        # Step 3: rho ** 0.84 = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16))
        # Step 4: rho = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16)) ** (1.0 / 0.84)
        rho = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16)) ** (1.0 / 0.84)
        return [rho]
    @kwasak
    def eqn_10_11(self, T_c=None, T_s=None):
        return
    def eqn_10_11__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 10
        result = []
        T_c = T_s + 10
        result.append(T_c)
        return result
    def eqn_10_11__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 10
        result = []
        T_s = T_c - 10
        result.append(T_s)
        return result
    @kwasak
    def eqn_10_12(self, T_c=None, T_s=None):
        return
    def eqn_10_12__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 5
        result = []
        T_c = T_s + 5
        result.append(T_c)
        return result
    def eqn_10_12__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 5
        result = []
        T_s = T_c - 5
        result.append(T_s)
        return result
    @kwasak
    def eqn_10_13(self, T_c=None, T_s=None):
        return
    def eqn_10_13__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 25
        result = []
        T_c = T_s + 25
        result.append(T_c)
        return result
    def eqn_10_13__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 25
        result = []
        T_s = T_c - 25
        result.append(T_s)
        return result
    @kwasak
    def eqn_10_14(self, T_c=None, T_s=None):
        return
    def eqn_10_14__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 12
        result = []
        T_c = T_s + 12
        result.append(T_c)
        return result
    def eqn_10_14__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 12
        result = []
        T_s = T_c - 12
        result.append(T_s)
        return result
    @kwasak
    def eqn_10_15(self, P=None, S_Th=None, S_p=None, p_s=None):
        return
    def eqn_10_15__P(self, S_Th: float, S_p: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        P = S_Th*p_s/(S_Th - S_p)
        result.append(P)
        return result
    def eqn_10_15__S_Th(self, P: float, S_p: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        S_Th = P*S_p/(P - p_s)
        result.append(S_Th)
        return result
    def eqn_10_15__S_p(self, P: float, S_Th: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        S_p = S_Th*(P - p_s)/P
        result.append(S_p)
        return result
    def eqn_10_15__p_s(self, P: float, S_Th: float, S_p: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        p_s = P*(S_Th - S_p)/S_Th
        result.append(p_s)
        return result
    @kwasak
    def eqn_10_16(self, P=None, S_0=None, S_Th=None, p_0=None):
        return
    def eqn_10_16__P(self, S_0: float, S_Th: float, p_0: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        P = p_0*(S_Th/S_0)**(5/3)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5/((-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5/((-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result
    def eqn_10_16__S_0(self, P: float, S_Th: float, p_0: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/(P/(P - p_0))**(3/5)
        result.append(S_0)
        return result
    def eqn_10_16__S_Th(self, P: float, S_0: float, p_0: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*(P/(P - p_0))**(3/5)
        result.append(S_Th)
        return result
    def eqn_10_16__p_0(self, P: float, S_0: float, S_Th: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        p_0 = P - P/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = P - P/(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = P - P/(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result
    @kwasak
    def eqn_10_17(self, P=None, S_0=None, S_Th=None, p_0=None, p_s=None):
        return
    def eqn_10_17__P(self, S_0: float, S_Th: float, p_0: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        P = (p_0*(S_Th/S_0)**(5/3) - p_s)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = (p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/((-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = (p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/((-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result
    def eqn_10_17__S_0(self, P: float, S_Th: float, p_0: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/((P - p_s)/(P - p_0))**(3/5)
        result.append(S_0)
        return result
    def eqn_10_17__S_Th(self, P: float, S_0: float, p_0: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*((P - p_s)/(P - p_0))**(3/5)
        result.append(S_Th)
        return result
    def eqn_10_17__p_0(self, P: float, S_0: float, S_Th: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_0 = (P*(S_Th/S_0)**(5/3) - P + p_s)/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = (P*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = (P*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result
    def eqn_10_17__p_s(self, P: float, S_0: float, S_Th: float, p_0: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_s = -P*(S_Th/S_0)**(5/3) + P + p_0*(S_Th/S_0)**(5/3)
        result.append(p_s)
        p_s = -P*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 + P + p_0*(-0.5*(S_Th/S_0)**0.333333333333333 - 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        p_s = -P*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5 + P + p_0*(-0.5*(S_Th/S_0)**0.333333333333333 + 0.866025403784439*I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        return result
    @kwasak
    def eqn_10_18(self, P=None, S_Th=None, S_p=None, T_e=None, T_i=None, p_c=None, p_s=None):
        """
        T_i := inlet  temperature of load
        """
        return
    def eqn_10_18__P(self, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        P = (S_Th*T_i*p_s + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*T_i + 460*S_Th - S_p*T_e - 460*S_p)
        result.append(P)
        return result
    def eqn_10_18__S_Th(self, P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_Th = S_p*(P*T_e + 460*P - T_e*p_c - 460*p_c)/(P*T_i + 460*P - T_i*p_s - 460*p_s)
        result.append(S_Th)
        return result
    def eqn_10_18__S_p(self, P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_p = S_Th*(P*T_i + 460*P - T_i*p_s - 460*p_s)/(P*T_e + 460*P - T_e*p_c - 460*p_c)
        result.append(S_p)
        return result
    def eqn_10_18__T_e(self, P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_e = (P*S_Th*T_i + 460*P*S_Th - 460*P*S_p - S_Th*T_i*p_s - 460*S_Th*p_s + 460*S_p*p_c)/(S_p*(P - p_c))
        result.append(T_e)
        return result
    def eqn_10_18__T_i(self, P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_i = (-460*P*S_Th + P*S_p*T_e + 460*P*S_p + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*(P - p_s))
        result.append(T_i)
        return result
    def eqn_10_18__p_c(self, P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_c = (-P*S_Th*T_i - 460*P*S_Th + P*S_p*T_e + 460*P*S_p + S_Th*T_i*p_s + 460*S_Th*p_s)/(S_p*(T_e + 460))
        result.append(p_c)
        return result
    def eqn_10_18__p_s(self, P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, **kwargs):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_s = (P*S_Th*T_i + 460*P*S_Th - P*S_p*T_e - 460*P*S_p + S_p*T_e*p_c + 460*S_p*p_c)/(S_Th*(T_i + 460))
        result.append(p_s)
        return result
    @kwasak
    def eqn_10_19(self, P=None, S_Th=None, S_p=None, T_e=None, T_i=None, p_c=None, p_s=None):
        return
    def eqn_10_19__P(self, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        # Solve for P:
        # Step 1: (S_p / S_Th) ** (1.666667) = (P - p_s)*(460 + T_i) / ( (P - p_c)*(460 + T_e) )
        R = (S_p / (S_Th)) ** (1.666667)
        # Step 2: R * ((460 + T_e)) * (P - p_c) = ((460 + T_i)) * (P - p_s)
        # Step 3: P * (R * ((460 + T_e)) - ((460 + T_i))) = R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s
        # Step 4: P = (R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s) / (R * ((460 + T_e)) - ((460 + T_i)))
        P = (R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s) / (R * ((460 + T_e)) - ((460 + T_i)))
        return [P]
    def eqn_10_19__S_Th(self, P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_Th = S_p/((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_Th)
        return result
    def eqn_10_19__S_p(self, P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_p = S_Th*((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_p)
        return result
    def eqn_10_19__T_e(self, P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_e = (P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P - T_i*p_s + 460.0*p_c*(S_p/S_Th)**(5/3) - 460.0*p_s)/((S_p/S_Th)**(5/3)*(P - p_c))
        result.append(T_e)
        T_e = (P*T_i - 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        T_e = (P*T_i - 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        return result
    def eqn_10_19__T_i(self, P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float, **kwargs):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_i = (P*T_e*(S_p/S_Th)**(5/3) + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P - T_e*p_c*(S_p/S_Th)**(5/3) - 460.0*p_c*(S_p/S_Th)**(5/3) + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - T_e*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - T_e*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        return result
    def eqn_10_19__p_c(self, P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float, **kwargs):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_c = (P*T_e*(S_p/S_Th)**(5/3) - P*T_i + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P + T_i*p_s + 460.0*p_s)/((S_p/S_Th)**(5/3)*(T_e + 460.0))
        result.append(p_c)
        p_c = (P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        p_c = (P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        return result
    def eqn_10_19__p_s(self, P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, **kwargs):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_s = (-P*T_e*(S_p/S_Th)**(5/3) + P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P + T_e*p_c*(S_p/S_Th)**(5/3) + 460.0*p_c*(S_p/S_Th)**(5/3))/(T_i + 460.0)
        result.append(p_s)
        p_s = (-P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + T_e*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        p_s = (-P*T_e*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 460.0*P*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + T_e*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        return result
    @kwasak
    def eqn_10_2(self, PS=None, Q_gas=None, V=None, dP=None, dt=None):
        return
    def eqn_10_2__PS(self, Q_gas: float, V: float, dP: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        PS = Q_gas - V*dP/dt
        result.append(PS)
        return result
    def eqn_10_2__Q_gas(self, PS: float, V: float, dP: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        Q_gas = PS + V*dP/dt
        result.append(Q_gas)
        return result
    def eqn_10_2__V(self, PS: float, Q_gas: float, dP: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        V = dt*(-PS + Q_gas)/dP
        result.append(V)
        return result
    def eqn_10_2__dP(self, PS: float, Q_gas: float, V: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        dP = dt*(-PS + Q_gas)/V
        result.append(dP)
        return result
    def eqn_10_2__dt(self, PS: float, Q_gas: float, V: float, dP: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        dt = -V*dP/(PS - Q_gas)
        result.append(dt)
        return result
    @kwasak
    def eqn_10_20(self, P=None, S_0=None, S_p=None, T_e=None, T_i=None, p_0=None, p_c=None, p_s=None):
        return
    def eqn_10_20__P(self, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # P appears 4 times — use numerical solver
        from scipy.optimize import brentq
        def _res(P_val):
            return S_p * ((P_val - p_0)*(460 + T_i) * (P_val - p_c) / (P_val * (P_val - p_s)*(460 + T_e) ) )**0.6 - S_0
        lo, hi = None, None
        prev = _res(0.01)
        for i in range(1, 100000):
            x = i * 0.01
            try:
                cur = _res(x)
            except Exception:
                continue
            if prev * cur < 0:
                lo, hi = x - 0.01, x
                break
            prev = cur
        if lo is None:
            raise UnsolvedException("No sign change found for P")
        P = brentq(_res, lo, hi)
        return [P]
    def eqn_10_20__S_0(self, P: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_0 = S_p*((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_0)
        return result
    def eqn_10_20__S_p(self, P: float, S_0: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_p = S_0/((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_p)
        return result
    def eqn_10_20__T_e(self, P: float, S_0: float, S_p: float, T_i: float, p_0: float, p_c: float, p_s: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for T_e:
        R = (S_0 / (S_p)) ** (1.666667)
        # (460 + T_e) = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (P - p_s)))
        # T_e = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (P - p_s))) - 460
        T_e = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (P - p_s))) - 460
        return [T_e]
    def eqn_10_20__T_i(self, P: float, S_0: float, S_p: float, T_e: float, p_0: float, p_c: float, p_s: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for T_i:
        R = (S_0 / (S_p)) ** (1.666667)
        # (460 + T_i) = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(P - p_c))
        # T_i = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(P - p_c)) - 460
        T_i = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(P - p_c)) - 460
        return [T_i]
    def eqn_10_20__p_0(self, P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for p_0:
        R = (S_0 / (S_p)) ** (1.666667)
        # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) 
        # (P - p_0) = R * (P * (P - p_s)*(460 + T_e) ) / ((460 + T_i) * (P - p_c))
        # p_0 = P - R * (P * (P - p_s)*(460 + T_e) ) / ((460 + T_i) * (P - p_c))
        p_0 = P - R * (P * (P - p_s)*(460 + T_e) ) / ((460 + T_i) * (P - p_c))
        return [p_0]
    def eqn_10_20__p_c(self, P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_s: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for p_c:
        R = (S_0 / (S_p)) ** (1.666667)
        # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) 
        # (P - p_c) = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
        # p_c = P - R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
        p_c = P - R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
        return [p_c]
    def eqn_10_20__p_s(self, P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, **kwargs):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for p_s:
        R = (S_0 / (S_p)) ** (1.666667)
        # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) 
        # (P - p_s) = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
        # p_s = P - ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
        p_s = P - ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
        return [p_s]
    @kwasak
    def eqn_10_21(self, P=None, P_d=None, P_prime=None):
        """
        P_prime := pseudo suction pressure
        P_d := actual pump discharge pressure
        """
        return
    def eqn_10_21__P(self, P_d: float, P_prime: float, **kwargs):
        # P_prime = P / P_d * 760
        result = []
        P = P_d*P_prime/760
        result.append(P)
        return result
    def eqn_10_21__P_d(self, P: float, P_prime: float, **kwargs):
        # P_prime = P / P_d * 760
        result = []
        P_d = 760*P/P_prime
        result.append(P_d)
        return result
    def eqn_10_21__P_prime(self, P: float, P_d: float, **kwargs):
        # P_prime = P / P_d * 760
        result = []
        P_prime = 760*P/P_d
        result.append(P_prime)
        return result
    @kwasak
    def eqn_10_3(self, N_mfw=None, Q_gas=None, T=None):
        return
    def eqn_10_3__N_mfw(self, Q_gas: float, T: float, **kwargs):
        # Q_gas = 9.25 * N_mfw * T
        result = []
        N_mfw = 0.108108108108108*Q_gas/T
        result.append(N_mfw)
        return result
    def eqn_10_3__Q_gas(self, N_mfw: float, T: float, **kwargs):
        # Q_gas = 9.25 * N_mfw * T
        result = []
        Q_gas = 9.25*N_mfw*T
        result.append(Q_gas)
        return result
    def eqn_10_3__T(self, N_mfw: float, Q_gas: float, **kwargs):
        # Q_gas = 9.25 * N_mfw * T
        result = []
        T = 0.108108108108108*Q_gas/N_mfw
        result.append(T)
        return result
    @kwasak
    def eqn_10_4(self, Q_gas=None, SP_1=None, SP_2=None, S_p=None, V=None, t=None):
        return
    def eqn_10_4__Q_gas(self, SP_1: float, SP_2: float, S_p: float, V: float, t: float, **kwargs):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        Q_gas = -(SP_1 - SP_2*exp(S_p*t/V))/(exp(S_p*t/V) - 1)
        result.append(Q_gas)
        return result
    def eqn_10_4__SP_1(self, Q_gas: float, SP_2: float, S_p: float, V: float, t: float, **kwargs):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_1 = Q_gas + (-Q_gas + SP_2)*exp(S_p*t/V)
        result.append(SP_1)
        return result
    def eqn_10_4__SP_2(self, Q_gas: float, SP_1: float, S_p: float, V: float, t: float, **kwargs):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_2 = (Q_gas*exp(S_p*t/V) - Q_gas + SP_1)*exp(-S_p*t/V)
        result.append(SP_2)
        return result
    def eqn_10_4__S_p(self, Q_gas: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        S_p = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/t
        result.append(S_p)
        return result
    def eqn_10_4__V(self, Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float, **kwargs):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        V = S_p*t/log((Q_gas - SP_1)/(Q_gas - SP_2))
        result.append(V)
        return result
    def eqn_10_4__t(self, Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float, **kwargs):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        t = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/S_p
        result.append(t)
        return result
    @kwasak
    def eqn_10_5(self, P_1=None, P_2=None, S_p=None, V=None, t=None):
        return
    def eqn_10_5__P_1(self, P_2: float, S_p: float, V: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_p*t/V)
        result.append(P_1)
        return result
    def eqn_10_5__P_2(self, P_1: float, S_p: float, V: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_p*t/V)
        result.append(P_2)
        return result
    def eqn_10_5__S_p(self, P_1: float, P_2: float, V: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        S_p = V*log(P_1/P_2)/t
        result.append(S_p)
        return result
    def eqn_10_5__V(self, P_1: float, P_2: float, S_p: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        V = S_p*t/log(P_1/P_2)
        result.append(V)
        return result
    def eqn_10_5__t(self, P_1: float, P_2: float, S_p: float, V: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_p
        result.append(t)
        return result
    @kwasak
    def eqn_10_6(self, P_1=None, P_2=None, S_a=None, V=None, t=None):
        return
    def eqn_10_6__P_1(self, P_2: float, S_a: float, V: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_a*t/V)
        result.append(P_1)
        return result
    def eqn_10_6__P_2(self, P_1: float, S_a: float, V: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_a*t/V)
        result.append(P_2)
        return result
    def eqn_10_6__S_a(self, P_1: float, P_2: float, V: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        S_a = V*log(P_1/P_2)/t
        result.append(S_a)
        return result
    def eqn_10_6__V(self, P_1: float, P_2: float, S_a: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        V = S_a*t/log(P_1/P_2)
        result.append(V)
        return result
    def eqn_10_6__t(self, P_1: float, P_2: float, S_a: float, V: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_a
        result.append(t)
        return result
    @kwasak
    def eqn_10_8(self, bhp=None, c_p=None, delta_T=None, delta_h_i=None, f_a=None, rho=None, w_i=None):
        return
    def eqn_10_8__bhp(self, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        bhp = 0.00315127701375246*c_p*delta_T*f_a*rho - 0.000392927308447937*delta_h_i*w_i
        result.append(bhp)
        return result
    def eqn_10_8__c_p(self, bhp: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        c_p = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(delta_T*f_a*rho)
        result.append(c_p)
        return result
    def eqn_10_8__delta_T(self, bhp: float, c_p: float, delta_h_i: float, f_a: float, rho: float, w_i: float, **kwargs):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_T = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*f_a*rho)
        result.append(delta_T)
        return result
    def eqn_10_8__delta_h_i(self, bhp: float, c_p: float, delta_T: float, f_a: float, rho: float, w_i: float, **kwargs):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_h_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/w_i
        result.append(delta_h_i)
        return result
    def eqn_10_8__f_a(self, bhp: float, c_p: float, delta_T: float, delta_h_i: float, rho: float, w_i: float, **kwargs):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        f_a = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*rho)
        result.append(f_a)
        return result
    def eqn_10_8__rho(self, bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, w_i: float, **kwargs):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        rho = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*f_a)
        result.append(rho)
        return result
    def eqn_10_8__w_i(self, bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, **kwargs):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        w_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/delta_h_i
        result.append(w_i)
        return result
    @kwasak
    def eqn_10_9(self, T_c=None, T_s=None, delta_T=None):
        return
    def eqn_10_9__T_c(self, T_s: float, delta_T: float, **kwargs):
        # T_c = T_s + delta_T
        result = []
        T_c = T_s + delta_T
        result.append(T_c)
        return result
    def eqn_10_9__T_s(self, T_c: float, delta_T: float, **kwargs):
        # T_c = T_s + delta_T
        result = []
        T_s = T_c - delta_T
        result.append(T_s)
        return result
    def eqn_10_9__delta_T(self, T_c: float, T_s: float, **kwargs):
        # T_c = T_s + delta_T
        result = []
        delta_T = T_c - T_s
        result.append(delta_T)
        return result

class Precondensors:
    @kwasak
    def eqn_7_1(self, P=None, p_i=None, y_i=None):
        return
    def eqn_7_1__P(self, p_i: float, y_i: float, **kwargs):
        # y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return result
    def eqn_7_1__p_i(self, P: float, y_i: float, **kwargs):
        # y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
        return result
    def eqn_7_1__y_i(self, P: float, p_i: float, **kwargs):
        # y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return result
    @kwasak
    def eqn_7_10(self, L_c_P=None, Q_condensor_heat_duty=None, del_T=None):
        return
    def eqn_7_10__L_c_P(self, Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty/(500*del_T)
        result.append(L_c_P)
        return result
    def eqn_7_10__Q_condensor_heat_duty(self, L_c_P: float, del_T: float, **kwargs):
        # L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500*L_c_P*del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_10__del_T(self, L_c_P: float, Q_condensor_heat_duty: float, **kwargs):
        # L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(500*L_c_P)
        result.append(del_T)
        return result
    @kwasak
    def eqn_7_11(self, Q_condensor_heat_duty=None, U_v=None, V_c=None, del_T_LM=None):
        return
    def eqn_7_11__Q_condensor_heat_duty(self, U_v: float, V_c: float, del_T_LM: float, **kwargs):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v*V_c*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_11__U_v(self, Q_condensor_heat_duty: float, V_c: float, del_T_LM: float, **kwargs):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
        result.append(U_v)
        return result
    def eqn_7_11__V_c(self, Q_condensor_heat_duty: float, U_v: float, del_T_LM: float, **kwargs):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty/(U_v*del_T_LM)
        result.append(V_c)
        return result
    def eqn_7_11__del_T_LM(self, Q_condensor_heat_duty: float, U_v: float, V_c: float, **kwargs):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(U_v*V_c)
        result.append(del_T_LM)
        return result
    @kwasak
    def eqn_7_12(self, A=None, Q_condensor_heat_duty=None, U=None, del_T=None):
        return
    def eqn_7_12__A(self, Q_condensor_heat_duty: float, U: float, del_T: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty/(U*del_T)
        result.append(A)
        return result
    def eqn_7_12__Q_condensor_heat_duty(self, A: float, U: float, del_T: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A*U*del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_12__U(self, A: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty/(A*del_T)
        result.append(U)
        return result
    def eqn_7_12__del_T(self, A: float, Q_condensor_heat_duty: float, U: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty/(A*U)
        result.append(del_T)
        return result
    @kwasak
    def eqn_7_14a(self, A=None, Q_condensor_heat_duty=None, U=None, del_T_LM=None):
        return
    def eqn_7_14a__A(self, Q_condensor_heat_duty: float, U: float, del_T_LM: float, **kwargs):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty/(U*del_T_LM)
        result.append(A)
        return result
    def eqn_7_14a__Q_condensor_heat_duty(self, A: float, U: float, del_T_LM: float, **kwargs):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A*U*del_T_LM
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_14a__U(self, A: float, Q_condensor_heat_duty: float, del_T_LM: float, **kwargs):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty/(A*del_T_LM)
        result.append(U)
        return result
    def eqn_7_14a__del_T_LM(self, A: float, Q_condensor_heat_duty: float, U: float, **kwargs):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(A*U)
        result.append(del_T_LM)
        return result
    @kwasak
    def eqn_7_14b(self, A=None, Q_condensor_heat_duty=None, U=None, del_T_1=None, del_T_2=None):
        return
    def eqn_7_14b__A(self, Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float, **kwargs):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(A)
        return result
    def eqn_7_14b__Q_condensor_heat_duty(self, A: float, U: float, del_T_1: float, del_T_2: float, **kwargs):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A*U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_14b__U(self, A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float, **kwargs):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty/(A*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(U)
        return result
    def eqn_7_14b__del_T_1(self, A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float, **kwargs):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_1)
        return result
    def eqn_7_14b__del_T_2(self, A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float, **kwargs):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_2)
        return result
    @kwasak
    def eqn_7_15(self, U=None, sum_R=None):
        return
    def eqn_7_15__U(self, sum_R: float, **kwargs):
        # 1 / U = sum_R
        result = []
        U = 1/sum_R
        result.append(U)
        return result
    def eqn_7_15__sum_R(self, U: float, **kwargs):
        # 1 / U = sum_R
        result = []
        sum_R = 1/U
        result.append(sum_R)
        return result
    @kwasak
    def eqn_7_16(self, D_0=None, D_LM=None, D_i=None, R_f_0=None, R_fi=None, U_0=None, h_0=None, h_i=None, k_w=None, x_w=None):
        return
    def eqn_7_16__D_0(self, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_f_0*U_0*h_0 - U_0 + h_0)/(U_0*h_0*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result
    def eqn_7_16__D_LM(self, D_0: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_0*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(D_LM)
        return result
    def eqn_7_16__D_i(self, D_0: float, D_LM: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_0*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_0*x_w + D_LM*R_f_0*U_0*h_0*k_w + D_LM*U_0*k_w - D_LM*h_0*k_w))
        result.append(D_i)
        return result
    def eqn_7_16__R_f_0(self, D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
        result.append(R_f_0)
        return result
    def eqn_7_16__R_fi(self, D_0: float, D_LM: float, D_i: float, R_f_0: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_f_0/D_0 - D_i/(D_0*h_0) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result
    def eqn_7_16__U_0(self, D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, h_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_0*h_i*k_w/(D_0*D_LM*R_fi*h_0*h_i*k_w + D_0*D_LM*h_0*k_w + D_0*D_i*h_0*h_i*x_w + D_LM*D_i*R_f_0*h_0*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result
    def eqn_7_16__h_0(self, D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_0 = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_f_0*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_0)
        return result
    def eqn_7_16__h_i(self, D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_0*k_w/(D_0*D_LM*R_fi*U_0*h_0*k_w + D_0*D_i*U_0*h_0*x_w + D_LM*D_i*R_f_0*U_0*h_0*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_0*k_w)
        result.append(h_i)
        return result
    def eqn_7_16__k_w(self, D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, x_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_0*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(k_w)
        return result
    def eqn_7_16__x_w(self, D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, **kwargs):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_f_0*k_w/D_0 - D_LM*k_w/(D_0*h_0) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result
    @kwasak
    def eqn_7_17(self, R_0=None, R_nc=None, h_c=None):
        return
    def eqn_7_17__R_0(self, R_nc: float, h_c: float, **kwargs):
        # R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1/h_c
        result.append(R_0)
        return result
    def eqn_7_17__R_nc(self, R_0: float, h_c: float, **kwargs):
        # R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1/h_c
        result.append(R_nc)
        return result
    def eqn_7_17__h_c(self, R_0: float, R_nc: float, **kwargs):
        # R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1/(R_0 - R_nc)
        result.append(h_c)
        return result
    @kwasak
    def eqn_7_18(self, D_0=None, D_LM=None, D_i=None, R_fi=None, R_fo=None, R_nc=None, U_0=None, h_c=None, h_i=None, k_w=None, x_w=None):
        return
    def eqn_7_18__D_0(self, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_fo*U_0*h_c - R_nc*U_0*h_c - U_0 + h_c)/(U_0*h_c*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return result
    def eqn_7_18__D_LM(self, D_0: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_c*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(D_LM)
        return result
    def eqn_7_18__D_i(self, D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_c*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_c*x_w + D_LM*R_fo*U_0*h_c*k_w + D_LM*R_nc*U_0*h_c*k_w + D_LM*U_0*k_w - D_LM*h_c*k_w))
        result.append(D_i)
        return result
    def eqn_7_18__R_fi(self, D_0: float, D_LM: float, D_i: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_fo/D_0 - D_i*R_nc/D_0 - D_i/(D_0*h_c) + D_i/(D_0*U_0)
        result.append(R_fi)
        return result
    def eqn_7_18__R_fo(self, D_0: float, D_LM: float, D_i: float, R_fi: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fo = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_nc - 1/h_c + 1/U_0
        result.append(R_fo)
        return result
    def eqn_7_18__R_nc(self, D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_nc = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_fo - 1/h_c + 1/U_0
        result.append(R_nc)
        return result
    def eqn_7_18__U_0(self, D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, h_c: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_c*h_i*k_w/(D_0*D_LM*R_fi*h_c*h_i*k_w + D_0*D_LM*h_c*k_w + D_0*D_i*h_c*h_i*x_w + D_LM*D_i*R_fo*h_c*h_i*k_w + D_LM*D_i*R_nc*h_c*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return result
    def eqn_7_18__h_c(self, D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_i: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_c = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_fo*U_0*h_i*k_w + D_LM*D_i*R_nc*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_c)
        return result
    def eqn_7_18__h_i(self, D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, k_w: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_c*k_w/(D_0*D_LM*R_fi*U_0*h_c*k_w + D_0*D_i*U_0*h_c*x_w + D_LM*D_i*R_fo*U_0*h_c*k_w + D_LM*D_i*R_nc*U_0*h_c*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_c*k_w)
        result.append(h_i)
        return result
    def eqn_7_18__k_w(self, D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, x_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_c*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(k_w)
        return result
    def eqn_7_18__x_w(self, D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, **kwargs):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_fo*k_w/D_0 - D_LM*R_nc*k_w/D_0 - D_LM*k_w/(D_0*h_c) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return result
    @kwasak
    def eqn_7_2(self, P_i_0=None, p_i=None, x_i=None):
        return
    def eqn_7_2__P_i_0(self, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i/x_i
        result.append(P_i_0)
        return result
    def eqn_7_2__p_i(self, P_i_0: float, x_i: float, **kwargs):
        # p_i = x_i * P_i_0
        result = []
        p_i = P_i_0*x_i
        result.append(p_i)
        return result
    def eqn_7_2__x_i(self, P_i_0: float, p_i: float, **kwargs):
        # p_i = x_i * P_i_0
        result = []
        x_i = p_i/P_i_0
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_3(self, P_i_0=None, epsilon_i=None, p_i=None, x_i=None):
        return
    def eqn_7_3__P_i_0(self, epsilon_i: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i/(epsilon_i*x_i)
        result.append(P_i_0)
        return result
    def eqn_7_3__epsilon_i(self, P_i_0: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i/(P_i_0*x_i)
        result.append(epsilon_i)
        return result
    def eqn_7_3__p_i(self, P_i_0: float, epsilon_i: float, x_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0*epsilon_i*x_i
        result.append(p_i)
        return result
    def eqn_7_3__x_i(self, P_i_0: float, epsilon_i: float, p_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i/(P_i_0*epsilon_i)
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_4a(self, P=None, p_c=None, p_nc=None):
        return
    def eqn_7_4a__P(self, p_c: float, p_nc: float, **kwargs):
        # p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return result
    def eqn_7_4a__p_c(self, P: float, p_nc: float, **kwargs):
        # p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return result
    def eqn_7_4a__p_nc(self, P: float, p_c: float, **kwargs):
        # p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return result
    @kwasak
    def eqn_7_4aa(self, n_i=None, n_nc=None, p_i=None, p_nc=None):
        return
    def eqn_7_4aa__n_i(self, n_nc: float, p_i: float, p_nc: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        n_i = n_nc*p_i/p_nc
        result.append(n_i)
        return result
    def eqn_7_4aa__n_nc(self, n_i: float, p_i: float, p_nc: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        n_nc = n_i*p_nc/p_i
        result.append(n_nc)
        return result
    def eqn_7_4aa__p_i(self, n_i: float, n_nc: float, p_nc: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        p_i = n_i*p_nc/n_nc
        result.append(p_i)
        return result
    def eqn_7_4aa__p_nc(self, n_i: float, n_nc: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        p_nc = n_nc*p_i/n_i
        result.append(p_nc)
        return result
    @kwasak
    def eqn_7_4ab(self, P_c=None, p=None, p_i=None, p_nc=None):
        return
    def eqn_7_4ab__P_c(self, p: float, p_i: float, p_nc: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
        return result
    def eqn_7_4ab__p(self, P_c: float, p_i: float, p_nc: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result
    def eqn_7_4ab__p_i(self, P_c: float, p: float, p_nc: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = 0
        result.append(p_i)
        return result
    def eqn_7_4ab__p_nc(self, P_c: float, p: float, p_i: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = -P_c + p
        result.append(p_nc)
        return result
    @kwasak
    def eqn_7_4ac(self, P_c=None, n_i=None, n_nc=None, p=None, p_i=None):
        return
    def eqn_7_4ac__P_c(self, n_i: float, n_nc: float, p: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc*p_i/n_i
        result.append(P_c)
        return result
    def eqn_7_4ac__n_i(self, P_c: float, n_nc: float, p: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/(-P_c + p)
        result.append(n_i)
        return result
    def eqn_7_4ac__n_nc(self, P_c: float, n_i: float, p: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*(-P_c + p)/p_i
        result.append(n_nc)
        return result
    def eqn_7_4ac__p(self, P_c: float, n_i: float, n_nc: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc*p_i/n_i
        result.append(p)
        return result
    def eqn_7_4ac__p_i(self, P_c: float, n_i: float, n_nc: float, p: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*(-P_c + p)/n_nc
        result.append(p_i)
        return result
    @kwasak
    def eqn_7_5(self, N_i=None, N_nc=None, P=None, P_c=None, p_i=None):
        return
    def eqn_7_5__N_i(self, N_nc: float, P: float, P_c: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
        return result
    def eqn_7_5__N_nc(self, N_i: float, P: float, P_c: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return result
    def eqn_7_5__P(self, N_i: float, N_nc: float, P_c: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc*p_i/N_i
        result.append(P)
        return result
    def eqn_7_5__P_c(self, N_i: float, N_nc: float, P: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc*p_i/N_i
        result.append(P_c)
        return result
    def eqn_7_5__p_i(self, N_i: float, N_nc: float, P: float, P_c: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i*(P - P_c)/N_nc
        result.append(p_i)
        return result
    @kwasak
    def eqn_7_6(self, M=None, P=None, P_i_0=None, W_air=None, W_i=None, p_c=None, x_i=None):
        return
    def eqn_7_6__M(self, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*x_i)
        result.append(M)
        return result
    def eqn_7_6__P(self, M: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*x_i/(29*W_i) + p_c
        result.append(P)
        return result
    def eqn_7_6__P_i_0(self, M: float, P: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*x_i)
        result.append(P_i_0)
        return result
    def eqn_7_6__W_air(self, M: float, P: float, P_i_0: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*x_i)
        result.append(W_air)
        return result
    def eqn_7_6__W_i(self, M: float, P: float, P_i_0: float, W_air: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*x_i/(29*(P - p_c))
        result.append(W_i)
        return result
    def eqn_7_6__p_c(self, M: float, P: float, P_i_0: float, W_air: float, W_i: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*x_i/(29*W_i) + P
        result.append(p_c)
        return result
    def eqn_7_6__x_i(self, M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, **kwargs):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air)
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_7(self, M=None, P=None, P_i_0=None, W_air=None, W_i=None, epsilon_i=None, p_c=None, x_i=None):
        return
    def eqn_7_7__M(self, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
        result.append(M)
        return result
    def eqn_7_7__P(self, M: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + p_c
        result.append(P)
        return result
    def eqn_7_7__P_i_0(self, M: float, P: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*epsilon_i*x_i)
        result.append(P_i_0)
        return result
    def eqn_7_7__W_air(self, M: float, P: float, P_i_0: float, W_i: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*epsilon_i*x_i)
        result.append(W_air)
        return result
    def eqn_7_7__W_i(self, M: float, P: float, P_i_0: float, W_air: float, epsilon_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*epsilon_i*x_i/(29*(P - p_c))
        result.append(W_i)
        return result
    def eqn_7_7__epsilon_i(self, M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*x_i)
        result.append(epsilon_i)
        return result
    def eqn_7_7__p_c(self, M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, x_i: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + P
        result.append(p_c)
        return result
    def eqn_7_7__x_i(self, M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, **kwargs):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*epsilon_i)
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_8(self, L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None):
        return
    def eqn_7_8__L_c(self, Q_condensor_heat_duty: float, c_p: float, del_T: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return result
    def eqn_7_8__Q_condensor_heat_duty(self, L_c: float, c_p: float, del_T: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_8__c_p(self, L_c: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return result
    def eqn_7_8__del_T(self, L_c: float, Q_condensor_heat_duty: float, c_p: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return result
    @kwasak
    def eqn_7_9(self, L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, rho=None):
        return
    def eqn_7_9__L_c(self, Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return result
    def eqn_7_9__Q_condensor_heat_duty(self, L_c: float, c_p: float, del_T: float, rho: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_9__c_p(self, L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return result
    def eqn_7_9__del_T(self, L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return result
    def eqn_7_9__rho(self, L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float, **kwargs):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return result

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
        K_i = y_i/x_i
        result.append(K_i)
        return result
    def eqn_5_1__x_i(self, K_i: float, y_i: float, **kwargs):
        # K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
        return result
    def eqn_5_1__y_i(self, K_i: float, x_i: float, **kwargs):
        # K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
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
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return result
    def eqn_5_10b__R(self, L_0: float, V_1: float, **kwargs):
        # L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return result
    def eqn_5_10b__V_1(self, L_0: float, R: float, **kwargs):
        # L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return result
    @kwasak
    def eqn_5_10c(self, D=None, L_0=None, R=None):
        return
    def eqn_5_10c__D(self, L_0: float, R: float, **kwargs):
        # (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result
    def eqn_5_10c__L_0(self, D: float, R: float, **kwargs):
        # (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
        return result
    def eqn_5_10c__R(self, D: float, L_0: float, **kwargs):
        # (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
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
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result
    def eqn_5_12__N_ES(self, Eff: float, N_t: float, T: float, **kwargs):
        # N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result
    def eqn_5_12__N_t(self, Eff: float, N_ES: float, T: float, **kwargs):
        # N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result
    def eqn_5_12__T(self, Eff: float, N_ES: float, N_t: float, **kwargs):
        # N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result
    @kwasak
    def eqn_5_13(self, HETP=None, H_p=None, N_ES=None):
        return
    def eqn_5_13__HETP(self, H_p: float, N_ES: float, **kwargs):
        # H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result
    def eqn_5_13__H_p(self, HETP: float, N_ES: float, **kwargs):
        # H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result
    def eqn_5_13__N_ES(self, HETP: float, H_p: float, **kwargs):
        # H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result
    @kwasak
    def eqn_5_14(self, M=None, P_0=None, T=None, W_E=None):
        return
    def eqn_5_14__M(self, P_0: float, T: float, W_E: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return result
    def eqn_5_14__P_0(self, M: float, T: float, W_E: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return result
    def eqn_5_14__T(self, M: float, P_0: float, W_E: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return result
    def eqn_5_14__W_E(self, M: float, P_0: float, T: float, **kwargs):
        # W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return result
    @kwasak
    def eqn_5_15(self, M_1=None, M_2=None, P_0_1=None, P_0_2=None, a_M_12=None):
        return
    def eqn_5_15__M_1(self, M_2: float, P_0_1: float, P_0_2: float, a_M_12: float, **kwargs):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        M_1 = M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        return result
    def eqn_5_15__M_2(self, M_1: float, P_0_1: float, P_0_2: float, a_M_12: float, **kwargs):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_2 = -M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        M_2 = M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        return result
    def eqn_5_15__P_0_1(self, M_1: float, M_2: float, P_0_2: float, a_M_12: float, **kwargs):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
        result.append(P_0_1)
        return result
    def eqn_5_15__P_0_2(self, M_1: float, M_2: float, P_0_1: float, a_M_12: float, **kwargs):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1*(M_2/M_1)**(2/5)/a_M_12
        result.append(P_0_2)
        return result
    def eqn_5_15__a_M_12(self, M_1: float, M_2: float, P_0_1: float, P_0_2: float, **kwargs):
        # a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
        result.append(a_M_12)
        return result
    @kwasak
    def eqn_5_16(self, H_i=None, p_i=None, x_i=None):
        return
    def eqn_5_16__H_i(self, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result
    def eqn_5_16__p_i(self, H_i: float, x_i: float, **kwargs):
        # p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result
    def eqn_5_16__x_i(self, H_i: float, p_i: float, **kwargs):
        # p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result
    @kwasak
    def eqn_5_17(self, H_2_1=None, H_2_3=None, H_2_mi=None, x_1=None, x_3=None):
        return
    def eqn_5_17__H_2_1(self, H_2_3: float, H_2_mi: float, x_1: float, x_3: float, **kwargs):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3*log(H_2_3) + log(H_2_mi))/x_1)
        result.append(H_2_1)
        return result
    def eqn_5_17__H_2_3(self, H_2_1: float, H_2_mi: float, x_1: float, x_3: float, **kwargs):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1*log(H_2_1) + log(H_2_mi))/x_3)
        result.append(H_2_3)
        return result
    def eqn_5_17__H_2_mi(self, H_2_1: float, H_2_3: float, x_1: float, x_3: float, **kwargs):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1*log(H_2_1) + x_3*log(H_2_3))
        result.append(H_2_mi)
        return result
    def eqn_5_17__x_1(self, H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float, **kwargs):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return result
    def eqn_5_17__x_3(self, H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float, **kwargs):
        # log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return result
    @kwasak
    def eqn_5_2a(self, K_1=None, K_2=None, alpha_1_2=None):
        return
    def eqn_5_2a__K_1(self, K_2: float, alpha_1_2: float, **kwargs):
        # alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2*alpha_1_2
        result.append(K_1)
        return result
    def eqn_5_2a__K_2(self, K_1: float, alpha_1_2: float, **kwargs):
        # alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1/alpha_1_2
        result.append(K_2)
        return result
    def eqn_5_2a__alpha_1_2(self, K_1: float, K_2: float, **kwargs):
        # alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1/K_2
        result.append(alpha_1_2)
        return result
    @kwasak
    def eqn_5_2b(self, K_1=None, K_2=None, x_1=None, x_2=None, y_1=None, y_2=None):
        return
    def eqn_5_2b__K_1(self, K_2: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result
    def eqn_5_2b__K_2(self, K_1: float, x_1: float, x_2: float, y_1: float, y_2: float, **kwargs):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result
    def eqn_5_2b__x_1(self, K_1: float, K_2: float, x_2: float, y_1: float, y_2: float, **kwargs):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result
    def eqn_5_2b__x_2(self, K_1: float, K_2: float, x_1: float, y_1: float, y_2: float, **kwargs):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result
    def eqn_5_2b__y_1(self, K_1: float, K_2: float, x_1: float, x_2: float, y_2: float, **kwargs):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result
    def eqn_5_2b__y_2(self, K_1: float, K_2: float, x_1: float, x_2: float, y_1: float, **kwargs):
        # K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
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
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return result
    def eqn_5_3__p_i(self, P_0_i: float, x_i: float, **kwargs):
        # p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return result
    def eqn_5_3__x_i(self, P_0_i: float, p_i: float, **kwargs):
        # p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return result
    @kwasak
    def eqn_5_4(self, P=None, P_0_i=None, x_i=None, y_i=None):
        return
    def eqn_5_4__P(self, P_0_i: float, x_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return result
    def eqn_5_4__P_0_i(self, P: float, x_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        P_0_i = P*y_i/x_i
        result.append(P_0_i)
        return result
    def eqn_5_4__x_i(self, P: float, P_0_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        x_i = P*y_i/P_0_i
        result.append(x_i)
        return result
    def eqn_5_4__y_i(self, P: float, P_0_i: float, x_i: float, **kwargs):
        # y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
        return result
    @kwasak
    def eqn_5_5(self, P_0_1=None, P_0_2=None, alpha_12=None):
        return
    def eqn_5_5__P_0_1(self, P_0_2: float, alpha_12: float, **kwargs):
        # alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2*alpha_12
        result.append(P_0_1)
        return result
    def eqn_5_5__P_0_2(self, P_0_1: float, alpha_12: float, **kwargs):
        # alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1/alpha_12
        result.append(P_0_2)
        return result
    def eqn_5_5__alpha_12(self, P_0_1: float, P_0_2: float, **kwargs):
        # alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return result
    @kwasak
    def eqn_5_6(self, P_0_i=None, gamma_i=None, p_i=None, x_i=None):
        return
    def eqn_5_6__P_0_i(self, gamma_i: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result
    def eqn_5_6__gamma_i(self, P_0_i: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result
    def eqn_5_6__p_i(self, P_0_i: float, gamma_i: float, x_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return result
    def eqn_5_6__x_i(self, P_0_i: float, gamma_i: float, p_i: float, **kwargs):
        # p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result
    @kwasak
    def eqn_5_7(self, P=None, P_0_i=None, gamma_i=None, x_i=None, y_i=None):
        return
    def eqn_5_7__P(self, P_0_i: float, gamma_i: float, x_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i*gamma_i*x_i/y_i
        result.append(P)
        return result
    def eqn_5_7__P_0_i(self, P: float, gamma_i: float, x_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P*y_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result
    def eqn_5_7__gamma_i(self, P: float, P_0_i: float, x_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P*y_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result
    def eqn_5_7__x_i(self, P: float, P_0_i: float, gamma_i: float, y_i: float, **kwargs):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result
    def eqn_5_7__y_i(self, P: float, P_0_i: float, gamma_i: float, x_i: float, **kwargs):
        # y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return result
    @kwasak
    def eqn_5_8(self, P_0_1=None, P_0_2=None, alpha_12=None, gamma_1=None, gamma_2=None):
        return
    def eqn_5_8__P_0_1(self, P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float, **kwargs):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result
    def eqn_5_8__P_0_2(self, P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float, **kwargs):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result
    def eqn_5_8__alpha_12(self, P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float, **kwargs):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result
    def eqn_5_8__gamma_1(self, P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float, **kwargs):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return result
    def eqn_5_8__gamma_2(self, P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float, **kwargs):
        # alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
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

class ProcessApp2:
    @kwasak
    def eqn_6_1(self, T_1=None, T_2=None, T_R=None, c_p=None, del_h_v=None, w_1=None, w_2=None, w_v=None):
        return
    def eqn_6_1__T_1(self, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R) + del_h_v*w_v)/(c_p*w_1)
        result.append(T_1)
        return result
    def eqn_6_1__T_2(self, T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_2 = (T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R) + del_h_v*w_v)/(c_p*w_2)
        result.append(T_2)
        return result
    def eqn_6_1__T_R(self, T_1: float, T_2: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        T_R = (T_1*c_p*w_1 + T_2*c_p*w_2 - del_h_v*w_v)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result
    def eqn_6_1__c_p(self, T_1: float, T_2: float, T_R: float, del_h_v: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = del_h_v*w_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result
    def eqn_6_1__del_h_v(self, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, w_v: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/w_v
        result.append(del_h_v)
        return result
    def eqn_6_1__w_1(self, T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_2: float, w_v: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (-T_2*c_p*w_2 + T_R*c_p*w_2 + del_h_v*w_v)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result
    def eqn_6_1__w_2(self, T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_v: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_2 = (-T_1*c_p*w_1 + T_R*c_p*w_1 + del_h_v*w_v)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result
    def eqn_6_1__w_v(self, T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = w_v * del_h_v
        result = []
        w_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/del_h_v
        result.append(w_v)
        return result
    @kwasak
    def eqn_6_10(self, A=None, dV_dt=None, delta_P=None, mu=None, r_c=None, s=None, tau=None):
        return
    def eqn_6_10__A(self, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float, **kwargs):
        # dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt*delta_P**(s - 1)*mu*r_c*tau
        result.append(A)
        return result
    def eqn_6_10__dV_dt(self, A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float, **kwargs):
        # dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A*delta_P**(1 - s)/(mu*r_c*tau)
        result.append(dV_dt)
        return result
    def eqn_6_10__delta_P(self, A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float, **kwargs):
        # dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
        result.append(delta_P)
        return result
    def eqn_6_10__mu(self, A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float, **kwargs):
        # dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A*delta_P**(1 - s)/(dV_dt*r_c*tau)
        result.append(mu)
        return result
    def eqn_6_10__r_c(self, A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float, **kwargs):
        # dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
        result.append(r_c)
        return result
    def eqn_6_10__s(self, A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float, **kwargs):
        # dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A*delta_P/(dV_dt*mu*r_c*tau))/log(delta_P)
        result.append(s)
        return result
    def eqn_6_10__tau(self, A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, **kwargs):
        # dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A*delta_P**(1 - s)/(dV_dt*mu*r_c)
        result.append(tau)
        return result
    @kwasak
    def eqn_6_11a(self, A_d=None, delta_T=None, delta_h_i=None, delta_m=None, h_d=None, m_b=None, t_R=None):
        return
    def eqn_6_11a__A_d(self, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
        result.append(A_d)
        return result
    def eqn_6_11a__delta_T(self, A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
        result.append(delta_T)
        return result
    def eqn_6_11a__delta_h_i(self, A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d*delta_T*h_d*t_R/(delta_m*m_b)
        result.append(delta_h_i)
        return result
    def eqn_6_11a__delta_m(self, A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float, **kwargs):
        # t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d*delta_T*h_d*t_R/(delta_h_i*m_b)
        result.append(delta_m)
        return result
    def eqn_6_11a__h_d(self, A_d: float, delta_T: float, delta_h_i: float, delta_m: float, m_b: float, t_R: float, **kwargs):
        # t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i*delta_m*m_b/(A_d*delta_T*t_R)
        result.append(h_d)
        return result
    def eqn_6_11a__m_b(self, A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float, **kwargs):
        # t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
        result.append(m_b)
        return result
    def eqn_6_11a__t_R(self, A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, **kwargs):
        # t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        t_R = delta_h_i*delta_m*m_b/(A_d*delta_T*h_d)
        result.append(t_R)
        return result
    @kwasak
    def eqn_6_2(self, Q_v=None, T_1=None, T_2=None, T_R=None, c_p=None, w_1=None, w_2=None):
        return
    def eqn_6_2__Q_v(self, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = c_p*(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)/12000
        result.append(Q_v)
        return result
    def eqn_6_2__T_1(self, Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000*Q_v + T_R*c_p*w_1 + c_p*w_2*(-T_2 + T_R))/(c_p*w_1)
        result.append(T_1)
        return result
    def eqn_6_2__T_2(self, Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_2 = (12000*Q_v + T_R*c_p*w_2 + c_p*w_1*(-T_1 + T_R))/(c_p*w_2)
        result.append(T_2)
        return result
    def eqn_6_2__T_R(self, Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        T_R = (-12000*Q_v + T_1*c_p*w_1 + T_2*c_p*w_2)/(c_p*(w_1 + w_2))
        result.append(T_R)
        return result
    def eqn_6_2__c_p(self, Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        c_p = 12000*Q_v/(T_1*w_1 + T_2*w_2 - T_R*w_1 - T_R*w_2)
        result.append(c_p)
        return result
    def eqn_6_2__w_1(self, Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_2: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_1 = (12000*Q_v - T_2*c_p*w_2 + T_R*c_p*w_2)/(c_p*(T_1 - T_R))
        result.append(w_1)
        return result
    def eqn_6_2__w_2(self, Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, **kwargs):
        # w_1 * c_p * (T_1 - T_R) + w_2 * c_p * (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000*Q_v - T_1*c_p*w_1 + T_R*c_p*w_1)/(c_p*(T_2 - T_R))
        result.append(w_2)
        return result
    @kwasak
    def eqn_6_4(self, Q_v=None, delta_h_v=None, w_v=None):
        return
    def eqn_6_4__Q_v(self, delta_h_v: float, w_v: float, **kwargs):
        # w_v = 12000 * Q_v / delta_h_v
        result = []
        Q_v = delta_h_v*w_v/12000
        result.append(Q_v)
        return result
    def eqn_6_4__delta_h_v(self, Q_v: float, w_v: float, **kwargs):
        # w_v = 12000 * Q_v / delta_h_v
        result = []
        delta_h_v = 12000*Q_v/w_v
        result.append(delta_h_v)
        return result
    def eqn_6_4__w_v(self, Q_v: float, delta_h_v: float, **kwargs):
        # w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return result
    @kwasak
    def eqn_6_7(self, C_1=None, C_2=None, T_1=None, T_2=None, c_p=None, delta_h_c=None, delta_h_v=None, m_b=None, m_v=None):
        return
    def eqn_6_7__C_1(self, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result
    def eqn_6_7__C_2(self, C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result
    def eqn_6_7__T_1(self, C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*m_v)/(c_p*m_b)
        result.append(T_1)
        return result
    def eqn_6_7__T_2(self, C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*m_v)/(c_p*m_b)
        result.append(T_2)
        return result
    def eqn_6_7__c_p(self, C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*m_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result
    def eqn_6_7__delta_h_c(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*m_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result
    def eqn_6_7__delta_h_v(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/m_v
        result.append(delta_h_v)
        return result
    def eqn_6_7__m_b(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_v: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_b = delta_h_v*m_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result
    def eqn_6_7__m_v(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, **kwargs):
        # m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/delta_h_v
        result.append(m_v)
        return result
    @kwasak
    def eqn_6_8(self, C_1=None, C_2=None, T_1=None, T_2=None, c_p=None, delta_h_c=None, delta_h_v=None, delta_t=None, m_b=None, w_v=None):
        return
    def eqn_6_8__C_1(self, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_1)
        return result
    def eqn_6_8__C_2(self, C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_2)
        return result
    def eqn_6_8__T_1(self, C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_1)
        return result
    def eqn_6_8__T_2(self, C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_2)
        return result
    def eqn_6_8__c_p(self, C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*delta_t*w_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return result
    def eqn_6_8__delta_h_c(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*delta_t*w_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return result
    def eqn_6_8__delta_h_v(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_t: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_t*w_v)
        result.append(delta_h_v)
        return result
    def eqn_6_8__delta_t(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_t = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*w_v)
        result.append(delta_t)
        return result
    def eqn_6_8__m_b(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, w_v: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        m_b = delta_h_v*delta_t*w_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return result
    def eqn_6_8__w_v(self, C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, delta_t: float, m_b: float, **kwargs):
        # w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        w_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*delta_t)
        result.append(w_v)
        return result
    @kwasak
    def eqn_6_9(self, A=None, dV_dt=None, delta_P=None, m=None, mu=None, r=None, r_M=None):
        return
    def eqn_6_9__A(self, dV_dt: float, delta_P: float, m: float, mu: float, r: float, r_M: float, **kwargs):
        # dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        A = (dV_dt*r_M - sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        A = (dV_dt*r_M + sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        return result
    def eqn_6_9__dV_dt(self, A: float, delta_P: float, m: float, mu: float, r: float, r_M: float, **kwargs):
        # dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
        result.append(dV_dt)
        return result
    def eqn_6_9__delta_P(self, A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float, **kwargs):
        # dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        delta_P = A*dV_dt*r_M/(A**2 - dV_dt*m*mu*r)
        result.append(delta_P)
        return result
    def eqn_6_9__m(self, A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float, **kwargs):
        # dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
        result.append(m)
        return result
    def eqn_6_9__mu(self, A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float, **kwargs):
        # dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
        result.append(mu)
        return result
    def eqn_6_9__r(self, A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float, **kwargs):
        # dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
        result.append(r)
        return result
    def eqn_6_9__r_M(self, A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float, **kwargs):
        # dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
        result.append(r_M)
        return result

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

class SelectingPump:
    @kwasak
    def eqn_8_1(self, NC=None, NS=None, SCON=None, installation_cost=None):
        """
        NS:= number ejector stages
        NC:= number of condensors
        SCON:=steam consumption based on 100-psig motive steam, lb/hr
        """
        return
    def eqn_8_1__NC(self, NS: float, SCON: float, installation_cost: float, **kwargs):
        # installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return result
    def eqn_8_1__NS(self, NC: float, SCON: float, installation_cost: float, **kwargs):
        # installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return result
    def eqn_8_1__SCON(self, NC: float, NS: float, installation_cost: float, **kwargs):
        # installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        # Solve for SCON:
        # Step 1: (SCON / 1000) ** 0.35 = installation_cost / (16000 * (NS + 2 * NC))
        # Step 2: SCON / 1000 = (installation_cost / (16000 * (NS + 2 * NC))) ** (1.0 / 0.35)
        # Step 3: SCON = 1000 * (installation_cost / (16000 * (NS + 2 * NC))) ** (1.0 / 0.35)
        SCON = 1000 * (installation_cost / (16000 * (NS + 2 * NC))) ** (1.0 / 0.35)
        return [SCON]
    def eqn_8_1__installation_cost(self, NC: float, NS: float, SCON: float, **kwargs):
        # installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        installation_cost = 1426.00150101399*SCON**(7/20)*(2.0*NC + NS)
        result.append(installation_cost)
        return result
    @kwasak
    def eqn_8_2(self, hp=None, installed_costs=None):
        """
        hp:= horse power of pump
        """
        return
    def eqn_8_2__hp(self, installed_costs: float, **kwargs):
        # installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        hp = 9.18273645546364e-9*installed_costs**2
        result.append(hp)
        return result
    def eqn_8_2__installed_costs(self, hp: float, **kwargs):
        # installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        installed_costs = 10435.5162785557*sqrt(hp)
        result.append(installed_costs)
        return result
    @kwasak
    def eqn_8_3(self, hp=None, installed_costs=None):
        return
    def eqn_8_3__hp(self, installed_costs: float, **kwargs):
        # installed_costs = 38000 * (hp / 10) ** 0.45
        # Solve for hp:
        # Step 1: (hp / 10) ** 0.45 = installed_costs / (38000)
        # Step 2: hp / 10 = (installed_costs / (38000)) ** (1.0 / 0.45)
        # Step 3: hp = 10 * (installed_costs / (38000)) ** (1.0 / 0.45)
        hp = 10 * (installed_costs / (38000)) ** (1.0 / 0.45)
        return [hp]
    def eqn_8_3__installed_costs(self, hp: float, **kwargs):
        # installed_costs = 38000 * (hp / 10) ** 0.45
        result = []
        installed_costs = 13482.9087908759*hp**(9/20)
        result.append(installed_costs)
        return result
    @kwasak
    def eqn_8_4(self, hp=None, installed_costs=None):
        return
    def eqn_8_4__hp(self, installed_costs: float, **kwargs):
        # installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        hp = -9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        hp = 9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        return result
    def eqn_8_4__installed_costs(self, hp: float, **kwargs):
        # installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        installed_costs = 10350.7864343909*hp**(2/5)
        result.append(installed_costs)
        return result
    @kwasak
    def eqn_8_5(self, Eff=None, actual_brake_horsepower=None, theoretical_adiabatic_horsepower=None):
        """
        Eff:= thermal efficiency
        """
        return
    def eqn_8_5__Eff(self, actual_brake_horsepower: float, theoretical_adiabatic_horsepower: float, **kwargs):
        # Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        Eff = theoretical_adiabatic_horsepower/actual_brake_horsepower
        result.append(Eff)
        return result
    def eqn_8_5__actual_brake_horsepower(self, Eff: float, theoretical_adiabatic_horsepower: float, **kwargs):
        # Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        actual_brake_horsepower = theoretical_adiabatic_horsepower/Eff
        result.append(actual_brake_horsepower)
        return result
    def eqn_8_5__theoretical_adiabatic_horsepower(self, Eff: float, actual_brake_horsepower: float, **kwargs):
        # Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        theoretical_adiabatic_horsepower = Eff*actual_brake_horsepower
        result.append(theoretical_adiabatic_horsepower)
        return result
    @kwasak
    def eqn_8_6(self, M=None, P_1=None, P_2=None, R=None, T=None, adiabatic_hp=None, k=None, w=None):
        """
        deg_R:=absolute temperature
        M:=molecular weight
        R:=gas constant, 1544 ft*lb_f / (lb*mol) * deg_R
        T:= absolute temperature, deg_R
        P:= absolute pressure, torr
        """
        return
    def eqn_8_6__M(self, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        M = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*adiabatic_hp*(k - 1))
        result.append(M)
        return result
    def eqn_8_6__P_1(self, M: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_1 = P_2/(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_1)
        return result
    def eqn_8_6__P_2(self, M: float, P_1: float, R: float, T: float, adiabatic_hp: float, k: float, w: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_2 = P_1*(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_2)
        return result
    def eqn_8_6__R(self, M: float, P_1: float, P_2: float, T: float, adiabatic_hp: float, k: float, w: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        R = 1980000*M*adiabatic_hp*(k - 1)/(T*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(R)
        return result
    def eqn_8_6__T(self, M: float, P_1: float, P_2: float, R: float, adiabatic_hp: float, k: float, w: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        T = 1980000*M*adiabatic_hp*(k - 1)/(R*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(T)
        return result
    def eqn_8_6__adiabatic_hp(self, M: float, P_1: float, P_2: float, R: float, T: float, k: float, w: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        adiabatic_hp = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*M*(k - 1))
        result.append(adiabatic_hp)
        return result
    def eqn_8_6__k(self, M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, w: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        # k appears in the exponent — use numerical solver
        from scipy.optimize import brentq
        def _res(k_val):
            return (k_val / (k_val - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k_val - 1) / k_val) - 1)) - adiabatic_hp
        lo, hi = None, None
        prev = _res(1.01)
        for i in range(1, 100000):
            x = 1.01 + i * 0.01
            try:
                cur = _res(x)
            except Exception:
                continue
            if prev * cur < 0:
                lo, hi = x - 0.01, x
                break
            prev = cur
        if lo is None:
            raise UnsolvedException("No sign change found for k")
        k = brentq(_res, lo, hi)
        return [k]
    def eqn_8_6__w(self, M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, **kwargs):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        w = 1980000*M*adiabatic_hp*(k - 1)/(R*T*k*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(w)
        return result
    @kwasak
    def eqn_8_7(self, P_1=None, P_2=None, adiabatic_hp=None, w=None):
        return
    def eqn_8_7__P_1(self, P_2: float, adiabatic_hp: float, w: float, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # Solve for P_1:
        # Step 1: (P_2 / P_1) ** 0.286 - 1 = adiabatic_hp / ((w / 20))
        # Step 2: (P_2 / P_1) ** 0.286 = adiabatic_hp / ((w / 20)) + 1
        # Step 3: P_1 = P_2 / (adiabatic_hp / ((w / 20)) + 1) ** (1.0 / 0.286)
        P_1 = P_2 / (adiabatic_hp / ((w / 20)) + 1) ** (1.0 / 0.286)
        return [P_1]
    def eqn_8_7__P_2(self, P_1: float, adiabatic_hp: float, w: float, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # Solve for P_2:
        # Step 1: (P_2 / P_1) ** 0.286 - 1 = adiabatic_hp / ((w / 20))
        # Step 2: (P_2 / P_1) ** 0.286 = adiabatic_hp / ((w / 20)) + 1
        # Step 3: P_2 / P_1 = (adiabatic_hp / ((w / 20)) + 1) ** (1.0 / 0.286)
        # Step 4: P_2 = P_1 * (adiabatic_hp / ((w / 20)) + 1) ** (1.0 / 0.286)
        P_2 = P_1 * (adiabatic_hp / ((w / 20)) + 1) ** (1.0 / 0.286)
        return [P_2]
    def eqn_8_7__adiabatic_hp(self, P_1: float, P_2: float, w: float, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_hp = 0.05*w*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_hp)
        return result
    def eqn_8_7__w(self, P_1: float, P_2: float, adiabatic_hp: float, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        w = 20.0*adiabatic_hp/((P_2/P_1)**0.286 - 1.0)
        result.append(w)
        return result
    @kwasak
    def eqn_8_8(self, P_1=None, P_2=None, adiabatic_power_watts=None, f=None):
        return
    def eqn_8_8__P_1(self, P_2: float, adiabatic_power_watts: float, f: float, **kwargs):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # Solve for P_1:
        # Step 1: (P_2 / P_1) ** 0.286 - 1 = adiabatic_power_watts / (f / 12)
        # Step 2: (P_2 / P_1) ** 0.286 = adiabatic_power_watts / (f / 12) + 1
        # Step 3: P_1 = P_2 / (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
        P_1 = P_2 / (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
        return [P_1]
    def eqn_8_8__P_2(self, P_1: float, adiabatic_power_watts: float, f: float, **kwargs):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # Solve for P_2:
        # Step 1: (P_2 / P_1) ** 0.286 - 1 = adiabatic_power_watts / (f / 12)
        # Step 2: (P_2 / P_1) ** 0.286 = adiabatic_power_watts / (f / 12) + 1
        # Step 3: P_2 / P_1 = (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
        # Step 4: P_2 = P_1 * (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
        P_2 = P_1 * (adiabatic_power_watts / (f / 12) + 1) ** (1.0 / 0.286)
        return [P_2]
    def eqn_8_8__adiabatic_power_watts(self, P_1: float, P_2: float, f: float, **kwargs):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_power_watts = 0.0833333333333333*f*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_power_watts)
        return result
    def eqn_8_8__f(self, P_1: float, P_2: float, adiabatic_power_watts: float, **kwargs):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        f = 12.0*adiabatic_power_watts/((P_2/P_1)**0.286 - 1.0)
        result.append(f)
        return result
    @kwasak
    def eqn_8_9(self, E_j=None, E_m=None, e=None, r=None, s=None):
        """
        E_j:=ejector thermal efficiency
        e:=electrical cost, cents per kWh
        s:=steam cost, dollar per 1000 lb
        E_m:=mechanical pump thermal efficiency
        """
        return
    def eqn_8_9__E_j(self, E_m: float, e: float, r: float, s: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_j = 0.341296928327645*E_m*r*s/e
        result.append(E_j)
        return result
    def eqn_8_9__E_m(self, E_j: float, e: float, r: float, s: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_m = 2.93*E_j*e/(r*s)
        result.append(E_m)
        return result
    def eqn_8_9__e(self, E_j: float, E_m: float, r: float, s: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        e = 0.341296928327645*E_m*r*s/E_j
        result.append(e)
        return result
    def eqn_8_9__r(self, E_j: float, E_m: float, e: float, s: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        r = 2.93*E_j*e/(E_m*s)
        result.append(r)
        return result
    def eqn_8_9__s(self, E_j: float, E_m: float, e: float, r: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        s = 2.93*E_j*e/(E_m*r)
        result.append(s)
        return result

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
        A = w_s/(rho_s*v)
        result.append(A)
        return result
    def eqn_9_1__rho_s(self, A: float, v: float, w_s: float, **kwargs):
        # w_s = v * A * rho_s
        result = []
        rho_s = w_s/(A*v)
        result.append(rho_s)
        return result
    def eqn_9_1__v(self, A: float, rho_s: float, w_s: float, **kwargs):
        # w_s = v * A * rho_s
        result = []
        v = w_s/(A*rho_s)
        result.append(v)
        return result
    def eqn_9_1__w_s(self, A: float, rho_s: float, v: float, **kwargs):
        # w_s = v * A * rho_s
        result = []
        w_s = A*rho_s*v
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
        P_m = 1.334027668054e-6*w_s**2/(d_n**4*rho_s)
        result.append(P_m)
        return result
    def eqn_9_2__d_n(self, P_m: float, rho_s: float, w_s: float, **kwargs):
        # w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        d_n = -0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        d_n = 0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        return result
    def eqn_9_2__rho_s(self, P_m: float, d_n: float, w_s: float, **kwargs):
        # w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        rho_s = 1.334027668054e-6*w_s**2/(P_m*d_n**4)
        result.append(rho_s)
        return result
    def eqn_9_2__w_s(self, P_m: float, d_n: float, rho_s: float, **kwargs):
        # w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        w_s = 865.8*d_n**2*sqrt(P_m*rho_s)
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
        P_s = 33.3333333333333*(23.0*V - 10.0*t_e*w_j)/V
        result.append(P_s)
        return result
    def eqn_9_3__V(self, P_s: float, t_e: float, w_j: float, **kwargs):
        # t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        V = -1000.0*t_e*w_j/(3.0*P_s - 2300.0)
        result.append(V)
        return result
    def eqn_9_3__t_e(self, P_s: float, V: float, w_j: float, **kwargs):
        # t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        t_e = 0.001*V*(2300.0 - 3.0*P_s)/w_j
        result.append(t_e)
        return result
    def eqn_9_3__w_j(self, P_s: float, V: float, t_e: float, **kwargs):
        # t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        w_j = 0.001*V*(2300.0 - 3.0*P_s)/t_e
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
        AEL = w_s/(SC*r)
        result.append(AEL)
        return result
    def eqn_9_4__SC(self, AEL: float, r: float, w_s: float, **kwargs):
        # w_s = AEL * r * SC
        result = []
        SC = w_s/(AEL*r)
        result.append(SC)
        return result
    def eqn_9_4__r(self, AEL: float, SC: float, w_s: float, **kwargs):
        # w_s = AEL * r * SC
        result = []
        r = w_s/(AEL*SC)
        result.append(r)
        return result
    def eqn_9_4__w_s(self, AEL: float, SC: float, r: float, **kwargs):
        # w_s = AEL * r * SC
        result = []
        w_s = AEL*SC*r
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
        V = t_h*w_h/r_h
        result.append(V)
        return result
    def eqn_9_5__r_h(self, V: float, t_h: float, w_h: float, **kwargs):
        # w_h = r_h * V / t_h
        result = []
        r_h = t_h*w_h/V
        result.append(r_h)
        return result
    def eqn_9_5__t_h(self, V: float, r_h: float, w_h: float, **kwargs):
        # w_h = r_h * V / t_h
        result = []
        t_h = V*r_h/w_h
        result.append(t_h)
        return result
    def eqn_9_5__w_h(self, V: float, r_h: float, t_h: float, **kwargs):
        # w_h = r_h * V / t_h
        result = []
        w_h = V*r_h/t_h
        result.append(w_h)
        return result

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
