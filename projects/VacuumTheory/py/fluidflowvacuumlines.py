from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


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
