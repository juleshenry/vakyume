from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


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
