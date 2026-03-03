from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException


def eqn_10_20__P(
    self,
    S_0: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_0: float,
    p_c: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    # After raising to 5/3: K = (P-p_0)*(P-p_c) / (P*(P-p_s))
    # where K = (S_0/S_p)^(5/3) * (460+T_e)/(460+T_i)
    # Expand: K*P*(P-p_s) = (P-p_0)*(P-p_c)
    # K*P^2 - K*p_s*P = P^2 - (p_0+p_c)*P + p_0*p_c
    # (K-1)*P^2 + (p_0+p_c - K*p_s)*P - p_0*p_c = 0  (note sign: -(K*p_s) - (-(p_0+p_c)) )
    # Wait, let me redo carefully:
    # K*P^2 - K*p_s*P = P^2 - p_0*P - p_c*P + p_0*p_c
    # (K-1)*P^2 + (-K*p_s + p_0 + p_c)*P - p_0*p_c = 0
    K = (S_0 / S_p) ** (5.0 / 3.0) * (460.0 + T_e) / (460.0 + T_i)
    a = K - 1.0
    b = -K * p_s + p_0 + p_c
    c = -p_0 * p_c
    disc = b**2 - 4.0 * a * c
    P_1 = (-b + sqrt(disc)) / (2.0 * a)
    P_2 = (-b - sqrt(disc)) / (2.0 * a)
    return [P_1, P_2]
