from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_20__T_e(self, P: float, S_0: float, S_p: float, T_i: float, p_0: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    try:
        _v = (P**2*T_i - 460.0*P**2*(S_0/S_p)**(5/3) + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + 460.0*P*p_s*(S_0/S_p)**(5/3) + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(S_0/S_p)**(5/3)*(P - p_s))
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6: result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError): pass
    except Exception: pass
    try:
        _v = (P**2*T_i - 460.0*P**2*(-0.5*(S_0/S_p)**0.333333333333333 - 0.866025403784439*I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + 460.0*P*p_s*(-0.5*(S_0/S_p)**0.333333333333333 - 0.866025403784439*I*(S_0/S_p)**0.333333333333333)**5 + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P - p_s)*(-0.5*(S_0/S_p)**0.333333333333333 - 0.866025403784439*I*(S_0/S_p)**0.333333333333333)**5)
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6: result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError): pass
    except Exception: pass
    try:
        _v = (P**2*T_i - 460.0*P**2*(-0.5*(S_0/S_p)**0.333333333333333 + 0.866025403784439*I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + 460.0*P*p_s*(-0.5*(S_0/S_p)**0.333333333333333 + 0.866025403784439*I*(S_0/S_p)**0.333333333333333)**5 + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P - p_s)*(-0.5*(S_0/S_p)**0.333333333333333 + 0.866025403784439*I*(S_0/S_p)**0.333333333333333)**5)
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6: result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError): pass
    except Exception: pass
    if not result:
        def _residual(T_e):
            return (S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6) - (S_0)
        result.append(safe_brentq(_residual))
    return result
