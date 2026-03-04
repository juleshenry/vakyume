from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

def eqn_10_19__P(self, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float, **kwargs):
    # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
    result = []
    try:
        _v = (T_e*p_c*(S_p/S_Th)**(5/3) - T_i*p_s + 460.0*p_c*(S_p/S_Th)**(5/3) - 460.0*p_s)/(T_e*(S_p/S_Th)**1.66666666666667 - T_i + 460.0*(S_p/S_Th)**1.66666666666667 - 460.0)
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6: result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError): pass
    except Exception: pass
    try:
        _v = (T_e*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - T_i*p_s + 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/(T_e*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - T_i + 460.0*(-0.5*(S_p/S_Th)**0.333333333333333 - 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0)
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6: result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError): pass
    except Exception: pass
    try:
        _v = (T_e*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - T_i*p_s + 460.0*p_c*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/(T_e*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - T_i + 460.0*(-0.5*(S_p/S_Th)**0.333333333333333 + 0.866025403784439*I*(S_p/S_Th)**0.333333333333333)**5 - 460.0)
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6: result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError): pass
    except Exception: pass
    if not result:
        def _residual(P):
            return (S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6) - (S_p)
        result.append(safe_brentq(_residual))
    return result
