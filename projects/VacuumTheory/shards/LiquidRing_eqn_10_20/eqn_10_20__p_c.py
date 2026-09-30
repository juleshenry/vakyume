from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq


def eqn_10_20__p_c(
    self,
    P: float,
    S_0: float,
    S_p: float,
    T_e: float,
    T_i: float,
    p_0: float,
    p_s: float,
    **kwargs,
):
    # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
    result = []
    try:
        _v = (
            P
            * (
                -P * T_e * (S_0 / S_p) ** (5 / 3)
                + P * T_i
                - 460.0 * P * (S_0 / S_p) ** (5 / 3)
                + 460.0 * P
                + T_e * p_s * (S_0 / S_p) ** (5 / 3)
                - T_i * p_0
                - 460.0 * p_0
                + 460.0 * p_s * (S_0 / S_p) ** (5 / 3)
            )
            / (P * T_i + 460.0 * P - T_i * p_0 - 460.0 * p_0)
        )
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6:
                result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError):
                pass
    except Exception:
        pass
    try:
        _v = (
            P
            * (
                -P
                * T_e
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
                + P * T_i
                - 460.0
                * P
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
                + 460.0 * P
                + T_e
                * p_s
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
                - T_i * p_0
                - 460.0 * p_0
                + 460.0
                * p_s
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    - 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
            )
            / (P * T_i + 460.0 * P - T_i * p_0 - 460.0 * p_0)
        )
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6:
                result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError):
                pass
    except Exception:
        pass
    try:
        _v = (
            P
            * (
                -P
                * T_e
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
                + P * T_i
                - 460.0
                * P
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
                + 460.0 * P
                + T_e
                * p_s
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
                - T_i * p_0
                - 460.0 * p_0
                + 460.0
                * p_s
                * (
                    -0.5 * (S_0 / S_p) ** 0.333333333333333
                    + 0.866025403784439 * I * (S_0 / S_p) ** 0.333333333333333
                )
                ** 5
            )
            / (P * T_i + 460.0 * P - T_i * p_0 - 460.0 * p_0)
        )
        if isinstance(_v, complex):
            if abs(_v.imag) < 1e-6:
                result.append(_v.real)
        else:
            try:
                _fv = float(_v)
                result.append(_fv)
            except (TypeError, ValueError):
                pass
    except Exception:
        pass
    if not result:

        def _residual(p_c):
            return (
                S_p
                * ((P - p_0) * (460 + T_i) * (P - p_c) / (P * (P - p_s) * (460 + T_e)))
                ** 0.6
            ) - (S_0)

        result.append(safe_brentq(_residual))
    return result
