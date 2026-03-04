from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class ElectricCurrent:
    @kwasak
    def eqn_19_1(self, I=None, Q=None, t=None):
        """
        I := electric current
        Q := electric charge
        t := time
        """
        return

    def eqn_19_1__I(self, Q: float, t: float, **kwargs):
        # I = Q / t
        result = []
        I = Q / t
        result.append(I)
        return result

    def eqn_19_1__Q(self, I: float, t: float, **kwargs):
        # I = Q / t
        result = []
        Q = I * t
        result.append(Q)
        return result

    def eqn_19_1__t(self, I: float, Q: float, **kwargs):
        # I = Q / t
        result = []
        t = Q / I
        result.append(t)
        return result

    @kwasak
    def eqn_19_6(self, R1=None, R2=None, Reff=None):
        """
        Reff := effective resistance
        """
        return

    def eqn_19_6__R1(self, R2: float, Reff: float, **kwargs):
        # Reff = 1 / (1 / R1 + 1 / R2 + ...)
        def _residual(R1):
            return (1 / (1 / R1 + 1 / R2 + ...)) - (Reff)

        from scipy.optimize import brentq as _brentq
        import math as _math

        def _rf(x):
            v = _residual(x)
            return v.real if isinstance(v, complex) else float(v)

        _sings = sorted(set([0]))
        _intervals = []
        if _sings[0] > 1e-12:
            _e0 = min(0.01, _sings[0] * 0.001)
            _intervals.append((1e-12, _sings[0] - _e0))
        for _i in range(len(_sings) - 1):
            _gap = _sings[_i + 1] - _sings[_i]
            _eg = min(0.01, _gap * 0.001)
            if _gap > 2 * _eg:
                _intervals.append((_sings[_i] + _eg, _sings[_i + 1] - _eg))
        _top = _sings[-1] + min(0.01, max(abs(_sings[-1]) * 0.001, 1e-10))
        for _hi in [_top + 1, _top + 10, _top * 100, _top + 1000, 1e6]:
            _intervals.append((_top, _hi))
        _roots = []
        for _lo, _hi in _intervals:
            if _lo >= _hi:
                continue
            _N = 50
            _step = (_hi - _lo) / _N
            for _j in range(_N):
                _a = _lo + _j * _step
                _b = _lo + (_j + 1) * _step
                try:
                    _fa = _rf(_a)
                    _fb = _rf(_b)
                    if _math.isfinite(_fa) and _math.isfinite(_fb) and _fa * _fb < 0:
                        _roots.append(_brentq(_rf, _a, _b))
                except Exception:
                    continue
        if _roots:
            _pos = [r for r in _roots if r > 0]
            return list(set(round(r, 10) for r in (_pos if _pos else _roots)))
        return [safe_brentq(_residual)]

    def eqn_19_6__R2(self, R1: float, Reff: float, **kwargs):
        # Reff = 1 / (1 / R1 + 1 / R2 + ...)
        def _residual(R2):
            return (1 / (1 / R1 + 1 / R2 + ...)) - (Reff)

        from scipy.optimize import brentq as _brentq
        import math as _math

        def _rf(x):
            v = _residual(x)
            return v.real if isinstance(v, complex) else float(v)

        _sings = sorted(set([0]))
        _intervals = []
        if _sings[0] > 1e-12:
            _e0 = min(0.01, _sings[0] * 0.001)
            _intervals.append((1e-12, _sings[0] - _e0))
        for _i in range(len(_sings) - 1):
            _gap = _sings[_i + 1] - _sings[_i]
            _eg = min(0.01, _gap * 0.001)
            if _gap > 2 * _eg:
                _intervals.append((_sings[_i] + _eg, _sings[_i + 1] - _eg))
        _top = _sings[-1] + min(0.01, max(abs(_sings[-1]) * 0.001, 1e-10))
        for _hi in [_top + 1, _top + 10, _top * 100, _top + 1000, 1e6]:
            _intervals.append((_top, _hi))
        _roots = []
        for _lo, _hi in _intervals:
            if _lo >= _hi:
                continue
            _N = 50
            _step = (_hi - _lo) / _N
            for _j in range(_N):
                _a = _lo + _j * _step
                _b = _lo + (_j + 1) * _step
                try:
                    _fa = _rf(_a)
                    _fb = _rf(_b)
                    if _math.isfinite(_fa) and _math.isfinite(_fb) and _fa * _fb < 0:
                        _roots.append(_brentq(_rf, _a, _b))
                except Exception:
                    continue
        if _roots:
            _pos = [r for r in _roots if r > 0]
            return list(set(round(r, 10) for r in (_pos if _pos else _roots)))
        return [safe_brentq(_residual)]

    def eqn_19_6__Reff(self, R1: float, R2: float, **kwargs):
        # Reff = 1 / (1 / R1 + 1 / R2 + ...)
        def _residual(Reff):
            return (1 / (1 / R1 + 1 / R2 + ...)) - (Reff)

        return [safe_brentq(_residual)]
