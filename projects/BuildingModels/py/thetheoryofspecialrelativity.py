from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


class TheTheoryOfSpecialRelativity:
    @kwasak
    def eqn_24_11(self, E=None, c=None, m_0=None, p=None):
        """
        E := energy
        p := momentum
        c := speed of light
        """
        return

    def eqn_24_11__E(self, c: float, m_0: float, p: float, **kwargs):
        # E ** 2 = p ** 2 * c ** 2 + m_0 ** 2
        result = []
        E = -sqrt(c**2 * p**2 + m_0**2)
        result.append(E)
        E = sqrt(c**2 * p**2 + m_0**2)
        result.append(E)
        return result

    def eqn_24_11__c(self, E: float, m_0: float, p: float, **kwargs):
        # E ** 2 = p ** 2 * c ** 2 + m_0 ** 2
        result = []
        c = -sqrt((E - m_0) * (E + m_0)) / p
        result.append(c)
        c = sqrt((E - m_0) * (E + m_0)) / p
        result.append(c)
        return result

    def eqn_24_11__m_0(self, E: float, c: float, p: float, **kwargs):
        # E ** 2 = p ** 2 * c ** 2 + m_0 ** 2
        result = []
        m_0 = -sqrt((E - c * p) * (E + c * p))
        result.append(m_0)
        m_0 = sqrt((E - c * p) * (E + c * p))
        result.append(m_0)
        return result

    def eqn_24_11__p(self, E: float, c: float, m_0: float, **kwargs):
        # E ** 2 = p ** 2 * c ** 2 + m_0 ** 2
        result = []
        p = -sqrt((E - m_0) * (E + m_0)) / c
        result.append(p)
        p = sqrt((E - m_0) * (E + m_0)) / c
        result.append(p)
        return result

    @kwasak
    def eqn_24_2(self, c=None, u_0=None, u_x=None, v=None, x=None):
        """
        FE := force per unit length
        l := length
        x_0 := initial position
        v := velocity
        t := time
        x := position
        x_0 := initial position
        v_0 := initial velocity
        c := speed of light
        x_0 := initial position
        v_0 := initial velocity
        t := time
        x := position
        x_0 := initial position
        v_0 := initial velocity
        c := speed of light
        u_0 := initial velocity
        v := velocity
        x := position
        t := time
        c := speed of light
        u_0 := initial velocity
        v := velocity
        x := position
        t := time
        c := speed of light
        """
        return

    def eqn_24_2__FE(self, l: float, r: float, **kwargs):
        # FE = 2 * l / (2 * r)
        result = []
        FE = l / r
        result.append(FE)
        return result

    def eqn_24_2__c(self, t: float, v: float, x: float, **kwargs):
        # x = 0 * (x + v * t) / (t * v / c ** 2)
        def _residual(c):
            return (0 * (x + v * t) / (t * v / c**2)) - (x)

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

    def eqn_24_2__l(self, FE: float, r: float, **kwargs):
        # FE = 2 * l / (2 * r)
        result = []
        l = FE * r
        result.append(l)
        return result

    def eqn_24_2__r(self, FE: float, l: float, **kwargs):
        # FE = 2 * l / (2 * r)
        result = []
        r = l / FE
        result.append(r)
        return result

    def eqn_24_2__t(self, c: float, v: float, x: float, **kwargs):
        # x = 0 * (x + v * t) / (t * v / c ** 2)
        # Placeholder for numerical solver
        raise UnsolvedException("Pending LLM/Manual Repair")

    def eqn_24_2__u_0(self, c: float, u_x: float, v: float, x: float, **kwargs):
        # u_x = u_0 * x / (1 - v * u_x / c ** 2)
        result = []
        u_0 = u_x * (c**2 - u_x * v) / (c**2 * x)
        result.append(u_0)
        return result

    def eqn_24_2__u_x(self, c: float, u_0: float, v: float, x: float, **kwargs):
        # u_x = u_0 * x / (1 - v * u_x / c ** 2)
        result = []
        u_x = c * (c - sqrt(c**2 - 4 * u_0 * v * x)) / (2 * v)
        result.append(u_x)
        u_x = c * (c + sqrt(c**2 - 4 * u_0 * v * x)) / (2 * v)
        result.append(u_x)
        return result

    def eqn_24_2__v(self, c: float, t: float, x: float, **kwargs):
        # x = 0 * (x + v * t) / (t * v / c ** 2)
        def _residual(v):
            return (0 * (x + v * t) / (t * v / c**2)) - (x)

        return [safe_brentq(_residual)]

    def eqn_24_2__vx(self, c: float, t: float, v: float, x: float, **kwargs):
        # x = 0 * (x + v * t) + vx / c ** 2
        result = []
        vx = c**2 * x
        result.append(vx)
        return result

    def eqn_24_2__x(self, c: float, t: float, v: float, **kwargs):
        # x = 0 * (x + v * t) / (t * v / c ** 2)
        result = []
        x = 0
        result.append(x)
        return result

    @kwasak
    def eqn_24_6(self, c=None, u=None, ux=None, v=None):
        """
        ux := relative velocity
        u := initial velocity
        v := velocity of the train
        c := speed of light
        """
        return

    def eqn_24_6__c(self, u: float, ux: float, v: float, **kwargs):
        # ux = u + v / (1 + (v * u) / c ** 2)
        result = []
        c = sqrt(u * v * (-u + ux) / (u - ux + v))
        result.append(c)
        c = -sqrt(-u * v * (u - ux) / (u - ux + v))
        result.append(c)
        return result

    def eqn_24_6__u(self, c: float, ux: float, v: float, **kwargs):
        # ux = u + v / (1 + (v * u) / c ** 2)
        result = []
        u = (
            -(c**2)
            + ux * v
            - sqrt((c**2 - 2 * c * v + ux * v) * (c**2 + 2 * c * v + ux * v))
        ) / (2 * v)
        result.append(u)
        u = (
            -(c**2)
            + ux * v
            + sqrt((c**2 - 2 * c * v + ux * v) * (c**2 + 2 * c * v + ux * v))
        ) / (2 * v)
        result.append(u)
        return result

    def eqn_24_6__ux(self, c: float, u: float, v: float, **kwargs):
        # ux = u + v / (1 + (v * u) / c ** 2)
        result = []
        ux = (c**2 * u + c**2 * v + u**2 * v) / (c**2 + u * v)
        result.append(ux)
        return result

    def eqn_24_6__v(self, c: float, u: float, ux: float, **kwargs):
        # ux = u + v / (1 + (v * u) / c ** 2)
        result = []
        v = c**2 * (-u + ux) / (c**2 + u**2 - u * ux)
        result.append(v)
        return result

    @kwasak
    def eqn_24_9(self, K=None, c=None, m_0=None, u=None):
        """
        m_0 := rest mass
        c := speed of light
        """
        return

    def eqn_24_9__K(self, c: float, m_0: float, u: float, **kwargs):
        # K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2
        result = []
        K = m_0 * (c**2 - u**2)
        result.append(K)
        return result

    def eqn_24_9__c(self, K: float, m_0: float, u: float, **kwargs):
        # K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2
        result = []
        c = -sqrt((K + m_0 * u**2) / m_0)
        result.append(c)
        c = sqrt((K + m_0 * u**2) / m_0)
        result.append(c)
        return result

    def eqn_24_9__m_0(self, K: float, c: float, u: float, **kwargs):
        # K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2
        result = []
        m_0 = K / (c**2 - u**2)
        result.append(m_0)
        return result

    def eqn_24_9__u(self, K: float, c: float, m_0: float, **kwargs):
        # K = (1 - u ** 2 / c ** 2) * m_0 * c ** 2
        result = []
        u = -sqrt(-K / m_0 + c**2)
        result.append(u)
        u = sqrt(-K / m_0 + c**2)
        result.append(u)
        return result
