from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException, safe_brentq
import numpy as np


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
        NC = -0.5 * NS + 0.000350630766969363 * installation_cost / SCON ** (7 / 20)
        result.append(NC)
        return result

    def eqn_8_1__NS(self, NC: float, SCON: float, installation_cost: float, **kwargs):
        # installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0 * NC + 0.000701261533938727 * installation_cost / SCON ** (7 / 20)
        result.append(NS)
        return result

    def eqn_8_1__SCON(self, NC, NS, installation_cost, **kwargs):
        # installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        def _residual(SCON):
            return (16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35) - (installation_cost)

        return [safe_brentq(_residual)]

    def eqn_8_1__installation_cost(self, NC: float, NS: float, SCON: float, **kwargs):
        # installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        installation_cost = 1426.00150101399 * SCON ** (7 / 20) * (2.0 * NC + NS)
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
        hp = 9.18273645546364e-9 * installed_costs**2
        result.append(hp)
        return result

    def eqn_8_2__installed_costs(self, hp: float, **kwargs):
        # installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        installed_costs = 10435.5162785557 * sqrt(hp)
        result.append(installed_costs)
        return result

    @kwasak
    def eqn_8_3(self, hp=None, installed_costs=None):
        return

    def eqn_8_3__hp(self, installed_costs, **kwargs):
        # installed_costs = 38000 * (hp / 10) ** 0.45
        result = []
        hp = 6.64818534458494e-10 * installed_costs ** (20 / 9)
        result.append(hp)
        hp = (
            -0.326678815618226 * installed_costs**0.111111111111111
            - 0.118901365050371 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        hp = (
            -0.326678815618226 * installed_costs**0.111111111111111
            + 0.118901365050371 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        hp = (
            -0.173822167159837 * installed_costs**0.111111111111111
            - 0.301068825002567 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        hp = (
            -0.173822167159837 * installed_costs**0.111111111111111
            + 0.301068825002567 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        hp = (
            0.0603678051308443 * installed_costs**0.111111111111111
            - 0.342362835728782 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        hp = (
            0.0603678051308443 * installed_costs**0.111111111111111
            + 0.342362835728782 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        hp = (
            0.266311010487382 * installed_costs**0.111111111111111
            - 0.223461470678411 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        hp = (
            0.266311010487382 * installed_costs**0.111111111111111
            + 0.223461470678411 * I * installed_costs**0.111111111111111
        ) ** 20
        result.append(hp)
        return result

    def eqn_8_3__installed_costs(self, hp: float, **kwargs):
        # installed_costs = 38000 * (hp / 10) ** 0.45
        result = []
        installed_costs = 13482.9087908759 * hp ** (9 / 20)
        result.append(installed_costs)
        return result

    @kwasak
    def eqn_8_4(self, hp=None, installed_costs=None):
        return

    def eqn_8_4__hp(self, installed_costs: float, **kwargs):
        # installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        hp = -9.1741667595569e-11 * installed_costs ** (5 / 2)
        result.append(hp)
        hp = 9.1741667595569e-11 * installed_costs ** (5 / 2)
        result.append(hp)
        return result

    def eqn_8_4__installed_costs(self, hp: float, **kwargs):
        # installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        installed_costs = 10350.7864343909 * hp ** (2 / 5)
        result.append(installed_costs)
        return result

    @kwasak
    def eqn_8_5(
        self,
        Eff=None,
        actual_brake_horsepower=None,
        theoretical_adiabatic_horsepower=None,
    ):
        """
        Eff:= thermal efficiency
        """
        return

    def eqn_8_5__Eff(
        self,
        actual_brake_horsepower: float,
        theoretical_adiabatic_horsepower: float,
        **kwargs,
    ):
        # Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result.append(Eff)
        return result

    def eqn_8_5__actual_brake_horsepower(
        self, Eff: float, theoretical_adiabatic_horsepower: float, **kwargs
    ):
        # Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        actual_brake_horsepower = theoretical_adiabatic_horsepower / Eff
        result.append(actual_brake_horsepower)
        return result

    def eqn_8_5__theoretical_adiabatic_horsepower(
        self, Eff: float, actual_brake_horsepower: float, **kwargs
    ):
        # Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        theoretical_adiabatic_horsepower = Eff * actual_brake_horsepower
        result.append(theoretical_adiabatic_horsepower)
        return result

    @kwasak
    def eqn_8_6(
        self,
        M=None,
        P_1=None,
        P_2=None,
        R=None,
        T=None,
        adiabatic_hp=None,
        k=None,
        w=None,
    ):
        """
        deg_R:=absolute temperature
        M:=molecular weight
        R:=gas constant, 1544 ft*lb_f / (lb*mol) * deg_R
        T:= absolute temperature, deg_R
        P:= absolute pressure, torr
        """
        return

    def eqn_8_6__M(
        self,
        P_1: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        M = (
            R
            * T
            * k
            * w
            * ((P_2 / P_1) ** ((k - 1) / k) - 1)
            / (1980000 * adiabatic_hp * (k - 1))
        )
        result.append(M)
        return result

    def eqn_8_6__P_1(
        self,
        M: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_1 = P_2 / (
            1980000 * M * adiabatic_hp / (R * T * w)
            - 1980000 * M * adiabatic_hp / (R * T * k * w)
            + 1
        ) ** (k / (k - 1))
        result.append(P_1)
        return result

    def eqn_8_6__P_2(
        self,
        M: float,
        P_1: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_2 = P_1 * (
            1980000 * M * adiabatic_hp / (R * T * w)
            - 1980000 * M * adiabatic_hp / (R * T * k * w)
            + 1
        ) ** (k / (k - 1))
        result.append(P_2)
        return result

    def eqn_8_6__R(
        self,
        M: float,
        P_1: float,
        P_2: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        R = (
            1980000
            * M
            * adiabatic_hp
            * (k - 1)
            / (T * k * w * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        )
        result.append(R)
        return result

    def eqn_8_6__T(
        self,
        M: float,
        P_1: float,
        P_2: float,
        R: float,
        adiabatic_hp: float,
        k: float,
        w: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        T = (
            1980000
            * M
            * adiabatic_hp
            * (k - 1)
            / (R * k * w * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        )
        result.append(T)
        return result

    def eqn_8_6__adiabatic_hp(
        self,
        M: float,
        P_1: float,
        P_2: float,
        R: float,
        T: float,
        k: float,
        w: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        adiabatic_hp = (
            R * T * k * w * ((P_2 / P_1) ** ((k - 1) / k) - 1) / (1980000 * M * (k - 1))
        )
        result.append(adiabatic_hp)
        return result

    def eqn_8_6__k(
        self,
        M: float,
        P_1: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        w: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        def _residual(k):
            return (
                (
                    k
                    / (k - 1)
                    * (w * R * T)
                    / (M * 550 * 3600)
                    * ((P_2 / P_1) ** ((k - 1) / k) - 1)
                )
            ) - (adiabatic_hp)

        lo = max(0, 1, 1) + 0.01
        hi_candidates = [lo + 0.5, lo + 2, lo + 10, lo * 100, lo + 1000, 1e6]
        from scipy.optimize import brentq as _brentq
        import math as _math

        for hi in hi_candidates:
            try:
                fa = _residual(lo)
                fb = _residual(hi)
                if isinstance(fa, complex):
                    fa = fa.real
                if isinstance(fb, complex):
                    fb = fb.real
                if _math.isfinite(fa) and _math.isfinite(fb) and fa * fb < 0:

                    def _rf(x):
                        v = _residual(x)
                        return v.real if isinstance(v, complex) else float(v)

                    return [_brentq(_rf, lo, hi)]
            except Exception:
                continue
        for _lo2, _hi2 in [
            (lo, lo * 1e3),
            (lo, lo * 1e6),
            (lo, lo * 1e8),
            (lo * 0.5 + 0.001, lo * 1e4),
        ]:
            try:
                _fa2 = _residual(_lo2)
                _fb2 = _residual(_hi2)
                if isinstance(_fa2, complex):
                    _fa2 = _fa2.real
                if isinstance(_fb2, complex):
                    _fb2 = _fb2.real
                if _math.isfinite(_fa2) and _math.isfinite(_fb2) and _fa2 * _fb2 < 0:

                    def _rf2(x):
                        v = _residual(x)
                        return v.real if isinstance(v, complex) else float(v)

                    return [_brentq(_rf2, _lo2, _hi2)]
            except Exception:
                continue
        return [safe_brentq(_residual)]

    def eqn_8_6__w(
        self,
        M: float,
        P_1: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        **kwargs,
    ):
        # adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        w = (
            1980000
            * M
            * adiabatic_hp
            * (k - 1)
            / (R * T * k * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        )
        result.append(w)
        return result

    @kwasak
    def eqn_8_7(self, P_1=None, P_2=None, adiabatic_hp=None, w=None):
        return

    def eqn_8_7__P_1(self, P_2, adiabatic_hp, w, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        def _residual(P_1):
            return ((w / 20) * ((P_2 / P_1) ** 0.286 - 1)) - (adiabatic_hp)

        return [safe_brentq(_residual)]

    def eqn_8_7__P_2(self, P_1: float, adiabatic_hp: float, w: float, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        def _residual(P_2):
            return ((w / 20) * ((P_2 / P_1) ** 0.286 - 1)) - (adiabatic_hp)

        return [safe_brentq(_residual)]

    def eqn_8_7__adiabatic_hp(self, P_1: float, P_2: float, w: float, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_hp = 0.05 * w * ((P_2 / P_1) ** (143 / 500) - 1.0)
        result.append(adiabatic_hp)
        return result

    def eqn_8_7__w(self, P_1: float, P_2: float, adiabatic_hp: float, **kwargs):
        # adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        w = 20.0 * adiabatic_hp / ((P_2 / P_1) ** 0.286 - 1.0)
        result.append(w)
        return result

    @kwasak
    def eqn_8_8(self, P_1=None, P_2=None, adiabatic_power_watts=None, f=None):
        return

    def eqn_8_8__P_1(self, P_2, adiabatic_power_watts, f, **kwargs):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        def _residual(P_1):
            return (f / 12 * ((P_2 / P_1) ** 0.286 - 1)) - (adiabatic_power_watts)

        return [safe_brentq(_residual)]

    def eqn_8_8__P_2(
        self, P_1: float, adiabatic_power_watts: float, f: float, **kwargs
    ):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        def _residual(P_2):
            return (f / 12 * ((P_2 / P_1) ** 0.286 - 1)) - (adiabatic_power_watts)

        return [safe_brentq(_residual)]

    def eqn_8_8__adiabatic_power_watts(
        self, P_1: float, P_2: float, f: float, **kwargs
    ):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_power_watts = (
            0.0833333333333333 * f * ((P_2 / P_1) ** (143 / 500) - 1.0)
        )
        result.append(adiabatic_power_watts)
        return result

    def eqn_8_8__f(
        self, P_1: float, P_2: float, adiabatic_power_watts: float, **kwargs
    ):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        f = 12.0 * adiabatic_power_watts / ((P_2 / P_1) ** 0.286 - 1.0)
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
        E_j = 0.341296928327645 * E_m * r * s / e
        result.append(E_j)
        return result

    def eqn_8_9__E_m(self, E_j: float, e: float, r: float, s: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_m = 2.93 * E_j * e / (r * s)
        result.append(E_m)
        return result

    def eqn_8_9__e(self, E_j: float, E_m: float, r: float, s: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        e = 0.341296928327645 * E_m * r * s / E_j
        result.append(e)
        return result

    def eqn_8_9__r(self, E_j: float, E_m: float, e: float, s: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        r = 2.93 * E_j * e / (E_m * s)
        result.append(r)
        return result

    def eqn_8_9__s(self, E_j: float, E_m: float, e: float, r: float, **kwargs):
        # r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        s = 2.93 * E_j * e / (E_m * r)
        result.append(s)
        return result
