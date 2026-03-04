from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class SelectingPump:
    @kwasak
    def eqn_8_1(self, NC=None, NS=None, SCON=None, installation_cost=None):
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
        installation_cost = 1426.00150101399 * SCON ** (7 / 20) * (2.0 * NC + NS)
        result.append(installation_cost)
        return result
    @kwasak
    def eqn_8_2(self, hp=None, installed_costs=None):
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
        # k appears in the exponent — use numerical solver
        from scipy.optimize import brentq
        import numpy as np

        def _res(k_val):
            try:
                # Force complex evaluation to handle negative bases in fractional powers
                target_var_complex = complex(k_val, 0)
                val = (
                    target_var_complex
                    / (target_var_complex - 1)
                    * (w * R * T)
                    / (M * 550 * 3600)
                    * ((P_2 / P_1) ** ((target_var_complex - 1) / target_var_complex) - 1)
                ) - adiabatic_hp
                return val.real if hasattr(val, "real") else val
            except Exception:
                return float("nan")

        lo, hi = None, None
        # Expanded search: log-space from 1e-6 to 1e6 plus some linear steps
        search_points = np.logspace(-6, 6, 500)
        for i in range(len(search_points) - 1):
            p1, p2 = search_points[i], search_points[i + 1]
            r1, r2 = _res(p1), _res(p2)
            if np.isfinite(r1) and np.isfinite(r2) and r1 * r2 <= 0:
                lo, hi = p1, p2
                break
        if lo is None:
            # Fallback to a wider linear search if logspace fails
            for x in np.linspace(0.001, 10000, 1000):
                r = _res(x)
                if np.isfinite(r):
                    if lo is None:
                        lo_val, lo = r, x
                    if r * lo_val <= 0:
                        hi = x
                        break
        if lo is None or hi is None:
            raise UnsolvedException("No sign change found for k in expanded range")
        k = brentq(_res, lo, hi)
        return [k]
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
        adiabatic_power_watts = 0.0833333333333333 * f * ((P_2 / P_1) ** (143 / 500) - 1.0)
        result.append(adiabatic_power_watts)
        return result
    def eqn_8_8__f(self, P_1: float, P_2: float, adiabatic_power_watts: float, **kwargs):
        # adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        f = 12.0 * adiabatic_power_watts / ((P_2 / P_1) ** 0.286 - 1.0)
        result.append(f)
        return result
    @kwasak
    def eqn_8_9(self, E_j=None, E_m=None, e=None, r=None, s=None):
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
