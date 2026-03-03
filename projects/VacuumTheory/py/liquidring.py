from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class LiquidRing:
    @kwasak
    def eqn_10_1(self, D_r=None, sig_R=None, w=None):
        """
        sig_R := rotor tip speed ft/s
        D_r := rotor Diameter
        w := rotational speed
        """
        return

    def eqn_10_1__D_r(self, sig_R: float, w: float, **kwargs):
        # sig_R = 0.00436 * D_r * w
        result = []
        D_r = 229.357798165138 * sig_R / w
        result.append(D_r)
        return result

    def eqn_10_1__sig_R(self, D_r: float, w: float, **kwargs):
        # sig_R = 0.00436 * D_r * w
        result = []
        sig_R = 0.00436 * D_r * w
        result.append(sig_R)
        return result

    def eqn_10_1__w(self, D_r: float, sig_R: float, **kwargs):
        # sig_R = 0.00436 * D_r * w
        result = []
        w = 229.357798165138 * sig_R / D_r
        result.append(w)
        return result

    @kwasak
    def eqn_10_10(self, bhp=None, bhp_0=None, mu=None, rho=None):
        return

    def eqn_10_10__bhp(self, bhp_0: float, mu: float, rho: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp = 0.0005 * bhp_0 * (31.0 * mu ** (4 / 25) * rho ** (21 / 25) + 1000.0)
        result.append(bhp)
        return result

    def eqn_10_10__bhp_0(self, bhp: float, mu: float, rho: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp_0 = 2000.0 * bhp / (31.0 * mu**0.16 * rho**0.84 + 1000.0)
        result.append(bhp_0)
        return result

    def eqn_10_10__mu(self, bhp: float, bhp_0: float, rho: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        mu = (
            -4.7751763343393e-10
            * I
            * (2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84) ** (25 / 4)
        )
        result.append(mu)
        mu = (
            4.7751763343393e-10
            * I
            * (2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84) ** (25 / 4)
        )
        result.append(mu)
        mu = -4.7751763343393e-10 * (
            2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84
        ) ** (25 / 4)
        result.append(mu)
        mu = 4.7751763343393e-10 * (
            2000.0 * bhp / (bhp_0 * rho**0.84) - 1000.0 / rho**0.84
        ) ** (25 / 4)
        result.append(mu)
        return result

    def eqn_10_10__rho(self, bhp: float, bhp_0: float, mu: float, **kwargs):
        # bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        # Solve for rho:
        # Step 1: bhp / bhp_0 = 0.5 + 0.0155 *rho ** 0.84* mu ** 0.16
        # Step 2: (bhp / bhp_0 - 0.5) = 0.0155 * rho ** 0.84 * mu ** 0.16
        # Step 3: rho ** 0.84 = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16))
        # Step 4: rho = ((bhp / bhp_0 - 0.5) / (0.0155 * mu ** 0.16)) ** (1.0 / 0.84)
        rho = ((bhp / bhp_0 - 0.5) / (0.0155 * mu**0.16)) ** (1.0 / 0.84)
        return [rho]

    @kwasak
    def eqn_10_11(self, T_c=None, T_s=None):
        return

    def eqn_10_11__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 10
        result = []
        T_c = T_s + 10
        result.append(T_c)
        return result

    def eqn_10_11__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 10
        result = []
        T_s = T_c - 10
        result.append(T_s)
        return result

    @kwasak
    def eqn_10_12(self, T_c=None, T_s=None):
        return

    def eqn_10_12__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 5
        result = []
        T_c = T_s + 5
        result.append(T_c)
        return result

    def eqn_10_12__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 5
        result = []
        T_s = T_c - 5
        result.append(T_s)
        return result

    @kwasak
    def eqn_10_13(self, T_c=None, T_s=None):
        return

    def eqn_10_13__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 25
        result = []
        T_c = T_s + 25
        result.append(T_c)
        return result

    def eqn_10_13__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 25
        result = []
        T_s = T_c - 25
        result.append(T_s)
        return result

    @kwasak
    def eqn_10_14(self, T_c=None, T_s=None):
        return

    def eqn_10_14__T_c(self, T_s: float, **kwargs):
        # T_c = T_s + 12
        result = []
        T_c = T_s + 12
        result.append(T_c)
        return result

    def eqn_10_14__T_s(self, T_c: float, **kwargs):
        # T_c = T_s + 12
        result = []
        T_s = T_c - 12
        result.append(T_s)
        return result

    @kwasak
    def eqn_10_15(self, P=None, S_Th=None, S_p=None, p_s=None):
        return

    def eqn_10_15__P(self, S_Th: float, S_p: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        P = S_Th * p_s / (S_Th - S_p)
        result.append(P)
        return result

    def eqn_10_15__S_Th(self, P: float, S_p: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        S_Th = P * S_p / (P - p_s)
        result.append(S_Th)
        return result

    def eqn_10_15__S_p(self, P: float, S_Th: float, p_s: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        S_p = S_Th * (P - p_s) / P
        result.append(S_p)
        return result

    def eqn_10_15__p_s(self, P: float, S_Th: float, S_p: float, **kwargs):
        # S_p = S_Th * (P - p_s) / P
        result = []
        p_s = P * (S_Th - S_p) / S_Th
        result.append(p_s)
        return result

    @kwasak
    def eqn_10_16(self, P=None, S_0=None, S_Th=None, p_0=None):
        return

    def eqn_10_16__P(self, S_0: float, S_Th: float, p_0: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        P = p_0 * (S_Th / S_0) ** (5 / 3) / ((S_Th / S_0) ** 1.66666666666667 - 1.0)
        result.append(P)
        P = (
            p_0
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            / (
                (
                    -0.5 * (S_Th / S_0) ** 0.333333333333333
                    - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
                )
                ** 5
                - 1.0
            )
        )
        result.append(P)
        P = (
            p_0
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            / (
                (
                    -0.5 * (S_Th / S_0) ** 0.333333333333333
                    + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
                )
                ** 5
                - 1.0
            )
        )
        result.append(P)
        return result

    def eqn_10_16__S_0(self, P: float, S_Th: float, p_0: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th / (P / (P - p_0)) ** (3 / 5)
        result.append(S_0)
        return result

    def eqn_10_16__S_Th(self, P: float, S_0: float, p_0: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0 * (P / (P - p_0)) ** (3 / 5)
        result.append(S_Th)
        return result

    def eqn_10_16__p_0(self, P: float, S_0: float, S_Th: float, **kwargs):
        # S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        p_0 = P - P / (S_Th / S_0) ** (5 / 3)
        result.append(p_0)
        p_0 = (
            P
            - P
            / (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
        )
        result.append(p_0)
        p_0 = (
            P
            - P
            / (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
        )
        result.append(p_0)
        return result

    @kwasak
    def eqn_10_17(self, P=None, S_0=None, S_Th=None, p_0=None, p_s=None):
        return

    def eqn_10_17__P(self, S_0: float, S_Th: float, p_0: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        P = (p_0 * (S_Th / S_0) ** (5 / 3) - p_s) / (
            (S_Th / S_0) ** 1.66666666666667 - 1.0
        )
        result.append(P)
        P = (
            p_0
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            - p_s
        ) / (
            (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            - 1.0
        )
        result.append(P)
        P = (
            p_0
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            - p_s
        ) / (
            (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            - 1.0
        )
        result.append(P)
        return result

    def eqn_10_17__S_0(self, P: float, S_Th: float, p_0: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th / ((P - p_s) / (P - p_0)) ** (3 / 5)
        result.append(S_0)
        return result

    def eqn_10_17__S_Th(self, P: float, S_0: float, p_0: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0 * ((P - p_s) / (P - p_0)) ** (3 / 5)
        result.append(S_Th)
        return result

    def eqn_10_17__p_0(self, P: float, S_0: float, S_Th: float, p_s: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_0 = (P * (S_Th / S_0) ** (5 / 3) - P + p_s) / (S_Th / S_0) ** (5 / 3)
        result.append(p_0)
        p_0 = (
            P
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            - P
            + p_s
        ) / (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        ) ** 5
        result.append(p_0)
        p_0 = (
            P
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            - P
            + p_s
        ) / (
            -0.5 * (S_Th / S_0) ** 0.333333333333333
            + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
        ) ** 5
        result.append(p_0)
        return result

    def eqn_10_17__p_s(self, P: float, S_0: float, S_Th: float, p_0: float, **kwargs):
        # S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_s = -P * (S_Th / S_0) ** (5 / 3) + P + p_0 * (S_Th / S_0) ** (5 / 3)
        result.append(p_s)
        p_s = (
            -P
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            + P
            + p_0
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                - 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
        )
        result.append(p_s)
        p_s = (
            -P
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
            + P
            + p_0
            * (
                -0.5 * (S_Th / S_0) ** 0.333333333333333
                + 0.866025403784439 * I * (S_Th / S_0) ** 0.333333333333333
            )
            ** 5
        )
        result.append(p_s)
        return result

    @kwasak
    def eqn_10_18(
        self, P=None, S_Th=None, S_p=None, T_e=None, T_i=None, p_c=None, p_s=None
    ):
        """
        T_i := inlet  temperature of load
        """
        return

    def eqn_10_18__P(
        self,
        S_Th: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        P = (
            S_Th * T_i * p_s + 460 * S_Th * p_s - S_p * T_e * p_c - 460 * S_p * p_c
        ) / (S_Th * T_i + 460 * S_Th - S_p * T_e - 460 * S_p)
        result.append(P)
        return result

    def eqn_10_18__S_Th(
        self,
        P: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_Th = (
            S_p
            * (P * T_e + 460 * P - T_e * p_c - 460 * p_c)
            / (P * T_i + 460 * P - T_i * p_s - 460 * p_s)
        )
        result.append(S_Th)
        return result

    def eqn_10_18__S_p(
        self,
        P: float,
        S_Th: float,
        T_e: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_p = (
            S_Th
            * (P * T_i + 460 * P - T_i * p_s - 460 * p_s)
            / (P * T_e + 460 * P - T_e * p_c - 460 * p_c)
        )
        result.append(S_p)
        return result

    def eqn_10_18__T_e(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_e = (
            P * S_Th * T_i
            + 460 * P * S_Th
            - 460 * P * S_p
            - S_Th * T_i * p_s
            - 460 * S_Th * p_s
            + 460 * S_p * p_c
        ) / (S_p * (P - p_c))
        result.append(T_e)
        return result

    def eqn_10_18__T_i(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_e: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_i = (
            -460 * P * S_Th
            + P * S_p * T_e
            + 460 * P * S_p
            + 460 * S_Th * p_s
            - S_p * T_e * p_c
            - 460 * S_p * p_c
        ) / (S_Th * (P - p_s))
        result.append(T_i)
        return result

    def eqn_10_18__p_c(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_c = (
            -P * S_Th * T_i
            - 460 * P * S_Th
            + P * S_p * T_e
            + 460 * P * S_p
            + S_Th * T_i * p_s
            + 460 * S_Th * p_s
        ) / (S_p * (T_e + 460))
        result.append(p_c)
        return result

    def eqn_10_18__p_s(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_c: float,
        **kwargs,
    ):
        # S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_s = (
            P * S_Th * T_i
            + 460 * P * S_Th
            - P * S_p * T_e
            - 460 * P * S_p
            + S_p * T_e * p_c
            + 460 * S_p * p_c
        ) / (S_Th * (T_i + 460))
        result.append(p_s)
        return result

    @kwasak
    def eqn_10_19(
        self, P=None, S_Th=None, S_p=None, T_e=None, T_i=None, p_c=None, p_s=None
    ):
        return

    def eqn_10_19__P(
        self,
        S_Th: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        # Solve for P:
        # Step 1: (S_p / S_Th) ** (1.666667) = (P - p_s)*(460 + T_i) / ( (P - p_c)*(460 + T_e) )
        R = (S_p / (S_Th)) ** (1.666667)
        # Step 2: R * ((460 + T_e)) * (P - p_c) = ((460 + T_i)) * (P - p_s)
        # Step 3: P * (R * ((460 + T_e)) - ((460 + T_i))) = R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s
        # Step 4: P = (R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s) / (R * ((460 + T_e)) - ((460 + T_i)))
        P = (R * ((460 + T_e)) * p_c - ((460 + T_i)) * p_s) / (
            R * ((460 + T_e)) - ((460 + T_i))
        )
        return [P]

    def eqn_10_19__S_Th(
        self,
        P: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_Th = S_p / (
            (P * T_i + 460.0 * P - T_i * p_s - 460.0 * p_s)
            / (P * T_e + 460.0 * P - T_e * p_c - 460.0 * p_c)
        ) ** (3 / 5)
        result.append(S_Th)
        return result

    def eqn_10_19__S_p(
        self,
        P: float,
        S_Th: float,
        T_e: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_p = S_Th * (
            (P * T_i + 460.0 * P - T_i * p_s - 460.0 * p_s)
            / (P * T_e + 460.0 * P - T_e * p_c - 460.0 * p_c)
        ) ** (3 / 5)
        result.append(S_p)
        return result

    def eqn_10_19__T_e(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_e = (
            P * T_i
            - 460.0 * P * (S_p / S_Th) ** (5 / 3)
            + 460.0 * P
            - T_i * p_s
            + 460.0 * p_c * (S_p / S_Th) ** (5 / 3)
            - 460.0 * p_s
        ) / ((S_p / S_Th) ** (5 / 3) * (P - p_c))
        result.append(T_e)
        T_e = (
            P * T_i
            - 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0 * P
            - T_i * p_s
            + 460.0
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0 * p_s
        ) / (
            (P - p_c)
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
        )
        result.append(T_e)
        T_e = (
            P * T_i
            - 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0 * P
            - T_i * p_s
            + 460.0
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0 * p_s
        ) / (
            (P - p_c)
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
        )
        result.append(T_e)
        return result

    def eqn_10_19__T_i(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_e: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_i = (
            P * T_e * (S_p / S_Th) ** (5 / 3)
            + 460.0 * P * (S_p / S_Th) ** (5 / 3)
            - 460.0 * P
            - T_e * p_c * (S_p / S_Th) ** (5 / 3)
            - 460.0 * p_c * (S_p / S_Th) ** (5 / 3)
            + 460.0 * p_s
        ) / (P - p_s)
        result.append(T_i)
        T_i = (
            P
            * T_e
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0 * P
            - T_e
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0 * p_s
        ) / (P - p_s)
        result.append(T_i)
        T_i = (
            P
            * T_e
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0 * P
            - T_e
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0 * p_s
        ) / (P - p_s)
        result.append(T_i)
        return result

    def eqn_10_19__p_c(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_s: float,
        **kwargs,
    ):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_c = (
            P * T_e * (S_p / S_Th) ** (5 / 3)
            - P * T_i
            + 460.0 * P * (S_p / S_Th) ** (5 / 3)
            - 460.0 * P
            + T_i * p_s
            + 460.0 * p_s
        ) / ((S_p / S_Th) ** (5 / 3) * (T_e + 460.0))
        result.append(p_c)
        p_c = (
            P
            * T_e
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - P * T_i
            + 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0 * P
            + T_i * p_s
            + 460.0 * p_s
        ) / (
            (T_e + 460.0)
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
        )
        result.append(p_c)
        p_c = (
            P
            * T_e
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - P * T_i
            + 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            - 460.0 * P
            + T_i * p_s
            + 460.0 * p_s
        ) / (
            (T_e + 460.0)
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
        )
        result.append(p_c)
        return result

    def eqn_10_19__p_s(
        self,
        P: float,
        S_Th: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_c: float,
        **kwargs,
    ):
        # S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_s = (
            -P * T_e * (S_p / S_Th) ** (5 / 3)
            + P * T_i
            - 460.0 * P * (S_p / S_Th) ** (5 / 3)
            + 460.0 * P
            + T_e * p_c * (S_p / S_Th) ** (5 / 3)
            + 460.0 * p_c * (S_p / S_Th) ** (5 / 3)
        ) / (T_i + 460.0)
        result.append(p_s)
        p_s = (
            -P
            * T_e
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + P * T_i
            - 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0 * P
            + T_e
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                - 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
        ) / (T_i + 460.0)
        result.append(p_s)
        p_s = (
            -P
            * T_e
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + P * T_i
            - 460.0
            * P
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0 * P
            + T_e
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
            + 460.0
            * p_c
            * (
                -0.5 * (S_p / S_Th) ** 0.333333333333333
                + 0.866025403784439 * I * (S_p / S_Th) ** 0.333333333333333
            )
            ** 5
        ) / (T_i + 460.0)
        result.append(p_s)
        return result

    @kwasak
    def eqn_10_2(self, PS=None, Q_gas=None, V=None, dP=None, dt=None):
        return

    def eqn_10_2__PS(self, Q_gas: float, V: float, dP: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        PS = Q_gas - V * dP / dt
        result.append(PS)
        return result

    def eqn_10_2__Q_gas(self, PS: float, V: float, dP: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        Q_gas = PS + V * dP / dt
        result.append(Q_gas)
        return result

    def eqn_10_2__V(self, PS: float, Q_gas: float, dP: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        V = dt * (-PS + Q_gas) / dP
        result.append(V)
        return result

    def eqn_10_2__dP(self, PS: float, Q_gas: float, V: float, dt: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        dP = dt * (-PS + Q_gas) / V
        result.append(dP)
        return result

    def eqn_10_2__dt(self, PS: float, Q_gas: float, V: float, dP: float, **kwargs):
        # PS = - V * dP / dt + Q_gas
        result = []
        dt = -V * dP / (PS - Q_gas)
        result.append(dt)
        return result

    @kwasak
    def eqn_10_20(
        self,
        P=None,
        S_0=None,
        S_p=None,
        T_e=None,
        T_i=None,
        p_0=None,
        p_c=None,
        p_s=None,
    ):
        return

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
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # P appears 4 times — use numerical solver
        from scipy.optimize import brentq

        def _res(P_val):
            return (
                S_p
                * (
                    (P_val - p_0)
                    * (460 + T_i)
                    * (P_val - p_c)
                    / (P_val * (P_val - p_s) * (460 + T_e))
                )
                ** 0.6
                - S_0
            )

        lo, hi = None, None
        prev = _res(0.01)
        for i in range(1, 100000):
            x = i * 0.01
            try:
                cur = _res(x)
            except Exception:
                continue
            if prev * cur < 0:
                lo, hi = x - 0.01, x
                break
            prev = cur
        if lo is None:
            raise UnsolvedException("No sign change found for P")
        P = brentq(_res, lo, hi)
        return [P]

    def eqn_10_20__S_0(
        self,
        P: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_0: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_0 = S_p * (
            (
                P**2 * T_i
                + 460.0 * P**2
                - P * T_i * p_0
                - P * T_i * p_c
                - 460.0 * P * p_0
                - 460.0 * P * p_c
                + T_i * p_0 * p_c
                + 460.0 * p_0 * p_c
            )
            / (P * (P * T_e + 460.0 * P - T_e * p_s - 460.0 * p_s))
        ) ** (3 / 5)
        result.append(S_0)
        return result

    def eqn_10_20__S_p(
        self,
        P: float,
        S_0: float,
        T_e: float,
        T_i: float,
        p_0: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_p = S_0 / (
            (
                P**2 * T_i
                + 460.0 * P**2
                - P * T_i * p_0
                - P * T_i * p_c
                - 460.0 * P * p_0
                - 460.0 * P * p_c
                + T_i * p_0 * p_c
                + 460.0 * p_0 * p_c
            )
            / (P * (P * T_e + 460.0 * P - T_e * p_s - 460.0 * p_s))
        ) ** (3 / 5)
        result.append(S_p)
        return result

    def eqn_10_20__T_e(
        self,
        P: float,
        S_0: float,
        S_p: float,
        T_i: float,
        p_0: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for T_e:
        R = (S_0 / (S_p)) ** (1.666667)
        # (460 + T_e) = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (P - p_s)))
        # T_e = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (P - p_s))) - 460
        T_e = ((P - p_0) * (460 + T_i) * (P - p_c)) / (R * (P * (P - p_s))) - 460
        return [T_e]

    def eqn_10_20__T_i(
        self,
        P: float,
        S_0: float,
        S_p: float,
        T_e: float,
        p_0: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for T_i:
        R = (S_0 / (S_p)) ** (1.666667)
        # (460 + T_i) = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(P - p_c))
        # T_i = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(P - p_c)) - 460
        T_i = R * (P * (P - p_s) * (460 + T_e)) / ((P - p_0) * (P - p_c)) - 460
        return [T_i]

    def eqn_10_20__p_0(
        self,
        P: float,
        S_0: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_c: float,
        p_s: float,
        **kwargs,
    ):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for p_0:
        R = (S_0 / (S_p)) ** (1.666667)
        # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) )
        # (P - p_0) = R * (P * (P - p_s)*(460 + T_e) ) / ((460 + T_i) * (P - p_c))
        # p_0 = P - R * (P * (P - p_s)*(460 + T_e) ) / ((460 + T_i) * (P - p_c))
        p_0 = P - R * (P * (P - p_s) * (460 + T_e)) / ((460 + T_i) * (P - p_c))
        return [p_0]

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
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for p_c:
        R = (S_0 / (S_p)) ** (1.666667)
        # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) )
        # (P - p_c) = R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
        # p_c = P - R * (P * (P - p_s)*(460 + T_e) ) / ((P - p_0)*(460 + T_i))
        p_c = P - R * (P * (P - p_s) * (460 + T_e)) / ((P - p_0) * (460 + T_i))
        return [p_c]

    def eqn_10_20__p_s(
        self,
        P: float,
        S_0: float,
        S_p: float,
        T_e: float,
        T_i: float,
        p_0: float,
        p_c: float,
        **kwargs,
    ):
        # S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # Solve for p_s:
        R = (S_0 / (S_p)) ** (1.666667)
        # After clearing **0.6: R = (P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) )
        # (P - p_s) = ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
        # p_s = P - ((P - p_0)*(460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
        p_s = P - ((P - p_0) * (460 + T_i) * (P - p_c)) / (R * (P * (460 + T_e)))
        return [p_s]

    @kwasak
    def eqn_10_21(self, P=None, P_d=None, P_prime=None):
        """
        P_prime := pseudo suction pressure
        P_d := actual pump discharge pressure
        """
        return

    def eqn_10_21__P(self, P_d: float, P_prime: float, **kwargs):
        # P_prime = P / P_d * 760
        result = []
        P = P_d * P_prime / 760
        result.append(P)
        return result

    def eqn_10_21__P_d(self, P: float, P_prime: float, **kwargs):
        # P_prime = P / P_d * 760
        result = []
        P_d = 760 * P / P_prime
        result.append(P_d)
        return result

    def eqn_10_21__P_prime(self, P: float, P_d: float, **kwargs):
        # P_prime = P / P_d * 760
        result = []
        P_prime = 760 * P / P_d
        result.append(P_prime)
        return result

    @kwasak
    def eqn_10_3(self, N_mfw=None, Q_gas=None, T=None):
        return

    def eqn_10_3__N_mfw(self, Q_gas: float, T: float, **kwargs):
        # Q_gas = 9.25 * N_mfw * T
        result = []
        N_mfw = 0.108108108108108 * Q_gas / T
        result.append(N_mfw)
        return result

    def eqn_10_3__Q_gas(self, N_mfw: float, T: float, **kwargs):
        # Q_gas = 9.25 * N_mfw * T
        result = []
        Q_gas = 9.25 * N_mfw * T
        result.append(Q_gas)
        return result

    def eqn_10_3__T(self, N_mfw: float, Q_gas: float, **kwargs):
        # Q_gas = 9.25 * N_mfw * T
        result = []
        T = 0.108108108108108 * Q_gas / N_mfw
        result.append(T)
        return result

    @kwasak
    def eqn_10_4(self, Q_gas=None, SP_1=None, SP_2=None, S_p=None, V=None, t=None):
        return

    def eqn_10_4__Q_gas(
        self, SP_1: float, SP_2: float, S_p: float, V: float, t: float, **kwargs
    ):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        Q_gas = -(SP_1 - SP_2 * exp(S_p * t / V)) / (exp(S_p * t / V) - 1)
        result.append(Q_gas)
        return result

    def eqn_10_4__SP_1(
        self, Q_gas: float, SP_2: float, S_p: float, V: float, t: float, **kwargs
    ):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_1 = Q_gas + (-Q_gas + SP_2) * exp(S_p * t / V)
        result.append(SP_1)
        return result

    def eqn_10_4__SP_2(
        self, Q_gas: float, SP_1: float, S_p: float, V: float, t: float, **kwargs
    ):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_2 = (Q_gas * exp(S_p * t / V) - Q_gas + SP_1) * exp(-S_p * t / V)
        result.append(SP_2)
        return result

    def eqn_10_4__S_p(
        self, Q_gas: float, SP_1: float, SP_2: float, V: float, t: float, **kwargs
    ):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        S_p = V * log((Q_gas - SP_1) / (Q_gas - SP_2)) / t
        result.append(S_p)
        return result

    def eqn_10_4__V(
        self, Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float, **kwargs
    ):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        V = S_p * t / log((Q_gas - SP_1) / (Q_gas - SP_2))
        result.append(V)
        return result

    def eqn_10_4__t(
        self, Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float, **kwargs
    ):
        # t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        t = V * log((Q_gas - SP_1) / (Q_gas - SP_2)) / S_p
        result.append(t)
        return result

    @kwasak
    def eqn_10_5(self, P_1=None, P_2=None, S_p=None, V=None, t=None):
        return

    def eqn_10_5__P_1(self, P_2: float, S_p: float, V: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        P_1 = P_2 * exp(S_p * t / V)
        result.append(P_1)
        return result

    def eqn_10_5__P_2(self, P_1: float, S_p: float, V: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        P_2 = P_1 * exp(-S_p * t / V)
        result.append(P_2)
        return result

    def eqn_10_5__S_p(self, P_1: float, P_2: float, V: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        S_p = V * log(P_1 / P_2) / t
        result.append(S_p)
        return result

    def eqn_10_5__V(self, P_1: float, P_2: float, S_p: float, t: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        V = S_p * t / log(P_1 / P_2)
        result.append(V)
        return result

    def eqn_10_5__t(self, P_1: float, P_2: float, S_p: float, V: float, **kwargs):
        # t = V / S_p * log(P_1 / P_2)
        result = []
        t = V * log(P_1 / P_2) / S_p
        result.append(t)
        return result

    @kwasak
    def eqn_10_6(self, P_1=None, P_2=None, S_a=None, V=None, t=None):
        return

    def eqn_10_6__P_1(self, P_2: float, S_a: float, V: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        P_1 = P_2 * exp(S_a * t / V)
        result.append(P_1)
        return result

    def eqn_10_6__P_2(self, P_1: float, S_a: float, V: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        P_2 = P_1 * exp(-S_a * t / V)
        result.append(P_2)
        return result

    def eqn_10_6__S_a(self, P_1: float, P_2: float, V: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        S_a = V * log(P_1 / P_2) / t
        result.append(S_a)
        return result

    def eqn_10_6__V(self, P_1: float, P_2: float, S_a: float, t: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        V = S_a * t / log(P_1 / P_2)
        result.append(V)
        return result

    def eqn_10_6__t(self, P_1: float, P_2: float, S_a: float, V: float, **kwargs):
        # S_a = V / t * log(P_1 / P_2)
        result = []
        t = V * log(P_1 / P_2) / S_a
        result.append(t)
        return result

    @kwasak
    def eqn_10_8(
        self,
        bhp=None,
        c_p=None,
        delta_T=None,
        delta_h_i=None,
        f_a=None,
        rho=None,
        w_i=None,
    ):
        return

    def eqn_10_8__bhp(
        self,
        c_p: float,
        delta_T: float,
        delta_h_i: float,
        f_a: float,
        rho: float,
        w_i: float,
        **kwargs,
    ):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        bhp = (
            0.00315127701375246 * c_p * delta_T * f_a * rho
            - 0.000392927308447937 * delta_h_i * w_i
        )
        result.append(bhp)
        return result

    def eqn_10_8__c_p(
        self,
        bhp: float,
        delta_T: float,
        delta_h_i: float,
        f_a: float,
        rho: float,
        w_i: float,
        **kwargs,
    ):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        c_p = (
            0.124688279301746 * (2545.0 * bhp + delta_h_i * w_i) / (delta_T * f_a * rho)
        )
        result.append(c_p)
        return result

    def eqn_10_8__delta_T(
        self,
        bhp: float,
        c_p: float,
        delta_h_i: float,
        f_a: float,
        rho: float,
        w_i: float,
        **kwargs,
    ):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_T = (
            0.124688279301746 * (2545.0 * bhp + delta_h_i * w_i) / (c_p * f_a * rho)
        )
        result.append(delta_T)
        return result

    def eqn_10_8__delta_h_i(
        self,
        bhp: float,
        c_p: float,
        delta_T: float,
        f_a: float,
        rho: float,
        w_i: float,
        **kwargs,
    ):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_h_i = 0.02 * (-127250.0 * bhp + 401.0 * c_p * delta_T * f_a * rho) / w_i
        result.append(delta_h_i)
        return result

    def eqn_10_8__f_a(
        self,
        bhp: float,
        c_p: float,
        delta_T: float,
        delta_h_i: float,
        rho: float,
        w_i: float,
        **kwargs,
    ):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        f_a = (
            0.124688279301746 * (2545.0 * bhp + delta_h_i * w_i) / (c_p * delta_T * rho)
        )
        result.append(f_a)
        return result

    def eqn_10_8__rho(
        self,
        bhp: float,
        c_p: float,
        delta_T: float,
        delta_h_i: float,
        f_a: float,
        w_i: float,
        **kwargs,
    ):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        rho = (
            0.124688279301746 * (2545.0 * bhp + delta_h_i * w_i) / (c_p * delta_T * f_a)
        )
        result.append(rho)
        return result

    def eqn_10_8__w_i(
        self,
        bhp: float,
        c_p: float,
        delta_T: float,
        delta_h_i: float,
        f_a: float,
        rho: float,
        **kwargs,
    ):
        # delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        w_i = 0.02 * (-127250.0 * bhp + 401.0 * c_p * delta_T * f_a * rho) / delta_h_i
        result.append(w_i)
        return result

    @kwasak
    def eqn_10_9(self, T_c=None, T_s=None, delta_T=None):
        return

    def eqn_10_9__T_c(self, T_s: float, delta_T: float, **kwargs):
        # T_c = T_s + delta_T
        result = []
        T_c = T_s + delta_T
        result.append(T_c)
        return result

    def eqn_10_9__T_s(self, T_c: float, delta_T: float, **kwargs):
        # T_c = T_s + delta_T
        result = []
        T_s = T_c - delta_T
        result.append(T_s)
        return result

    def eqn_10_9__delta_T(self, T_c: float, T_s: float, **kwargs):
        # T_c = T_s + delta_T
        result = []
        delta_T = T_c - T_s
        result.append(delta_T)
        return result
