from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class Precondensors:
    @kwasak
    def eqn_7_1(self, P=None, p_i=None, y_i=None):
        return
    def eqn_7_1__P(self, p_i: float, y_i: float, **kwargs):
        # y_i = p_i / P
        result = []
        P = p_i / y_i
        result.append(P)
        return result
    def eqn_7_1__p_i(self, P: float, y_i: float, **kwargs):
        # y_i = p_i / P
        result = []
        p_i = P * y_i
        result.append(p_i)
        return result
    def eqn_7_1__y_i(self, P: float, p_i: float, **kwargs):
        # y_i = p_i / P
        result = []
        y_i = p_i / P
        result.append(y_i)
        return result
    @kwasak
    def eqn_7_10(self, L_c_P=None, Q_condensor_heat_duty=None, del_T=None):
        return
    def eqn_7_10__L_c_P(self, Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result.append(L_c_P)
        return result
    def eqn_7_10__Q_condensor_heat_duty(self, L_c_P: float, del_T: float, **kwargs):
        # L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500 * L_c_P * del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_10__del_T(self, L_c_P: float, Q_condensor_heat_duty: float, **kwargs):
        # L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty / (500 * L_c_P)
        result.append(del_T)
        return result
    @kwasak
    def eqn_7_11(self, Q_condensor_heat_duty=None, U_v=None, V_c=None, del_T_LM=None):
        return
    def eqn_7_11__Q_condensor_heat_duty(
        self, U_v: float, V_c: float, del_T_LM: float, **kwargs
    ):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v * V_c * del_T_LM
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_11__U_v(
        self, Q_condensor_heat_duty: float, V_c: float, del_T_LM: float, **kwargs
    ):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty / (V_c * del_T_LM)
        result.append(U_v)
        return result
    def eqn_7_11__V_c(
        self, Q_condensor_heat_duty: float, U_v: float, del_T_LM: float, **kwargs
    ):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result.append(V_c)
        return result
    def eqn_7_11__del_T_LM(
        self, Q_condensor_heat_duty: float, U_v: float, V_c: float, **kwargs
    ):
        # V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty / (U_v * V_c)
        result.append(del_T_LM)
        return result
    @kwasak
    def eqn_7_12(self, A=None, Q_condensor_heat_duty=None, U=None, del_T=None):
        return
    def eqn_7_12__A(self, Q_condensor_heat_duty: float, U: float, del_T: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty / (U * del_T)
        result.append(A)
        return result
    def eqn_7_12__Q_condensor_heat_duty(self, A: float, U: float, del_T: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A * U * del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_12__U(self, A: float, Q_condensor_heat_duty: float, del_T: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty / (A * del_T)
        result.append(U)
        return result
    def eqn_7_12__del_T(self, A: float, Q_condensor_heat_duty: float, U: float, **kwargs):
        # Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty / (A * U)
        result.append(del_T)
        return result
    @kwasak
    def eqn_7_14a(self, A=None, Q_condensor_heat_duty=None, U=None, del_T_LM=None):
        return
    def eqn_7_14a__A(
        self, Q_condensor_heat_duty: float, U: float, del_T_LM: float, **kwargs
    ):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty / (U * del_T_LM)
        result.append(A)
        return result
    def eqn_7_14a__Q_condensor_heat_duty(
        self, A: float, U: float, del_T_LM: float, **kwargs
    ):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A * U * del_T_LM
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_14a__U(
        self, A: float, Q_condensor_heat_duty: float, del_T_LM: float, **kwargs
    ):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty / (A * del_T_LM)
        result.append(U)
        return result
    def eqn_7_14a__del_T_LM(
        self, A: float, Q_condensor_heat_duty: float, U: float, **kwargs
    ):
        # A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty / (A * U)
        result.append(del_T_LM)
        return result
    @kwasak
    def eqn_7_14b(
        self, A=None, Q_condensor_heat_duty=None, U=None, del_T_1=None, del_T_2=None
    ):
        return
    def eqn_7_14b__A(
        self,
        Q_condensor_heat_duty: float,
        U: float,
        del_T_1: float,
        del_T_2: float,
        **kwargs,
    ):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty / (U * (del_T_1 - del_T_2) * log(del_T_1 - del_T_2))
        result.append(A)
        return result
    def eqn_7_14b__Q_condensor_heat_duty(
        self, A: float, U: float, del_T_1: float, del_T_2: float, **kwargs
    ):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A * U * (del_T_1 - del_T_2) * log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_14b__U(
        self,
        A: float,
        Q_condensor_heat_duty: float,
        del_T_1: float,
        del_T_2: float,
        **kwargs,
    ):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty / (A * (del_T_1 - del_T_2) * log(del_T_1 - del_T_2))
        result.append(U)
        return result
    def eqn_7_14b__del_T_1(
        self, A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float, **kwargs
    ):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty / (A * U)))
        result.append(del_T_1)
        return result
    def eqn_7_14b__del_T_2(
        self, A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float, **kwargs
    ):
        # A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty / (A * U)))
        result.append(del_T_2)
        return result
    @kwasak
    def eqn_7_15(self, U=None, sum_R=None):
        return
    def eqn_7_15__U(self, sum_R: float, **kwargs):
        # 1 / U = sum_R
        result = []
        U = 1 / sum_R
        result.append(U)
        return result
    def eqn_7_15__sum_R(self, U: float, **kwargs):
        # 1 / U = sum_R
        result = []
        sum_R = 1 / U
        result.append(sum_R)
        return result
    @kwasak
    def eqn_7_16(
        self,
        D_0=None,
        D_LM=None,
        D_i=None,
        R_f_0=None,
        R_fi=None,
        U_0=None,
        h_0=None,
        h_i=None,
        k_w=None,
        x_w=None,
    ):
        return
    def eqn_7_16__D_0(
        self,
        D_LM: float,
        D_i: float,
        R_f_0: float,
        R_fi: float,
        U_0: float,
        h_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = (
            D_LM
            * D_i
            * h_i
            * k_w
            * (-R_f_0 * U_0 * h_0 - U_0 + h_0)
            / (U_0 * h_0 * (D_LM * R_fi * h_i * k_w + D_LM * k_w + D_i * h_i * x_w))
        )
        result.append(D_0)
        return result
    def eqn_7_16__D_LM(
        self,
        D_0: float,
        D_i: float,
        R_f_0: float,
        R_fi: float,
        U_0: float,
        h_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = (
            -D_0
            * D_i
            * U_0
            * h_0
            * h_i
            * x_w
            / (
                k_w
                * (
                    D_0 * R_fi * U_0 * h_0 * h_i
                    + D_0 * U_0 * h_0
                    + D_i * R_f_0 * U_0 * h_0 * h_i
                    + D_i * U_0 * h_i
                    - D_i * h_0 * h_i
                )
            )
        )
        result.append(D_LM)
        return result
    def eqn_7_16__D_i(
        self,
        D_0: float,
        D_LM: float,
        R_f_0: float,
        R_fi: float,
        U_0: float,
        h_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = (
            -D_0
            * D_LM
            * U_0
            * h_0
            * k_w
            * (R_fi * h_i + 1)
            / (
                h_i
                * (
                    D_0 * U_0 * h_0 * x_w
                    + D_LM * R_f_0 * U_0 * h_0 * k_w
                    + D_LM * U_0 * k_w
                    - D_LM * h_0 * k_w
                )
            )
        )
        result.append(D_i)
        return result
    def eqn_7_16__R_f_0(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        U_0: float,
        h_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = (
            -D_0 * R_fi / D_i
            - D_0 / (D_i * h_i)
            - D_0 * x_w / (D_LM * k_w)
            - 1 / h_0
            + 1 / U_0
        )
        result.append(R_f_0)
        return result
    def eqn_7_16__R_fi(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_f_0: float,
        U_0: float,
        h_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = (
            -1 / h_i
            - D_i * x_w / (D_LM * k_w)
            - D_i * R_f_0 / D_0
            - D_i / (D_0 * h_0)
            + D_i / (D_0 * U_0)
        )
        result.append(R_fi)
        return result
    def eqn_7_16__U_0(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_f_0: float,
        R_fi: float,
        h_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = (
            D_LM
            * D_i
            * h_0
            * h_i
            * k_w
            / (
                D_0 * D_LM * R_fi * h_0 * h_i * k_w
                + D_0 * D_LM * h_0 * k_w
                + D_0 * D_i * h_0 * h_i * x_w
                + D_LM * D_i * R_f_0 * h_0 * h_i * k_w
                + D_LM * D_i * h_i * k_w
            )
        )
        result.append(U_0)
        return result
    def eqn_7_16__h_0(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_f_0: float,
        R_fi: float,
        U_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_0 = (
            -D_LM
            * D_i
            * U_0
            * h_i
            * k_w
            / (
                D_0 * D_LM * R_fi * U_0 * h_i * k_w
                + D_0 * D_LM * U_0 * k_w
                + D_0 * D_i * U_0 * h_i * x_w
                + D_LM * D_i * R_f_0 * U_0 * h_i * k_w
                - D_LM * D_i * h_i * k_w
            )
        )
        result.append(h_0)
        return result
    def eqn_7_16__h_i(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_f_0: float,
        R_fi: float,
        U_0: float,
        h_0: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = (
            -D_0
            * D_LM
            * U_0
            * h_0
            * k_w
            / (
                D_0 * D_LM * R_fi * U_0 * h_0 * k_w
                + D_0 * D_i * U_0 * h_0 * x_w
                + D_LM * D_i * R_f_0 * U_0 * h_0 * k_w
                + D_LM * D_i * U_0 * k_w
                - D_LM * D_i * h_0 * k_w
            )
        )
        result.append(h_i)
        return result
    def eqn_7_16__k_w(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_f_0: float,
        R_fi: float,
        U_0: float,
        h_0: float,
        h_i: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = (
            -D_0
            * D_i
            * U_0
            * h_0
            * h_i
            * x_w
            / (
                D_LM
                * (
                    D_0 * R_fi * U_0 * h_0 * h_i
                    + D_0 * U_0 * h_0
                    + D_i * R_f_0 * U_0 * h_0 * h_i
                    + D_i * U_0 * h_i
                    - D_i * h_0 * h_i
                )
            )
        )
        result.append(k_w)
        return result
    def eqn_7_16__x_w(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_f_0: float,
        R_fi: float,
        U_0: float,
        h_0: float,
        h_i: float,
        k_w: float,
        **kwargs,
    ):
        # 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = (
            -D_LM * R_fi * k_w / D_i
            - D_LM * k_w / (D_i * h_i)
            - D_LM * R_f_0 * k_w / D_0
            - D_LM * k_w / (D_0 * h_0)
            + D_LM * k_w / (D_0 * U_0)
        )
        result.append(x_w)
        return result
    @kwasak
    def eqn_7_17(self, R_0=None, R_nc=None, h_c=None):
        return
    def eqn_7_17__R_0(self, R_nc: float, h_c: float, **kwargs):
        # R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1 / h_c
        result.append(R_0)
        return result
    def eqn_7_17__R_nc(self, R_0: float, h_c: float, **kwargs):
        # R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1 / h_c
        result.append(R_nc)
        return result
    def eqn_7_17__h_c(self, R_0: float, R_nc: float, **kwargs):
        # R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1 / (R_0 - R_nc)
        result.append(h_c)
        return result
    @kwasak
    def eqn_7_18(
        self,
        D_0=None,
        D_LM=None,
        D_i=None,
        R_fi=None,
        R_fo=None,
        R_nc=None,
        U_0=None,
        h_c=None,
        h_i=None,
        k_w=None,
        x_w=None,
    ):
        return
    def eqn_7_18__D_0(
        self,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = (
            D_LM
            * D_i
            * h_i
            * k_w
            * (-R_fo * U_0 * h_c - R_nc * U_0 * h_c - U_0 + h_c)
            / (U_0 * h_c * (D_LM * R_fi * h_i * k_w + D_LM * k_w + D_i * h_i * x_w))
        )
        result.append(D_0)
        return result
    def eqn_7_18__D_LM(
        self,
        D_0: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = (
            -D_0
            * D_i
            * U_0
            * h_c
            * h_i
            * x_w
            / (
                k_w
                * (
                    D_0 * R_fi * U_0 * h_c * h_i
                    + D_0 * U_0 * h_c
                    + D_i * R_fo * U_0 * h_c * h_i
                    + D_i * R_nc * U_0 * h_c * h_i
                    + D_i * U_0 * h_i
                    - D_i * h_c * h_i
                )
            )
        )
        result.append(D_LM)
        return result
    def eqn_7_18__D_i(
        self,
        D_0: float,
        D_LM: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = (
            -D_0
            * D_LM
            * U_0
            * h_c
            * k_w
            * (R_fi * h_i + 1)
            / (
                h_i
                * (
                    D_0 * U_0 * h_c * x_w
                    + D_LM * R_fo * U_0 * h_c * k_w
                    + D_LM * R_nc * U_0 * h_c * k_w
                    + D_LM * U_0 * k_w
                    - D_LM * h_c * k_w
                )
            )
        )
        result.append(D_i)
        return result
    def eqn_7_18__R_fi(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = (
            -1 / h_i
            - D_i * x_w / (D_LM * k_w)
            - D_i * R_fo / D_0
            - D_i * R_nc / D_0
            - D_i / (D_0 * h_c)
            + D_i / (D_0 * U_0)
        )
        result.append(R_fi)
        return result
    def eqn_7_18__R_fo(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fo = (
            -D_0 * R_fi / D_i
            - D_0 / (D_i * h_i)
            - D_0 * x_w / (D_LM * k_w)
            - R_nc
            - 1 / h_c
            + 1 / U_0
        )
        result.append(R_fo)
        return result
    def eqn_7_18__R_nc(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        U_0: float,
        h_c: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_nc = (
            -D_0 * R_fi / D_i
            - D_0 / (D_i * h_i)
            - D_0 * x_w / (D_LM * k_w)
            - R_fo
            - 1 / h_c
            + 1 / U_0
        )
        result.append(R_nc)
        return result
    def eqn_7_18__U_0(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        h_c: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = (
            D_LM
            * D_i
            * h_c
            * h_i
            * k_w
            / (
                D_0 * D_LM * R_fi * h_c * h_i * k_w
                + D_0 * D_LM * h_c * k_w
                + D_0 * D_i * h_c * h_i * x_w
                + D_LM * D_i * R_fo * h_c * h_i * k_w
                + D_LM * D_i * R_nc * h_c * h_i * k_w
                + D_LM * D_i * h_i * k_w
            )
        )
        result.append(U_0)
        return result
    def eqn_7_18__h_c(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_i: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_c = (
            -D_LM
            * D_i
            * U_0
            * h_i
            * k_w
            / (
                D_0 * D_LM * R_fi * U_0 * h_i * k_w
                + D_0 * D_LM * U_0 * k_w
                + D_0 * D_i * U_0 * h_i * x_w
                + D_LM * D_i * R_fo * U_0 * h_i * k_w
                + D_LM * D_i * R_nc * U_0 * h_i * k_w
                - D_LM * D_i * h_i * k_w
            )
        )
        result.append(h_c)
        return result
    def eqn_7_18__h_i(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        k_w: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = (
            -D_0
            * D_LM
            * U_0
            * h_c
            * k_w
            / (
                D_0 * D_LM * R_fi * U_0 * h_c * k_w
                + D_0 * D_i * U_0 * h_c * x_w
                + D_LM * D_i * R_fo * U_0 * h_c * k_w
                + D_LM * D_i * R_nc * U_0 * h_c * k_w
                + D_LM * D_i * U_0 * k_w
                - D_LM * D_i * h_c * k_w
            )
        )
        result.append(h_i)
        return result
    def eqn_7_18__k_w(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        h_i: float,
        x_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = (
            -D_0
            * D_i
            * U_0
            * h_c
            * h_i
            * x_w
            / (
                D_LM
                * (
                    D_0 * R_fi * U_0 * h_c * h_i
                    + D_0 * U_0 * h_c
                    + D_i * R_fo * U_0 * h_c * h_i
                    + D_i * R_nc * U_0 * h_c * h_i
                    + D_i * U_0 * h_i
                    - D_i * h_c * h_i
                )
            )
        )
        result.append(k_w)
        return result
    def eqn_7_18__x_w(
        self,
        D_0: float,
        D_LM: float,
        D_i: float,
        R_fi: float,
        R_fo: float,
        R_nc: float,
        U_0: float,
        h_c: float,
        h_i: float,
        k_w: float,
        **kwargs,
    ):
        # 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = (
            -D_LM * R_fi * k_w / D_i
            - D_LM * k_w / (D_i * h_i)
            - D_LM * R_fo * k_w / D_0
            - D_LM * R_nc * k_w / D_0
            - D_LM * k_w / (D_0 * h_c)
            + D_LM * k_w / (D_0 * U_0)
        )
        result.append(x_w)
        return result
    @kwasak
    def eqn_7_2(self, P_i_0=None, p_i=None, x_i=None):
        return
    def eqn_7_2__P_i_0(self, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i / x_i
        result.append(P_i_0)
        return result
    def eqn_7_2__p_i(self, P_i_0: float, x_i: float, **kwargs):
        # p_i = x_i * P_i_0
        result = []
        p_i = P_i_0 * x_i
        result.append(p_i)
        return result
    def eqn_7_2__x_i(self, P_i_0: float, p_i: float, **kwargs):
        # p_i = x_i * P_i_0
        result = []
        x_i = p_i / P_i_0
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_3(self, P_i_0=None, epsilon_i=None, p_i=None, x_i=None):
        return
    def eqn_7_3__P_i_0(self, epsilon_i: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i / (epsilon_i * x_i)
        result.append(P_i_0)
        return result
    def eqn_7_3__epsilon_i(self, P_i_0: float, p_i: float, x_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i / (P_i_0 * x_i)
        result.append(epsilon_i)
        return result
    def eqn_7_3__p_i(self, P_i_0: float, epsilon_i: float, x_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0 * epsilon_i * x_i
        result.append(p_i)
        return result
    def eqn_7_3__x_i(self, P_i_0: float, epsilon_i: float, p_i: float, **kwargs):
        # p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i / (P_i_0 * epsilon_i)
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_4a(self, P=None, p_c=None, p_nc=None):
        return
    def eqn_7_4a__P(self, p_c: float, p_nc: float, **kwargs):
        # p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return result
    def eqn_7_4a__p_c(self, P: float, p_nc: float, **kwargs):
        # p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return result
    def eqn_7_4a__p_nc(self, P: float, p_c: float, **kwargs):
        # p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return result
    @kwasak
    def eqn_7_4aa(self, n_i=None, n_nc=None, p_i=None, p_nc=None):
        return
    def eqn_7_4aa__n_i(self, n_nc: float, p_i: float, p_nc: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        n_i = n_nc * p_i / p_nc
        result.append(n_i)
        return result
    def eqn_7_4aa__n_nc(self, n_i: float, p_i: float, p_nc: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        n_nc = n_i * p_nc / p_i
        result.append(n_nc)
        return result
    def eqn_7_4aa__p_i(self, n_i: float, n_nc: float, p_nc: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        p_i = n_i * p_nc / n_nc
        result.append(p_i)
        return result
    def eqn_7_4aa__p_nc(self, n_i: float, n_nc: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / p_nc
        result = []
        p_nc = n_nc * p_i / n_i
        result.append(p_nc)
        return result
    @kwasak
    def eqn_7_4ab(self, P_c=None, p=None, p_i=None, p_nc=None):
        return
    def eqn_7_4ab__P_c(self, p: float, p_i: float, p_nc: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
        return result
    def eqn_7_4ab__p(self, P_c: float, p_i: float, p_nc: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result
    def eqn_7_4ab__p_i(self, P_c: float, p: float, p_nc: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = 0
        result.append(p_i)
        return result
    def eqn_7_4ab__p_nc(self, P_c: float, p: float, p_i: float, **kwargs):
        # p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = -P_c + p
        result.append(p_nc)
        return result
    @kwasak
    def eqn_7_4ac(self, P_c=None, n_i=None, n_nc=None, p=None, p_i=None):
        return
    def eqn_7_4ac__P_c(self, n_i: float, n_nc: float, p: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        P_c = p - n_nc * p_i / n_i
        result.append(P_c)
        return result
    def eqn_7_4ac__n_i(self, P_c: float, n_nc: float, p: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc * p_i / (-P_c + p)
        result.append(n_i)
        return result
    def eqn_7_4ac__n_nc(self, P_c: float, n_i: float, p: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i * (-P_c + p) / p_i
        result.append(n_nc)
        return result
    def eqn_7_4ac__p(self, P_c: float, n_i: float, n_nc: float, p_i: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        p = P_c + n_nc * p_i / n_i
        result.append(p)
        return result
    def eqn_7_4ac__p_i(self, P_c: float, n_i: float, n_nc: float, p: float, **kwargs):
        # n_i / n_nc = p_i / (p - P_c)
        result = []
        p_i = n_i * (-P_c + p) / n_nc
        result.append(p_i)
        return result
    @kwasak
    def eqn_7_5(self, N_i=None, N_nc=None, P=None, P_c=None, p_i=None):
        return
    def eqn_7_5__N_i(self, N_nc: float, P: float, P_c: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc * p_i / (P - P_c)
        result.append(N_i)
        return result
    def eqn_7_5__N_nc(self, N_i: float, P: float, P_c: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i * (P - P_c) / p_i
        result.append(N_nc)
        return result
    def eqn_7_5__P(self, N_i: float, N_nc: float, P_c: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc * p_i / N_i
        result.append(P)
        return result
    def eqn_7_5__P_c(self, N_i: float, N_nc: float, P: float, p_i: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc * p_i / N_i
        result.append(P_c)
        return result
    def eqn_7_5__p_i(self, N_i: float, N_nc: float, P: float, P_c: float, **kwargs):
        # N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i * (P - P_c) / N_nc
        result.append(p_i)
        return result
    @kwasak
    def eqn_7_6(
        self, M=None, P=None, P_i_0=None, W_air=None, W_i=None, p_c=None, x_i=None
    ):
        return
    def eqn_7_6__M(
        self,
        P: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29 * W_i * (P - p_c) / (P_i_0 * W_air * x_i)
        result.append(M)
        return result
    def eqn_7_6__P(
        self,
        M: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M * P_i_0 * W_air * x_i / (29 * W_i) + p_c
        result.append(P)
        return result
    def eqn_7_6__P_i_0(
        self, M: float, P: float, W_air: float, W_i: float, p_c: float, x_i: float, **kwargs
    ):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29 * W_i * (P - p_c) / (M * W_air * x_i)
        result.append(P_i_0)
        return result
    def eqn_7_6__W_air(
        self, M: float, P: float, P_i_0: float, W_i: float, p_c: float, x_i: float, **kwargs
    ):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29 * W_i * (P - p_c) / (M * P_i_0 * x_i)
        result.append(W_air)
        return result
    def eqn_7_6__W_i(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_air: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M * P_i_0 * W_air * x_i / (29 * (P - p_c))
        result.append(W_i)
        return result
    def eqn_7_6__p_c(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M * P_i_0 * W_air * x_i / (29 * W_i) + P
        result.append(p_c)
        return result
    def eqn_7_6__x_i(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        p_c: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29 * W_i * (P - p_c) / (M * P_i_0 * W_air)
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_7(
        self,
        M=None,
        P=None,
        P_i_0=None,
        W_air=None,
        W_i=None,
        epsilon_i=None,
        p_c=None,
        x_i=None,
    ):
        return
    def eqn_7_7__M(
        self,
        P: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        epsilon_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29 * W_i * (P - p_c) / (P_i_0 * W_air * epsilon_i * x_i)
        result.append(M)
        return result
    def eqn_7_7__P(
        self,
        M: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        epsilon_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M * P_i_0 * W_air * epsilon_i * x_i / (29 * W_i) + p_c
        result.append(P)
        return result
    def eqn_7_7__P_i_0(
        self,
        M: float,
        P: float,
        W_air: float,
        W_i: float,
        epsilon_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29 * W_i * (P - p_c) / (M * W_air * epsilon_i * x_i)
        result.append(P_i_0)
        return result
    def eqn_7_7__W_air(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_i: float,
        epsilon_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29 * W_i * (P - p_c) / (M * P_i_0 * epsilon_i * x_i)
        result.append(W_air)
        return result
    def eqn_7_7__W_i(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_air: float,
        epsilon_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M * P_i_0 * W_air * epsilon_i * x_i / (29 * (P - p_c))
        result.append(W_i)
        return result
    def eqn_7_7__epsilon_i(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        p_c: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29 * W_i * (P - p_c) / (M * P_i_0 * W_air * x_i)
        result.append(epsilon_i)
        return result
    def eqn_7_7__p_c(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        epsilon_i: float,
        x_i: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M * P_i_0 * W_air * epsilon_i * x_i / (29 * W_i) + P
        result.append(p_c)
        return result
    def eqn_7_7__x_i(
        self,
        M: float,
        P: float,
        P_i_0: float,
        W_air: float,
        W_i: float,
        epsilon_i: float,
        p_c: float,
        **kwargs,
    ):
        # W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29 * W_i * (P - p_c) / (M * P_i_0 * W_air * epsilon_i)
        result.append(x_i)
        return result
    @kwasak
    def eqn_7_8(self, L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None):
        return
    def eqn_7_8__L_c(
        self, Q_condensor_heat_duty: float, c_p: float, del_T: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty / (c_p * del_T)
        result.append(L_c)
        return result
    def eqn_7_8__Q_condensor_heat_duty(
        self, L_c: float, c_p: float, del_T: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c * c_p * del_T
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_8__c_p(
        self, L_c: float, Q_condensor_heat_duty: float, del_T: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty / (L_c * del_T)
        result.append(c_p)
        return result
    def eqn_7_8__del_T(
        self, L_c: float, Q_condensor_heat_duty: float, c_p: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty / (L_c * c_p)
        result.append(del_T)
        return result
    @kwasak
    def eqn_7_9(
        self, L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, rho=None
    ):
        return
    def eqn_7_9__L_c(
        self, Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746 * Q_condensor_heat_duty / (c_p * del_T * rho)
        result.append(L_c)
        return result
    def eqn_7_9__Q_condensor_heat_duty(
        self, L_c: float, c_p: float, del_T: float, rho: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02 * L_c * c_p * del_T * rho
        result.append(Q_condensor_heat_duty)
        return result
    def eqn_7_9__c_p(
        self, L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746 * Q_condensor_heat_duty / (L_c * del_T * rho)
        result.append(c_p)
        return result
    def eqn_7_9__del_T(
        self, L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746 * Q_condensor_heat_duty / (L_c * c_p * rho)
        result.append(del_T)
        return result
    def eqn_7_9__rho(
        self, L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float, **kwargs
    ):
        # L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746 * Q_condensor_heat_duty / (L_c * c_p * del_T)
        result.append(rho)
        return result
