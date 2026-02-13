from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class Precondensors:
    @kwasak_static
    def eqn_7_14b(
        A: float = None,
        Q_condensor_heat_duty: float = None,
        U: float = None,
        del_T_1: float = None,
        del_T_2: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_7_14b__A(
        Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float
    ):
        result = []
        try:
            # LMTD is only physical for positive temperature differences
            t1, t2 = abs(del_T_1), abs(del_T_2)
            if t1 < 1e-12 or t2 < 1e-12: return []
            if abs(t1 - t2) < 1e-7:
                L = t1
            else:
                L = (t1 - t2) / log(t1 / t2)
            A = Q_condensor_heat_duty / (U * L)
            result.append(A)
        except:
            pass
        return result

    @staticmethod
    def eqn_7_14b__Q_condensor_heat_duty(
        A: float, U: float, del_T_1: float, del_T_2: float
    ):
        result = []
        try:
            t1, t2 = abs(del_T_1), abs(del_T_2)
            if t1 < 1e-12 or t2 < 1e-12: return []
            if abs(t1 - t2) < 1e-7:
                L = t1
            else:
                L = (t1 - t2) / log(t1 / t2)
            Q = U * A * L
            result.append(Q)
        except:
            pass
        return result

    @staticmethod
    def eqn_7_14b__U(
        A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float
    ):
        result = []
        try:
            t1, t2 = abs(del_T_1), abs(del_T_2)
            if t1 < 1e-12 or t2 < 1e-12: return []
            if abs(t1 - t2) < 1e-7:
                L = t1
            else:
                L = (t1 - t2) / log(t1 / t2)
            U = Q_condensor_heat_duty / (A * L)
            result.append(U)
        except:
            pass
        return result

    @staticmethod
    def eqn_7_14b__del_T_1(
        A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float
    ):
        try:
            from scipy.optimize import brentq
            K = Q_condensor_heat_duty / (U * A)
            t2 = abs(del_T_2)
            if t2 < 1e-12: return []
            
            def func(t1):
                if abs(t1 - t2) < 1e-9: return t1 - K
                lmtd = (t1 - t2) / log(t1 / t2)
                return lmtd - K

            if abs(K - t2) < 1e-7: return [t2]
            
            if K > t2:
                # Root is in (K, inf).
                upper = max(K * 2, t2 * 2, 100.0)
                while func(upper) < 0 and upper < 1e15:
                    upper *= 10
                return [float(brentq(func, t2 + 1e-12, upper))]
            else:
                # Root is in (0, t2).
                return [float(brentq(func, 1e-15, t2 - 1e-12))]
        except:
            return []

    @staticmethod
    def eqn_7_14b__del_T_2(
        A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float
    ):
        try:
            from scipy.optimize import brentq
            K = Q_condensor_heat_duty / (U * A)
            t1 = abs(del_T_1)
            if t1 < 1e-12: return []
            
            def func(t2):
                if abs(t1 - t2) < 1e-9: return t2 - K
                lmtd = (t1 - t2) / log(t1 / t2)
                return lmtd - K

            if abs(K - t1) < 1e-7: return [t1]
            
            if K > t1:
                # Root is in (K, inf).
                upper = max(K * 2, t1 * 2, 100.0)
                while func(upper) < 0 and upper < 1e15:
                    upper *= 10
                return [float(brentq(func, t1 + 1e-12, upper))]
            else:
                # Root is in (0, t1).
                return [float(brentq(func, 1e-15, t1 - 1e-12))]
        except:
            return []


