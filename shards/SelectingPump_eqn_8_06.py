from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SelectingPump:
    @kwasak_static
    def eqn_8_06(
        M: float = None,
        P_1: float = None,
        P_2: float = None,
        R: float = None,
        T: float = None,
        adiabatic_hp: float = None,
        k: float = None,
        w: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_8_06__M(
        P_1: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        try:
            M = (
                R
                * T
                * k
                * w
                * (abs(P_2 / P_1) ** ((k - 1) / k) - 1)
                / (1980000 * adiabatic_hp * (k - 1))
            )
            result.append(M)
        except:
            pass
        return result

    @staticmethod
    def eqn_8_06__P_1(
        M: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        try:
            base = 1980000 * M * adiabatic_hp / (R * T * w) - 1980000 * M * adiabatic_hp / (R * T * k * w) + 1
            P_1 = P_2 / (abs(base) ** (k / (k - 1)))
            result.append(P_1)
        except:
            pass
        return result

    @staticmethod
    def eqn_8_06__P_2(
        M: float,
        P_1: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        try:
            base = 1980000 * M * adiabatic_hp / (R * T * w) - 1980000 * M * adiabatic_hp / (R * T * k * w) + 1
            P_2 = P_1 * (abs(base) ** (k / (k - 1)))
            result.append(P_2)
        except:
            pass
        return result

    @staticmethod
    def eqn_8_06__R(
        M: float,
        P_1: float,
        P_2: float,
        T: float,
        adiabatic_hp: float,
        k: float,
        w: float,
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        try:
            R = (
                1980000
                * M
                * adiabatic_hp
                * (k - 1)
                / (T * k * w * (abs(P_2 / P_1) ** ((k - 1) / k) - 1))
            )
            result.append(R)
        except:
            pass
        return result

    @staticmethod
    def eqn_8_06__T(
        M: float,
        P_1: float,
        P_2: float,
        R: float,
        adiabatic_hp: float,
        k: float,
        w: float,
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        try:
            T = (
                1980000
                * M
                * adiabatic_hp
                * (k - 1)
                / (R * k * w * (abs(P_2 / P_1) ** ((k - 1) / k) - 1))
            )
            result.append(T)
        except:
            pass
        return result

    @staticmethod
    def eqn_8_06__adiabatic_hp(
        M: float, P_1: float, P_2: float, R: float, T: float, k: float, w: float
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        try:
            adiabatic_hp = (
                R * T * k * w * (abs(P_2 / P_1) ** ((k - 1) / k) - 1) / (1980000 * M * (k - 1))
            )
            result.append(adiabatic_hp)
        except:
            pass
        return result

    @staticmethod
    def eqn_8_06__k(
        M: float,
        P_1: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        w: float,
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        
        try:
            from scipy.optimize import fsolve
            ratio = abs(P_2 / P_1)
            if abs(ratio - 1.0) < 1e-9:
                return []
            
            # Constant part: (w*R*T) / (M*1980000)
            C_const = (w * R * T) / (M * 1980000.0)
            if abs(C_const) < 1e-18: return []
            
            # We need to solve: adiabatic_hp = C_const * (k/(k-1)) * (ratio**((k-1)/k) - 1)
            def func(k_val):
                if abs(k_val) < 1e-7 or abs(k_val - 1.0) < 1e-7:
                    return 1e12 # Penalize
                x = (k_val - 1.0) / k_val
                # Use a more stable form for (ratio**x - 1)/x if needed, but ratio**x - 1 should be fine
                return C_const * (1.0/x) * (ratio ** x - 1.0) - adiabatic_hp
            
            # Try a range of guesses covering typical adiabatic indices and some outliers
            best_sol = None
            min_err = 1e9
            
            for guess in [1.4, 1.3, 1.6, 1.1, 2.0, 0.7, 5.0, -0.5]:
                sol, info, ier, msg = fsolve(func, guess, full_output=True)
                if ier == 1:
                    err = abs(func(sol[0]))
                    if err < min_err:
                        min_err = err
                        best_sol = float(sol[0])
            
            if best_sol is not None and min_err < 1e-4:
                return [best_sol]
            
            # If fsolve failed, try a very broad search with brentq if we can find a bracket
            from scipy.optimize import brentq
            def func_root(k_v):
                return func(k_v)
            
            # Test some brackets for k
            for a, b in [(1.0001, 10.0), (0.001, 0.999), (-10.0, -0.001)]:
                try:
                    if func_root(a) * func_root(b) < 0:
                        return [float(brentq(func_root, a, b))]
                except:
                    continue
                    
            return []
        except:
            return []

    @staticmethod
    def eqn_8_06__w(
        M: float,
        P_1: float,
        P_2: float,
        R: float,
        T: float,
        adiabatic_hp: float,
        k: float,
    ):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        try:
            w = (
                1980000
                * M
                * adiabatic_hp
                * (k - 1)
                / (R * T * k * (abs(P_2 / P_1) ** ((k - 1) / k) - 1))
            )
            result.append(w)
        except:
            pass
        return result


