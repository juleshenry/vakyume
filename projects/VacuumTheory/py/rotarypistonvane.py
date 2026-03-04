from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
from vakyume.kwasak import kwasak
from vakyume.config import UnsolvedException
import numpy as np


class RotaryPistonVane:
    @kwasak
    def eqn_11_1(
        self,
        PS=None,
        Q_0=None,
        Q_external_gas_throughput=None,
        V=None,
        dP=None,
        dT=None,
    ):
        return
    def eqn_11_1__PS(
        self,
        Q_0: float,
        Q_external_gas_throughput: float,
        V: float,
        dP: float,
        dT: float,
        **kwargs,
    ):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V * dP / dT
        result.append(PS)
        return result
    def eqn_11_1__Q_0(
        self,
        PS: float,
        Q_external_gas_throughput: float,
        V: float,
        dP: float,
        dT: float,
        **kwargs,
    ):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V * dP / dT
        result.append(Q_0)
        return result
    def eqn_11_1__Q_external_gas_throughput(
        self, PS: float, Q_0: float, V: float, dP: float, dT: float, **kwargs
    ):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V * dP / dT
        result.append(Q_external_gas_throughput)
        return result
    def eqn_11_1__V(
        self,
        PS: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        dP: float,
        dT: float,
        **kwargs,
    ):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT * (-PS + Q_0 + Q_external_gas_throughput) / dP
        result.append(V)
        return result
    def eqn_11_1__dP(
        self,
        PS: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        V: float,
        dT: float,
        **kwargs,
    ):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT * (-PS + Q_0 + Q_external_gas_throughput) / V
        result.append(dP)
        return result
    def eqn_11_1__dT(
        self,
        PS: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        V: float,
        dP: float,
        **kwargs,
    ):
        # PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V * dP / (-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result
    @kwasak
    def eqn_11_2(
        self,
        Q=None,
        Q_0=None,
        Q_external_gas_throughput=None,
        SP_1=None,
        SP_2=None,
        S_vol_pump_speed=None,
        V=None,
        t=None,
    ):
        return
    def eqn_11_2__Q(
        self,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
        **kwargs,
    ):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q = (
            Q_0
            + Q_external_gas_throughput
            - SP_1
            + (-Q_0 + SP_2) * exp(S_vol_pump_speed * t / V)
        ) * exp(-S_vol_pump_speed * t / V)
        result.append(Q)
        return result
    def eqn_11_2__Q_0(
        self,
        Q: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
        **kwargs,
    ):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        Q_0 = -(
            Q * exp(S_vol_pump_speed * t / V)
            - Q_external_gas_throughput
            + SP_1
            - SP_2 * exp(S_vol_pump_speed * t / V)
        ) / (exp(S_vol_pump_speed * t / V) - 1)
        result.append(Q_0)
        return result
    def eqn_11_2__Q_external_gas_throughput(
        self,
        Q: float,
        Q_0: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
        **kwargs,
    ):
        # Solve for Q_external_gas_throughput by rearranging the equation
        try:
            result = []
            term1 = (Q - Q_0) * exp(S_vol_pump_speed * t / V) + SP_1
            term2 = (-SP_2 - Q_0) * exp(-S_vol_pump_speed * t / V)
            # Isolate the exponential terms and solve for e^(S_vol_pump_speed * t / V)
            exponent = (term1 + SP_2 - Q_0) / term2
            inv_exponent = 1 / exponent
            exp_solution = exp(inv_exponent)
            # Isolate the remaining terms and solve for e^(-S_vol_pump_speed * t / V)
            base_term = (Q - Q_0) + SP_2
            term3 = S_vol_pump_speed * t / V
            inv_base_term = 1 / base_term
            # Solve for the unknown variable by isolating it on one side of the equation
            result.append(SP_2 - exp_solution * (Q - Q_0) + term3 * inv_exponent)
        except ZeroDivisionError:
            print("Error: Division by zero encountered.")
        return [result[0]] if result else None
    def eqn_11_2__SP_1(
        self,
        Q: complex,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
        **kwargs,
    ):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        if isinstance(Q, complex):
            Q_extended = complex(
                Q.real + Q_0
            )  # Combine real and imaginary parts of Q if it's a complex number
        else:
            Q_extended = Q + Q_0

        numerator = V * t - (SP_2 - SP_1) * log(
            (-t / S_vol_pump_speed) + exp(-t / S_vol_pump_speed)
        )  # Rearrange the equation to isolate SP_1 on one side

        denominator = -(V / S_vol_pump_speed) - t / (exp(t / S_vol_pump_speed))

        try:
            SP_1 = numerator / denominator  # Calculate SP_1 using the rearranged equation
            result.append(SP_1)
        except ZeroDivisionError:
            print("Error: Division by zero encountered")
            return []

        return [result]
    def eqn_11_2__SP_2(
        self,
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
        **kwargs,
    ):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        result = []
        SP_2 = (
            -Q_0
            - Q_external_gas_throughput
            + SP_1
            + (Q + Q_0) * exp(S_vol_pump_speed * t / V)
        ) * exp(-S_vol_pump_speed * t / V)
        result.append(SP_2)
        return result
    def eqn_11_2__S_vol_pump_speed(
        self,
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        V: float,
        t: float,
        **kwargs,
    ):
        # Solve for S_vol_pump_speed by rearranging the equation

        s = symbols("s")
        Q_eq = (Q - Q_0) + (-Q_0 + SP_2) * exp(s * t / V) * exp(-s * t / V)
        S_vol_pump_speed = solve(Eq(Q_eq, 0), s)[0]

        return [S_vol_pump_speed.evalf()]
    def eqn_11_2__V(
        self,
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        t: float,
        **kwargs,
    ):
        # Solve for V by rearranging the equation
        result = []
        try:
            A = (SP_1 - Q_external_gas_throughput + Q_0) / SP_2 - Q - Q_0
            B = S_vol_pump_speed * t / V
            C = 1 / t
            D = exp(B)
            E = log(-C)
            F = A / (E - 1)
            if not isinstance(F, complex):
                raise ValueError("Invalid result from the calculation")

            # Check for division by zero or negative argument in exponential function which would lead to a math domain error.
            if B == 0 or E <= 0:
                return [None]
            V = -t / (S_vol_pump_speed * log(F))
        except ValueError as e:
            print(e)
            return [None]
        result.append(V)
        return [result]
    def eqn_11_2__t(
        self,
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        **kwargs,
    ):
        # t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        return [
            V
            / S_vol_pump_speed
            * ln((SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0)))
        ]
    @kwasak
    def eqn_11_3(self, F_s=None, t=None, t_c=None):
        return
    def eqn_11_3__F_s(self, t: float, t_c: float, **kwargs):
        # t = t_c * F_s
        result = []
        F_s = t / t_c
        result.append(F_s)
        return result
    def eqn_11_3__t(self, F_s: float, t_c: float, **kwargs):
        # t = t_c * F_s
        result = []
        t = F_s * t_c
        result.append(t)
        return result
    def eqn_11_3__t_c(self, F_s: float, t: float, **kwargs):
        # t = t_c * F_s
        result = []
        t_c = t / F_s
        result.append(t_c)
        return result
    @kwasak
    def eqn_11_4(self, p_g=None, p_s=None, p_v=None):
        return
    def eqn_11_4__p_g(self, p_s: float, p_v: float, **kwargs):
        # Solve for p_g by rearranging the equation p_v / (p_v + p_g) = 1/p_s
        def _res(p_g_val):
            try:
                # Force complex evaluation to handle negative bases in fractional powers
                target_var_complex = complex(p_g_val, 0)
                val = 1 / p_s - p_v / (p_v + p_g_val)
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
            raise UnsolvedException("No sign change found for p_g in expanded range")
        p_g = brentq(_res, lo, hi)
        return [p_g]
    def eqn_11_4__p_s(self, p_g: float, p_normalize=False, **kwargs):
        # p_v / (p_v + p_g) = p_v / p_s
        # Solve for p_s by rearranging the equation p_v / (p_v + p_g) = p_v / p_s
        if not p_normalize and "p_v" in kwargs:
            raise ValueError(
                "Keyword argument 'p_v' must be provided when normalization is False"
            )

        # Extract the value of p_v from keyword arguments or default to  end['p_s'] = float('inf') if not provided and p_normalize is True, else return [0.0]
        try:
            p_v = kwargs.get("p_v", None)
            if p_v is None:
                raise ValueError(
                    "Keyword argument 'p_v' must be provided when normalization is False"
                )

            # Rearrange the equation for p_s: (p_g / ((1 - p_v))) * p_v = p_v^2 => p_s = p_g / ((1 + p_v)) if 'p_v' is not zero, else return [0.0]
            p_s = float("inf") if p_normalize and p_v == 0 else (p_g / (1 - p_v))
        except ZeroDivisionError:
            raise ValueError("p_v must be non-zero when normalization is True")

        return [p_s]
    def eqn_11_4__p_v(self, p_g: float, p_s: float, **kwargs):
        # p_v / (p_v + p_g) = p_v / p_s
        # p_v appears 3 times — use numerical solver
        from scipy.optimize import brentq
        import numpy as np

        def _res(p_v_val):
            try:
                # Force complex evaluation to handle negative bases in fractional powers
                target_var_complex = complex(p_v_val, 0)
                val = target_var_complex / p_s - p_v / (p_v + p_g)
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
            raise UnsolvedException("No sign change found for p_v in expanded range")
        p_v = brentq(_res, lo, hi)
        return [p_v]
    @kwasak
    def eqn_11_5(self, P_0_v=None, P_D=None, p_g=None, p_v_max=None):
        return
    def eqn_11_5__P_0_v(self, P_D: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D * p_v_max / (p_g + p_v_max)
        result.append(P_0_v)
        return result
    def eqn_11_5__P_D(self, P_0_v: float, p_g: float, p_v_max: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v * (p_g + p_v_max) / p_v_max
        result.append(P_D)
        return result
    def eqn_11_5__p_g(self, P_0_v: float, P_D: float, p_v_max: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max * (-P_0_v + P_D) / P_0_v
        result.append(p_g)
        return result
    def eqn_11_5__p_v_max(self, P_0_v: float, P_D: float, p_g: float, **kwargs):
        # p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v * p_g / (P_0_v - P_D)
        result.append(p_v_max)
        return result
    @kwasak
    def eqn_11_6(
        self,
        P_0_V=None,
        P_D=None,
        P_v_0=None,
        S_B=None,
        S_D=None,
        p_b=None,
        p_g=None,
        p_v_max=None,
    ):
        return
    def eqn_11_6__P_0_V(
        self,
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (
            P_D * S_B * p_b
            + P_D * S_D * p_v_max
            - P_v_0 * S_D * p_g
            - P_v_0 * S_D * p_v_max
        ) / (P_D * S_B)
        result.append(P_0_V)
        return result
    def eqn_11_6__P_D(
        self,
        P_0_V: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0 * S_D * (p_g + p_v_max) / (-P_0_V * S_B + S_B * p_b + S_D * p_v_max)
        result.append(P_D)
        return result
    def eqn_11_6__P_v_0(
        self,
        P_0_V: float,
        P_D: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = P_D * (-P_0_V * S_B + S_B * p_b + S_D * p_v_max) / (S_D * (p_g + p_v_max))
        result.append(P_v_0)
        return result
    def eqn_11_6__S_B(
        self,
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = S_D * (P_D * p_v_max - P_v_0 * p_g - P_v_0 * p_v_max) / (P_D * (P_0_V - p_b))
        result.append(S_B)
        return result
    def eqn_11_6__S_D(
        self,
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = P_D * S_B * (P_0_V - p_b) / (P_D * p_v_max - P_v_0 * p_g - P_v_0 * p_v_max)
        result.append(S_D)
        return result
    def eqn_11_6__p_b(
        self,
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_g: float,
        p_v_max: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (
            P_0_V * P_D * S_B
            - P_D * S_D * p_v_max
            + P_v_0 * S_D * p_g
            + P_v_0 * S_D * p_v_max
        ) / (P_D * S_B)
        result.append(p_b)
        return result
    def eqn_11_6__p_g(
        self,
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_v_max: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (
            -P_0_V * P_D * S_B
            + P_D * S_B * p_b
            + P_D * S_D * p_v_max
            - P_v_0 * S_D * p_v_max
        ) / (P_v_0 * S_D)
        result.append(p_g)
        return result
    def eqn_11_6__p_v_max(
        self,
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
        **kwargs,
    ):
        # p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V * P_D * S_B - P_D * S_B * p_b + P_v_0 * S_D * p_g) / (
            S_D * (P_D - P_v_0)
        )
        result.append(p_v_max)
        return result
