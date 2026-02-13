from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class RotaryPistonVane:
    @kwasak_static
    def eqn_11_02(
        Q: float = None,
        Q_0: float = None,
        Q_external_gas_throughput: float = None,
        SP_1: float = None,
        SP_2: float = None,
        S_vol_pump_speed: float = None,
        V: float = None,
        t: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_11_02__Q(
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
    ):
        """
        Solve equation 11-02 for Q given all other parameters.

        Original equation:
        t = V / S_vol_pump_speed * ln((SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))

        Parameters:
        Q_0: base gas throughput
        Q_external_gas_throughput: external gas throughput
        SP_1: initial pressure
        SP_2: final pressure
        S_vol_pump_speed: volumetric pump speed
        V: volume
        t: time

        Returns:
        Q: gas throughput
        """
        # Rearrange the equation to solve for Q
        # exp(t * S_vol_pump_speed / V) = (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0))
        # (SP_2 - (Q + Q_0)) * exp(t * S_vol_pump_speed / V) = SP_1 - (Q_external_gas_throughput + Q_0)
        # SP_2 * exp(t * S_vol_pump_speed / V) - (Q + Q_0) * exp(t * S_vol_pump_speed / V) = SP_1 - Q_external_gas_throughput - Q_0
        # -Q * exp(t * S_vol_pump_speed / V) = SP_1 - Q_external_gas_throughput - Q_0 - SP_2 * exp(t * S_vol_pump_speed / V) + Q_0 * exp(t * S_vol_pump_speed / V)
        # Q = (SP_2 * exp(t * S_vol_pump_speed / V) - SP_1 + Q_external_gas_throughput + Q_0 - Q_0 * exp(t * S_vol_pump_speed / V)) / exp(t * S_vol_pump_speed / V)

        exp_term = exp(t * S_vol_pump_speed / V)

        Q = (
            SP_2 * exp_term - SP_1 + Q_external_gas_throughput + Q_0 - Q_0 * exp_term
        ) / exp_term

        # Simplify the equation
        Q = SP_2 - (SP_1 - Q_external_gas_throughput - Q_0) / exp_term - Q_0

        return [Q]

    @staticmethod
    def eqn_11_02__Q_0(
        Q: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
    ):
        """
        Solve equation 11-02 for Q_0 given all other parameters.

        Original equation:
        t = V / S_vol_pump_speed * ln((SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))

        Parameters:
        Q: gas throughput
        Q_external_gas_throughput: external gas throughput
        SP_1: initial pressure
        SP_2: final pressure
        S_vol_pump_speed: volumetric pump speed
        V: volume
        t: time

        Returns:
        Q_0: base gas throughput
        """
        exp_term = exp(t * S_vol_pump_speed / V)

        # Rearrange to solve for Q_0
        # exp_term = (SP_1 - Q_external_gas_throughput - Q_0) / (SP_2 - Q - Q_0)
        # exp_term * (SP_2 - Q - Q_0) = SP_1 - Q_external_gas_throughput - Q_0
        # exp_term * SP_2 - exp_term * Q - exp_term * Q_0 = SP_1 - Q_external_gas_throughput - Q_0
        # Q_0 - exp_term * Q_0 = SP_1 - Q_external_gas_throughput - exp_term * SP_2 + exp_term * Q
        # Q_0 * (1 - exp_term) = SP_1 - Q_external_gas_throughput - exp_term * SP_2 + exp_term * Q

        denom = 1 - exp_term
        if abs(denom) < 1e-9:
            return []

        Q_0 = (SP_1 - Q_external_gas_throughput - exp_term * SP_2 + exp_term * Q) / denom

        return [Q_0]

    @staticmethod
    def eqn_11_02__Q_external_gas_throughput(
        Q: float,
        Q_0: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
    ):
        """
        Solve equation 11-02 for Q_external_gas_throughput given all other parameters.

        Original equation:
        t = V / S_vol_pump_speed * ln((SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))

        Parameters:
        Q: gas throughput
        Q_0: base gas throughput
        SP_1: initial pressure
        SP_2: final pressure
        S_vol_pump_speed: volumetric pump speed
        V: volume
        t: time

        Returns:
        Q_external_gas_throughput: external gas throughput
        """
        exp_term = exp(t * S_vol_pump_speed / V)

        # Rearrange to solve for Q_external_gas_throughput
        # exp_term = (SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0))
        # exp_term * (SP_2 - Q - Q_0) = SP_1 - Q_external_gas_throughput - Q_0
        # Q_external_gas_throughput = SP_1 - Q_0 - exp_term * (SP_2 - Q - Q_0)

        Q_external_gas_throughput = SP_1 - Q_0 - exp_term * (SP_2 - Q - Q_0)

        return [Q_external_gas_throughput]

    @staticmethod
    def eqn_11_02__SP_1(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
    ):
        """
        Solve equation 11-02 for SP_1 given all other parameters.

        Original equation:
        t = V / S_vol_pump_speed * ln((SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))

        Parameters:
        Q: gas throughput
        Q_0: base gas throughput
        Q_external_gas_throughput: external gas throughput
        SP_2: final pressure
        S_vol_pump_speed: volumetric pump speed
        V: volume
        t: time

        Returns:
        SP_1: initial pressure
        """
        exp_term = exp(t * S_vol_pump_speed / V)

        # Rearrange to solve for SP_1
        # exp_term = (SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0))
        # SP_1 - (Q_external_gas_throughput + Q_0) = exp_term * (SP_2 - (Q + Q_0))
        # SP_1 = exp_term * (SP_2 - Q - Q_0) + Q_external_gas_throughput + Q_0

        SP_1 = exp_term * (SP_2 - Q - Q_0) + Q_external_gas_throughput + Q_0

        return [SP_1]

    @staticmethod
    def eqn_11_02__SP_2(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        S_vol_pump_speed: float,
        V: float,
        t: float,
    ):
        """
        Solve equation 11-02 for SP_2 given all other parameters.

        Original equation:
        t = V / S_vol_pump_speed * ln((SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))

        Parameters:
        Q: gas throughput
        Q_0: base gas throughput
        Q_external_gas_throughput: external gas throughput
        SP_1: initial pressure
        S_vol_pump_speed: volumetric pump speed
        V: volume
        t: time

        Returns:
        SP_2: final pressure
        """
        exp_term = exp(t * S_vol_pump_speed / V)

        # Rearrange to solve for SP_2
        # exp_term = (SP_1 - (Q_external_gas_throughput + Q_0)) / (SP_2 - (Q + Q_0))
        # exp_term * (SP_2 - Q - Q_0) = SP_1 - Q_external_gas_throughput - Q_0
        # SP_2 - Q - Q_0 = (SP_1 - Q_external_gas_throughput - Q_0) / exp_term
        # SP_2 = (SP_1 - Q_external_gas_throughput - Q_0) / exp_term + Q + Q_0

        SP_2 = (SP_1 - Q_external_gas_throughput - Q_0) / exp_term + Q + Q_0

        return [SP_2]

    @staticmethod
    def eqn_11_02__S_vol_pump_speed(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        V: float,
        t: float,
    ):
        # Rearrange to solve for S_vol_pump_speed
        # t = V / S_vol_pump_speed * ln(numerator / denominator)
        # S_vol_pump_speed = V * ln(numerator / denominator) / t

        try:
            numerator = SP_1 - (Q_external_gas_throughput + Q_0)
            denominator = SP_2 - (Q + Q_0)
            S_vol_pump_speed = V * log(abs(numerator / denominator)) / t
            return [S_vol_pump_speed]
        except:
            return []

    @staticmethod
    def eqn_11_02__V(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        t: float,
    ):
        # Rearrange to solve for V
        # t = V / S_vol_pump_speed * ln(numerator / denominator)
        # V = t * S_vol_pump_speed / ln(numerator / denominator)

        try:
            numerator = SP_1 - (Q_external_gas_throughput + Q_0)
            denominator = SP_2 - (Q + Q_0)
            log_term = log(abs(numerator / denominator))
            V = t * S_vol_pump_speed / log_term
            return [V]
        except:
            return []

    @staticmethod
    def eqn_11_02__t(
        Q: float,
        Q_0: float,
        Q_external_gas_throughput: float,
        SP_1: float,
        SP_2: float,
        S_vol_pump_speed: float,
        V: float,
    ):
        # Direct calculation using the original formula
        try:
            numerator = SP_1 - (Q_external_gas_throughput + Q_0)
            denominator = SP_2 - (Q + Q_0)
            t = V / S_vol_pump_speed * log(abs(numerator / denominator))
            return [t]
        except:
            return []


