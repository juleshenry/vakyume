from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton


class RotaryPistonVane:

    @staticmethod
    def eqn_11_1__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_1__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_1__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_1__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

    @staticmethod
    def eqn_11_1__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_1__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_2__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        """
        Solve the given equation for Q_0.
        
        Parameters:
            Q (float): Initially defined gas throughput value.
            Q_external_gas_throughput (float): External gas throughput value.
            SP_1 (float): Value 1 in the natural logarithm argument of the equation.
            SP_2 (float): Value 2 in the natural logarithm argument of the equation.
            S_vol_pump_speed (float): Volumetric pump speed value.
            V (float): Constant volume value.
            t (float): Time variable.
            
        Returns:
            Q_0 (float): Calculated gas throughput value.
        """
        
        # Step 1: Move -t to the left side of the equation and divide both sides by V/S_vol_pump_speed
        natural_log_arg = (SP_1 - (Q_external_gas_throughput + Q)) / (SP_2 - (Q + Q))
        
        # Step 2: Eliminate the natural logarithm by taking exponent on both sides
        exp_value = exp(t * S_vol_pump_speed / V)
        
        # Step 3: Solve for Q_0 by rearranging the equation and isolating it
        Q_0 = (SP_1 - SP_2 + Q_external_gas_throughput) / (exp_value - 1)
        
        return [ Q_0 ]


    @staticmethod
    def eqn_11_2__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        """
        Calculate the value of t based on the given equation.
        Parameters:
        Q (float): The variable in terms of which other parameters are expressed.
        Q_0 (float): A constant representing initial quantity or base value.
        Q_external_gas_throughput (float): External gas throughput quantity.
        SP_1 (float): First set point pressure.
        SP_2 (float): Second set point pressure.
        S_vol_pump_speed (float): Volume pump speed.
        V (float): A constant representing some physical property or initial value.
        
        Returns:
        float: The calculated time t based on the given equation.
        """
        numerator = SP_1 - (Q + Q_0)
        denominator = SP_2 - (Q + Q_0)
        
        # Ensure that denominator is not zero to avoid math error
        if denominator == 0:
            raise ValueError("SP_2 must be greater than SP_1 + Q_0 to avoid division by zero.")
            
        ln_ratio = log(numerator / denominator)
        
        t = -V / S_vol_pump_speed * ln_ratio
        
        return [ t ]


    @staticmethod
    def eqn_11_2__V(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        """
        Solve the equation for V using given parameters.
        
        Parameters:
        Q (float): The gas throughput rate in question units per time.
        Q_0 (float): Initial gas throughput rate at a reference state.
        Q_external_gas_throughput (float): External gas throughput rate.
        SP_1 (float): Related pressure point 1 value.
        SP_2 (float): Related pressure point 2 value.
        S_vol_pump_speed (float): Speed of the volume pump in units per time.
        t (float): Time period over which system operates, in appropriate time unit.
        
        Returns:
        float: The calculated value for V.
        """
        numerator = Q - (Q_external_gas_throughput + Q_0)
        denominator = SP_2 - (Q + Q_0)
        
        if denominator == 0:
            raise ValueError("SP_2 must be different from the sum of Q and Q_0 to avoid division by zero.")
        
        result = S_vol_pump_speed * log(numerator / denominator) * t
        return [ result ]


    @staticmethod
    def eqn_11_2__SP_2(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        numerator = SP_1 - (Q + Q_0) + exp(V / S_vol_pump_speed * t) * (Q + Q_0)
        denominator = 1 + exp(V / S_vol_pump_speed * t)
        
        SP_2 = numerator / denominator
        return [ SP_2 ]


    @staticmethod
    def eqn_11_2__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]


    @staticmethod
    def eqn_11_2__SP_1(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]


    @staticmethod
    def eqn_11_2__S_vol_pump_speed(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]


    @staticmethod
    def eqn_11_2__Q_external_gas_throughput(Q: float, Q_0: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        """
        Calculate the external gas throughput based on given parameters and an equation derived from isolating the variable.
        
        Parameters:
        - Q (float): An initial quantity or related variable to Q in terms of its derivative.
        - Q_0 (float): Constant term added to Q for specific conditions, assuming it affects throughput directly.
        - SP_1 (float): Upper static pressure influencing the logarithmic part of the equation.
        - SP_2 (float): Lower static pressure in relation to Q and its derivative or change over time.
        - S_vol_pump_speed (float): Volume pump speed affecting the rate at which Q is processed, related to how throughput changes.
        - V (float): A volume-related variable that could impact Q through its effect on processing rates or conditions.
        - t (float): Time parameter impacting the change in system behavior over time, directly influencing Q_throughput calculation.
        
        Returns:
        - float: Calculated external gas throughput (`Q_external_gas_throughput`) based on the given parameters and equation.
        """
        # Step 4 rephrased to use variables in terms of their relation to Q
        SP_1 -= Q + Q_0
        SP_2 -= V / S_vol_pump_speed * t
        
        # Calculate the term inside the exponent using adjusted parameters
        exponent_term = exp(t / (V / S_vol_pump_speed))
        
        # Isolate Q_external_gas_throughput in the equation
        Q_external_gas_throughput = SP_1 - SP_2 * exponent_term - Q_0
        
        return [ Q_external_gas_throughput ]


    @staticmethod
    def eqn_11_3__t_c(F_s: float, t: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t_c = t/F_s
        result.append(t_c)
        return result

    @staticmethod
    def eqn_11_3__t(F_s: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t = F_s*t_c
        result.append(t)
        return result

    @staticmethod
    def eqn_11_3__F_s(t: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        F_s = t/t_c
        result.append(F_s)
        return result

    @staticmethod
    def eqn_11_4__p_s(p_g: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_s = p_g + p_v
        result.append(p_s)
        return result

    @staticmethod
    def eqn_11_4__p_v(p_g: float, p_s: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_v = 0
        result.append(p_v)
        p_v = -p_g + p_s
        result.append(p_v)
        return result

    @staticmethod
    def eqn_11_4__p_g(p_s: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_g = p_s - p_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_5__P_0_v(P_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D*p_v_max/(p_g + p_v_max)
        result.append(P_0_v)
        return result

    @staticmethod
    def eqn_11_5__p_v_max(P_0_v: float, P_D: float, p_g: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v*p_g/(P_0_v - P_D)
        result.append(p_v_max)
        return result

    @staticmethod
    def eqn_11_5__P_D(P_0_v: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v*(p_g + p_v_max)/p_v_max
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_5__p_g(P_0_v: float, P_D: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max*(-P_0_v + P_D)/P_0_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_6__S_D(P_0_V: float, P_D: float, P_v_0: float, S_B: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = P_D*S_B*(P_0_V - p_b)/(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)
        result.append(S_D)
        return result

    @staticmethod
    def eqn_11_6__p_b(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (P_0_V*P_D*S_B - P_D*S_D*p_v_max + P_v_0*S_D*p_g + P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(p_b)
        return result

    @staticmethod
    def eqn_11_6__S_B(P_0_V: float, P_D: float, P_v_0: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = S_D*(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)/(P_D*(P_0_V - p_b))
        result.append(S_B)
        return result

    @staticmethod
    def eqn_11_6__P_0_V(P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_g - P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(P_0_V)
        return result

    @staticmethod
    def eqn_11_6__P_D(P_0_V: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0*S_D*(p_g + p_v_max)/(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_6__P_v_0(P_0_V: float, P_D: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = P_D*(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)/(S_D*(p_g + p_v_max))
        result.append(P_v_0)
        return result

    @staticmethod
    def eqn_11_6__p_v_max(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V*P_D*S_B - P_D*S_B*p_b + P_v_0*S_D*p_g)/(S_D*(P_D - P_v_0))
        result.append(p_v_max)
        return result

    @staticmethod
    def eqn_11_6__p_g(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (-P_0_V*P_D*S_B + P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_v_max)/(P_v_0*S_D)
        result.append(p_g)
        return result
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton


class RotaryPistonVane:

    @staticmethod
    def eqn_11_1__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_1__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_1__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_1__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_1__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_1__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

    @staticmethod
    def eqn_11_2__SP_2(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton


class RotaryPistonVane:

    @staticmethod
    def eqn_11_1__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_1__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_1__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_1__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

    @staticmethod
    def eqn_11_1__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_1__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_2__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton


class RotaryPistonVane:

    @staticmethod
    def eqn_11_1__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

    @staticmethod
    def eqn_11_1__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_1__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_1__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_1__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_1__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_2__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
