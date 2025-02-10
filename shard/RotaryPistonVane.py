from ..kwasak import kwasak_staticclass RotaryPistonVane:

    @kwasak_static
    def eqn_11_01(PS: float = None, Q_0: float = None, Q_external_gas_throughput: float = None, V: float = None, dP: float = None, dT: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_01__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_01__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_01__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_01__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_01__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_01__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

    @kwasak_static
    def eqn_11_02(Q: float = None, Q_0: float = None, Q_external_gas_throughput: float = None, SP_1: float = None, SP_2: float = None, S_vol_pump_speed: float = None, V: float = None, t: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_02__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_11_02__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_11_02__Q_external_gas_throughput(Q: float, Q_0: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_11_02__SP_1(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_11_02__SP_2(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_11_02__S_vol_pump_speed(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_11_02__V(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_11_02__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [Sympy Failover]
        pass # Ollama offline

    @kwasak_static
    def eqn_11_03(F_s: float = None, t: float = None, t_c: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_03__F_s(t: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        F_s = t/t_c
        result.append(F_s)
        return result

    @staticmethod
    def eqn_11_03__t(F_s: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t = F_s*t_c
        result.append(t)
        return result

    @staticmethod
    def eqn_11_03__t_c(F_s: float, t: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t_c = t/F_s
        result.append(t_c)
        return result

    @kwasak_static
    def eqn_11_04(p_g: float = None, p_s: float = None, p_v: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_04__p_g(p_s: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_g = p_s - p_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_04__p_s(p_g: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_s = p_g + p_v
        result.append(p_s)
        return result

    @staticmethod
    def eqn_11_04__p_v(p_g: float, p_s: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_v = 0
        result.append(p_v)
        p_v = -p_g + p_s
        result.append(p_v)
        return result

    @kwasak_static
    def eqn_11_05(P_0_v: float = None, P_D: float = None, p_g: float = None, p_v_max: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_05__P_0_v(P_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D*p_v_max/(p_g + p_v_max)
        result.append(P_0_v)
        return result

    @staticmethod
    def eqn_11_05__P_D(P_0_v: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v*(p_g + p_v_max)/p_v_max
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_05__p_g(P_0_v: float, P_D: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max*(-P_0_v + P_D)/P_0_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_05__p_v_max(P_0_v: float, P_D: float, p_g: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v*p_g/(P_0_v - P_D)
        result.append(p_v_max)
        return result

    @kwasak_static
    def eqn_11_06(P_0_V: float = None, P_D: float = None, P_v_0: float = None, S_B: float = None, S_D: float = None, p_b: float = None, p_g: float = None, p_v_max: float = None,**kwargs):
        return


    @staticmethod
    def eqn_11_06__P_0_V(P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_g - P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(P_0_V)
        return result

    @staticmethod
    def eqn_11_06__P_D(P_0_V: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0*S_D*(p_g + p_v_max)/(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_06__P_v_0(P_0_V: float, P_D: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = P_D*(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)/(S_D*(p_g + p_v_max))
        result.append(P_v_0)
        return result

    @staticmethod
    def eqn_11_06__S_B(P_0_V: float, P_D: float, P_v_0: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = S_D*(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)/(P_D*(P_0_V - p_b))
        result.append(S_B)
        return result

    @staticmethod
    def eqn_11_06__S_D(P_0_V: float, P_D: float, P_v_0: float, S_B: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = P_D*S_B*(P_0_V - p_b)/(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)
        result.append(S_D)
        return result

    @staticmethod
    def eqn_11_06__p_b(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (P_0_V*P_D*S_B - P_D*S_D*p_v_max + P_v_0*S_D*p_g + P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(p_b)
        return result

    @staticmethod
    def eqn_11_06__p_g(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (-P_0_V*P_D*S_B + P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_v_max)/(P_v_0*S_D)
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_06__p_v_max(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V*P_D*S_B - P_D*S_B*p_b + P_v_0*S_D*p_g)/(S_D*(P_D - P_v_0))
        result.append(p_v_max)
        return result

