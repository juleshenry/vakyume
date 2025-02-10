from ..kwasak import kwasak_static
    @kwasak_static
    def eqn_10_10(bhp: float = None, bhp_0: float = None, mu: float = None, rho: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_10__bhp(bhp_0: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp = 0.0005*bhp_0*(31.0*mu**(4/25)*rho**(21/25) + 1000.0)
        result.append(bhp)
        return result

    @staticmethod
    def eqn_10_10__bhp_0(bhp: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp_0 = 2000.0*bhp/(31.0*mu**0.16*rho**0.84 + 1000.0)
        result.append(bhp_0)
        return result

    @staticmethod
    def eqn_10_10__mu(bhp: float, bhp_0: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        mu = -204374584201.104*I*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = 204374584201.104*I*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = -204374584201.104*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        mu = 204374584201.104*(bhp/(bhp_0*rho**0.84) - 0.5/rho**0.84)**(25/4)
        result.append(mu)
        return result

    @staticmethod
    def eqn_10_10__rho(bhp: float, bhp_0: float, mu: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        # [Sympy Failover]
        pass # Ollama offline

    @kwasak_static
    def eqn_10_19(P: float = None, S_Th: float = None, S_p: float = None, T_e: float = None, T_i: float = None, p_c: float = None, p_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_19__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_10_19__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_Th = S_p/((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_19__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_p = S_Th*((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_19__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_e = (P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P - T_i*p_s + 460.0*p_c*(S_p/S_Th)**(5/3) - 460.0*p_s)/((S_p/S_Th)**(5/3)*(P - p_c))
        result.append(T_e)
        T_e = 2.05280095711867*(P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        T_e = 2.05280095711867*(P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        return result

    @staticmethod
    def eqn_10_19__T_i(P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_i = (P*T_e*(S_p/S_Th)**(5/3) + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P - T_e*p_c*(S_p/S_Th)**(5/3) - 460.0*p_c*(S_p/S_Th)**(5/3) + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        return result

    @staticmethod
    def eqn_10_19__p_c(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_c = (P*T_e*(S_p/S_Th)**(5/3) - P*T_i + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P + T_i*p_s + 460.0*p_s)/((S_p/S_Th)**(5/3)*(T_e + 460.0))
        result.append(p_c)
        p_c = 2.05280095711867*(0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        p_c = 2.05280095711867*(0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        return result

    @staticmethod
    def eqn_10_19__p_s(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_s = (-P*T_e*(S_p/S_Th)**(5/3) + P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P + T_e*p_c*(S_p/S_Th)**(5/3) + 460.0*p_c*(S_p/S_Th)**(5/3))/(T_i + 460.0)
        result.append(p_s)
        p_s = (-0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        p_s = (-0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        return result

    @kwasak_static
    def eqn_10_20(P: float = None, S_0: float = None, S_p: float = None, T_e: float = None, T_i: float = None, p_0: float = None, p_c: float = None, p_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_20__P(S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_10_20__S_0(P: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_0 = S_p*((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_20__S_p(P: float, S_0: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_p = S_0/((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_20__T_e(P: float, S_0: float, S_p: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_10_20__T_i(P: float, S_0: float, S_p: float, T_e: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_10_20__p_0(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_10_20__p_c(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_10_20__p_s(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @kwasak_static