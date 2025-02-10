from ..kwasak import kwasak_staticclass LiquidRing:

    @kwasak_static
    def eqn_10_01(D_r: float = None, sig_R: float = None, w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_01__D_r(sig_R: float, w: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        D_r = 229.357798165138*sig_R/w
        result.append(D_r)
        return result

    @staticmethod
    def eqn_10_01__sig_R(D_r: float, w: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        sig_R = 0.00436*D_r*w
        result.append(sig_R)
        return result

    @staticmethod
    def eqn_10_01__w(D_r: float, sig_R: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        w = 229.357798165138*sig_R/D_r
        result.append(w)
        return result

    @kwasak_static
    def eqn_10_02(PS: float = None, Q_gas: float = None, V: float = None, dP: float = None, dt: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_02__PS(Q_gas: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        PS = Q_gas - V*dP/dt
        result.append(PS)
        return result

    @staticmethod
    def eqn_10_02__Q_gas(PS: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        Q_gas = PS + V*dP/dt
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_02__V(PS: float, Q_gas: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        V = dt*(-PS + Q_gas)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_10_02__dP(PS: float, Q_gas: float, V: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dP = dt*(-PS + Q_gas)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_10_02__dt(PS: float, Q_gas: float, V: float, dP: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dt = -V*dP/(PS - Q_gas)
        result.append(dt)
        return result

    @kwasak_static
    def eqn_10_03(N_mfw: float = None, Q_gas: float = None, T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_03__N_mfw(Q_gas: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        N_mfw = 0.108108108108108*Q_gas/T
        result.append(N_mfw)
        return result

    @staticmethod
    def eqn_10_03__Q_gas(N_mfw: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        Q_gas = 9.25*N_mfw*T
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_03__T(N_mfw: float, Q_gas: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        T = 0.108108108108108*Q_gas/N_mfw
        result.append(T)
        return result

    @kwasak_static
    def eqn_10_04(Q_gas: float = None, SP_1: float = None, SP_2: float = None, S_p: float = None, V: float = None, t: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_04__Q_gas(SP_1: float, SP_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        Q_gas = -(SP_1 - SP_2*exp(S_p*t/V))/(exp(S_p*t/V) - 1)
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_04__SP_1(Q_gas: float, SP_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_1 = Q_gas + (-Q_gas + SP_2)*exp(S_p*t/V)
        result.append(SP_1)
        return result

    @staticmethod
    def eqn_10_04__SP_2(Q_gas: float, SP_1: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_2 = (Q_gas*exp(S_p*t/V) - Q_gas + SP_1)*exp(-S_p*t/V)
        result.append(SP_2)
        return result

    @staticmethod
    def eqn_10_04__S_p(Q_gas: float, SP_1: float, SP_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        S_p = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_04__V(Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        V = S_p*t/log((Q_gas - SP_1)/(Q_gas - SP_2))
        result.append(V)
        return result

    @staticmethod
    def eqn_10_04__t(Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        t = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/S_p
        result.append(t)
        return result

    @kwasak_static
    def eqn_10_05(P_1: float = None, P_2: float = None, S_p: float = None, V: float = None, t: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_05__P_1(P_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_p*t/V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_05__P_2(P_1: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_p*t/V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_05__S_p(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        S_p = V*log(P_1/P_2)/t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_05__V(P_1: float, P_2: float, S_p: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        V = S_p*t/log(P_1/P_2)
        result.append(V)
        return result

    @staticmethod
    def eqn_10_05__t(P_1: float, P_2: float, S_p: float, V: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_p
        result.append(t)
        return result

    @kwasak_static
    def eqn_10_06(P_1: float = None, P_2: float = None, S_a: float = None, V: float = None, t: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_06__P_1(P_2: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_a*t/V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_06__P_2(P_1: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_a*t/V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_06__S_a(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        S_a = V*log(P_1/P_2)/t
        result.append(S_a)
        return result

    @staticmethod
    def eqn_10_06__V(P_1: float, P_2: float, S_a: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        V = S_a*t/log(P_1/P_2)
        result.append(V)
        return result

    @staticmethod
    def eqn_10_06__t(P_1: float, P_2: float, S_a: float, V: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_a
        result.append(t)
        return result

    @kwasak_static
    def eqn_10_08(bhp: float = None, c_p: float = None, delta_T: float = None, delta_h_i: float = None, f_a: float = None, rho: float = None, w_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_08__bhp(c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        bhp = 0.00315127701375246*c_p*delta_T*f_a*rho - 0.000392927308447937*delta_h_i*w_i
        result.append(bhp)
        return result

    @staticmethod
    def eqn_10_08__c_p(bhp: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        c_p = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(delta_T*f_a*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_10_08__delta_T(bhp: float, c_p: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_T = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*f_a*rho)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_10_08__delta_h_i(bhp: float, c_p: float, delta_T: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_h_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/w_i
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_10_08__f_a(bhp: float, c_p: float, delta_T: float, delta_h_i: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        f_a = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*rho)
        result.append(f_a)
        return result

    @staticmethod
    def eqn_10_08__rho(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        rho = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*f_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_10_08__w_i(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        w_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/delta_h_i
        result.append(w_i)
        return result

    @kwasak_static
    def eqn_10_09(T_c: float = None, T_s: float = None, delta_T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_09__T_c(T_s: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_c = T_s + delta_T
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_09__T_s(T_c: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_s = T_c - delta_T
        result.append(T_s)
        return result

    @staticmethod
    def eqn_10_09__delta_T(T_c: float, T_s: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        delta_T = T_c - T_s
        result.append(delta_T)
        return result

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
    def eqn_10_11(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_11__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_c = T_s + 10
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_11__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_s = T_c - 10
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_12(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_12__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_c = T_s + 5
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_12__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_s = T_c - 5
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_13(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_13__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_c = T_s + 25
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_13__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_s = T_c - 25
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_14(T_c: float = None, T_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_14__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_c = T_s + 12
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_14__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_s = T_c - 12
        result.append(T_s)
        return result

    @kwasak_static
    def eqn_10_15(P: float = None, S_Th: float = None, S_p: float = None, p_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_15__P(S_Th: float, S_p: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        P = S_Th*p_s/(S_Th - S_p)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_15__S_Th(P: float, S_p: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_Th = P*S_p/(P - p_s)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_15__S_p(P: float, S_Th: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        S_p = S_Th*(P - p_s)/P
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_15__p_s(P: float, S_Th: float, S_p: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        p_s = P*(S_Th - S_p)/S_Th
        result.append(p_s)
        return result

    @kwasak_static
    def eqn_10_16(P: float = None, S_0: float = None, S_Th: float = None, p_0: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_16__P(S_0: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        P = p_0*(S_Th/S_0)**(5/3)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_16__S_0(P: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/(P/(P - p_0))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_16__S_Th(P: float, S_0: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*(P/(P - p_0))**(3/5)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_16__p_0(P: float, S_0: float, S_Th: float):
        # [.pyeqn] S_Th = S_0 * (P / (P - p_0)) ** 0.6
        result = []
        p_0 = P - P/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = P - 2.05280095711867*P/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = P - 2.05280095711867*P/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result

    @kwasak_static
    def eqn_10_17(P: float = None, S_0: float = None, S_Th: float = None, p_0: float = None, p_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_17__P(S_0: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        P = (p_0*(S_Th/S_0)**(5/3) - p_s)/((S_Th/S_0)**1.66666666666667 - 1.0)
        result.append(P)
        P = (0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        P = (0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - p_s)/(0.487139289628747*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - 1.0)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_17__S_0(P: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/((P - p_s)/(P - p_0))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_17__S_Th(P: float, S_0: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*((P - p_s)/(P - p_0))**(3/5)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_17__p_0(P: float, S_0: float, S_Th: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_0 = (P*(S_Th/S_0)**(5/3) - P + p_s)/(S_Th/S_0)**(5/3)
        result.append(p_0)
        p_0 = 2.05280095711867*(0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        p_0 = 2.05280095711867*(0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 - P + p_s)/(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_0)
        return result

    @staticmethod
    def eqn_10_17__p_s(P: float, S_0: float, S_Th: float, p_0: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        p_s = -P*(S_Th/S_0)**(5/3) + P + p_0*(S_Th/S_0)**(5/3)
        result.append(p_s)
        p_s = -0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5 + P + 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 - I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        p_s = -0.487139289628747*P*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5 + P + 0.487139289628747*p_0*(-0.577350269189626*(S_Th/S_0)**0.333333333333333 + I*(S_Th/S_0)**0.333333333333333)**5
        result.append(p_s)
        return result

    @kwasak_static
    def eqn_10_18(P: float = None, S_Th: float = None, S_p: float = None, T_e: float = None, T_i: float = None, p_c: float = None, p_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_18__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        P = (S_Th*T_i*p_s + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*T_i + 460*S_Th - S_p*T_e - 460*S_p)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_18__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_Th = S_p*(P*T_e + 460*P - T_e*p_c - 460*p_c)/(P*T_i + 460*P - T_i*p_s - 460*p_s)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_18__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_p = S_Th*(P*T_i + 460*P - T_i*p_s - 460*p_s)/(P*T_e + 460*P - T_e*p_c - 460*p_c)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_18__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_e = (P*S_Th*T_i + 460*P*S_Th - 460*P*S_p - S_Th*T_i*p_s - 460*S_Th*p_s + 460*S_p*p_c)/(S_p*(P - p_c))
        result.append(T_e)
        return result

    @staticmethod
    def eqn_10_18__T_i(P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_i = (-460*P*S_Th + P*S_p*T_e + 460*P*S_p + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*(P - p_s))
        result.append(T_i)
        return result

    @staticmethod
    def eqn_10_18__p_c(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_c = (-P*S_Th*T_i - 460*P*S_Th + P*S_p*T_e + 460*P*S_p + S_Th*T_i*p_s + 460*S_Th*p_s)/(S_p*(T_e + 460))
        result.append(p_c)
        return result

    @staticmethod
    def eqn_10_18__p_s(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_s = (P*S_Th*T_i + 460*P*S_Th - P*S_p*T_e - 460*P*S_p + S_p*T_e*p_c + 460*S_p*p_c)/(S_Th*(T_i + 460))
        result.append(p_s)
        return result

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
    def eqn_10_21(P: float = None, P_d: float = None, P_prime: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_21__P(P_d: float, P_prime: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P = P_d*P_prime/760
        result.append(P)
        return result

    @staticmethod
    def eqn_10_21__P_d(P: float, P_prime: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_d = 760*P/P_prime
        result.append(P_d)
        return result

    @staticmethod
    def eqn_10_21__P_prime(P: float, P_d: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_prime = 760*P/P_d
        result.append(P_prime)
        return result


