from ..kwasak import kwasak_staticclass ProcessAppI:

    @kwasak_static
    def eqn_5_01(K_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_01__K_i(x_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i/x_i
        result.append(K_i)
        return result

    @staticmethod
    def eqn_5_01__x_i(K_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_01__y_i(K_i: float, x_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_02a(K_1: float = None, K_2: float = None, alpha_1_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02a__K_1(K_2: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_1 = K_2*alpha_1_2
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02a__K_2(K_1: float, alpha_1_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        K_2 = K_1/alpha_1_2
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02a__alpha_1_2(K_1: float, K_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2
        result = []
        alpha_1_2 = K_1/K_2
        result.append(alpha_1_2)
        return result

    @kwasak_static
    def eqn_5_02b(K_1: float = None, K_2: float = None, x_1: float = None, x_2: float = None, y_1: float = None, y_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_02b__K_1(K_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*x_2*y_1/(x_1*y_2)
        result.append(K_1)
        return result

    @staticmethod
    def eqn_5_02b__K_2(K_1: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1*x_1*y_2/(x_2*y_1)
        result.append(K_2)
        return result

    @staticmethod
    def eqn_5_02b__x_1(K_1: float, K_2: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_1 = K_2*x_2*y_1/(K_1*y_2)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_02b__x_2(K_1: float, K_2: float, x_1: float, y_1: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        x_2 = K_1*x_1*y_2/(K_2*y_1)
        result.append(x_2)
        return result

    @staticmethod
    def eqn_5_02b__y_1(K_1: float, K_2: float, x_1: float, x_2: float, y_2: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_1 = K_1*x_1*y_2/(K_2*x_2)
        result.append(y_1)
        return result

    @staticmethod
    def eqn_5_02b__y_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float):
        # [.pyeqn] K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        y_2 = K_2*x_2*y_1/(K_1*x_1)
        result.append(y_2)
        return result

    @kwasak_static
    def eqn_5_03(P_0_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_03__P_0_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_03__p_i(P_0_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_03__x_i(P_0_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_04(P: float = None, P_0_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_04__P(P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_04__P_0_i(P: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P_0_i = P*y_i/x_i
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_04__x_i(P: float, P_0_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        x_i = P*y_i/P_0_i
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_04__y_i(P: float, P_0_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_05(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_05__P_0_1(P_0_2: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2*alpha_12
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_05__P_0_2(P_0_1: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1/alpha_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_05__alpha_12(P_0_1: float, P_0_2: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return result

    @kwasak_static
    def eqn_5_06(P_0_i: float = None, gamma_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_06__P_0_i(gamma_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_06__gamma_i(P_0_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_06__p_i(P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_06__x_i(P_0_i: float, gamma_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_07(P: float = None, P_0_i: float = None, gamma_i: float = None, x_i: float = None, y_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_07__P(P_0_i: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i*gamma_i*x_i/y_i
        result.append(P)
        return result

    @staticmethod
    def eqn_5_07__P_0_i(P: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P*y_i/(gamma_i*x_i)
        result.append(P_0_i)
        return result

    @staticmethod
    def eqn_5_07__gamma_i(P: float, P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P*y_i/(P_0_i*x_i)
        result.append(gamma_i)
        return result

    @staticmethod
    def eqn_5_07__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return result

    @staticmethod
    def eqn_5_07__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return result

    @kwasak_static
    def eqn_5_08(P_0_1: float = None, P_0_2: float = None, alpha_12: float = None, gamma_1: float = None, gamma_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_08__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_08__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_08__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return result

    @staticmethod
    def eqn_5_08__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return result

    @staticmethod
    def eqn_5_08__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return result

    @kwasak_static
    def eqn_5_09(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_09__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_09__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_09__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10a(D: float = None, L_0: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10a__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10a__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10a__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10b(L_0: float = None, R: float = None, V_1: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10b__L_0(R: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10b__R(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return result

    @staticmethod
    def eqn_5_10b__V_1(L_0: float, R: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return result

    @kwasak_static
    def eqn_5_10c(D: float = None, L_0: float = None, R: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_10c__D(L_0: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        D = L_0/R
        result.append(D)
        return result

    @staticmethod
    def eqn_5_10c__L_0(D: float, R: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        L_0 = D*R
        result.append(L_0)
        return result

    @staticmethod
    def eqn_5_10c__R(D: float, L_0: float):
        # [.pyeqn] (L_0 / D) / (L_0 / D + 1) = R / (R + 1)
        result = []
        R = L_0/D
        result.append(R)
        return result

    @kwasak_static
    def eqn_5_11(B: float = None, L_N: float = None, V_0: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_11__B(L_N: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        B = L_N - V_0
        result.append(B)
        return result

    @staticmethod
    def eqn_5_11__L_N(B: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        L_N = B + V_0
        result.append(L_N)
        return result

    @staticmethod
    def eqn_5_11__V_0(B: float, L_N: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        V_0 = -B + L_N
        result.append(V_0)
        return result

    @kwasak_static
    def eqn_5_12(Eff: float = None, N_ES: float = None, N_t: float = None, T: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_12__Eff(N_ES: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        Eff = (N_ES/N_t)**(1/T)
        result.append(Eff)
        return result

    @staticmethod
    def eqn_5_12__N_ES(Eff: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_ES = Eff**T*N_t
        result.append(N_ES)
        return result

    @staticmethod
    def eqn_5_12__N_t(Eff: float, N_ES: float, T: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        N_t = N_ES/Eff**T
        result.append(N_t)
        return result

    @staticmethod
    def eqn_5_12__T(Eff: float, N_ES: float, N_t: float):
        # [.pyeqn] N_t = N_ES / Eff ** T
        result = []
        T = log(N_ES/N_t)/log(Eff)
        result.append(T)
        return result

    @kwasak_static
    def eqn_5_13(HETP: float = None, H_p: float = None, N_ES: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_13__HETP(H_p: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return result

    @staticmethod
    def eqn_5_13__H_p(HETP: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return result

    @staticmethod
    def eqn_5_13__N_ES(HETP: float, H_p: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return result

    @kwasak_static
    def eqn_5_14(M: float = None, P_0: float = None, T: float = None, W_E: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_14__M(P_0: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return result

    @staticmethod
    def eqn_5_14__P_0(M: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return result

    @staticmethod
    def eqn_5_14__T(M: float, P_0: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return result

    @staticmethod
    def eqn_5_14__W_E(M: float, P_0: float, T: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return result

    @kwasak_static
    def eqn_5_15(M_1: float = None, M_2: float = None, P_0_1: float = None, P_0_2: float = None, a_M_12: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_15__M_1(M_2: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        M_1 = M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        return result

    @staticmethod
    def eqn_5_15__M_2(M_1: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_2 = -M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        M_2 = M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        return result

    @staticmethod
    def eqn_5_15__P_0_1(M_1: float, M_2: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
        result.append(P_0_1)
        return result

    @staticmethod
    def eqn_5_15__P_0_2(M_1: float, M_2: float, P_0_1: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1*(M_2/M_1)**(2/5)/a_M_12
        result.append(P_0_2)
        return result

    @staticmethod
    def eqn_5_15__a_M_12(M_1: float, M_2: float, P_0_1: float, P_0_2: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
        result.append(a_M_12)
        return result

    @kwasak_static
    def eqn_5_16(H_i: float = None, p_i: float = None, x_i: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_16__H_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return result

    @staticmethod
    def eqn_5_16__p_i(H_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return result

    @staticmethod
    def eqn_5_16__x_i(H_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return result

    @kwasak_static
    def eqn_5_17(H_2_1: float = None, H_2_3: float = None, H_2_mi: float = None, x_1: float = None, x_3: float = None,**kwargs):
        return


    @staticmethod
    def eqn_5_17__H_2_1(H_2_3: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_1 = exp((-x_3*log(H_2_3) + log(H_2_mi))/x_1)
        result.append(H_2_1)
        return result

    @staticmethod
    def eqn_5_17__H_2_3(H_2_1: float, H_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_3 = exp((-x_1*log(H_2_1) + log(H_2_mi))/x_3)
        result.append(H_2_3)
        return result

    @staticmethod
    def eqn_5_17__H_2_mi(H_2_1: float, H_2_3: float, x_1: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        H_2_mi = exp(x_1*log(H_2_1) + x_3*log(H_2_3))
        result.append(H_2_mi)
        return result

    @staticmethod
    def eqn_5_17__x_1(H_2_1: float, H_2_3: float, H_2_mi: float, x_3: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return result

    @staticmethod
    def eqn_5_17__x_3(H_2_1: float, H_2_3: float, H_2_mi: float, x_1: float):
        # [.pyeqn] log(H_2_mi) = x_1 * log(H_2_1) + x_3 * log(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return result


