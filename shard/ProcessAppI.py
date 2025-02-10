from kwasak import kwasak_static

class ProcessAppI:

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
