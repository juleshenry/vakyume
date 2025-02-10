from ..kwasak import kwasak_static
    @kwasak_static
    def eqn_7_04ab(P_c: float = None, p: float = None, p_i: float = None, p_nc: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_04ab__P_c(p: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        P_c = p - p_nc
        result.append(P_c)
        return result

    @staticmethod
    def eqn_7_04ab__p(P_c: float, p_i: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p = P_c + p_nc
        result.append(p)
        return result

    @staticmethod
    def eqn_7_04ab__p_i(P_c: float, p: float, p_nc: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = 0
        result.append(p_i)
        return result

    @staticmethod
    def eqn_7_04ab__p_nc(P_c: float, p: float, p_i: float):
        # [.pyeqn] p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = -P_c + p
        result.append(p_nc)
        return result

    @kwasak_static
    def eqn_7_14b(A: float = None, Q_condensor_heat_duty: float = None, U: float = None, del_T_1: float = None, del_T_2: float = None,**kwargs):
        return


    @staticmethod
    def eqn_7_14b__A(Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(A)
        return result

    @staticmethod
    def eqn_7_14b__Q_condensor_heat_duty(A: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A*U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
        return result

    @staticmethod
    def eqn_7_14b__U(A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty/(A*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(U)
        return result

    @staticmethod
    def eqn_7_14b__del_T_1(A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_1)
        return result

    @staticmethod
    def eqn_7_14b__del_T_2(A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_2)
        return result

    @kwasak_static