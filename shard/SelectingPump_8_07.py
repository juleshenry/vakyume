from kwasak import kwasak_static

class SelectingPump:

    @kwasak_static
    def eqn_8_07(P_1: float = None, P_2: float = None, adiabatic_hp: float = None, w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_07__P_1(P_2: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_07__P_2(P_1: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_07__adiabatic_hp(P_1: float, P_2: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_hp = 0.05*w*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_07__w(P_1: float, P_2: float, adiabatic_hp: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        w = 20.0*adiabatic_hp/((P_2/P_1)**0.286 - 1.0)
        result.append(w)
        return result
