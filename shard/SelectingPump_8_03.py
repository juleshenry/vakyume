from kwasak import kwasak_static

class SelectingPump:

    @kwasak_static
    def eqn_8_03(hp: float = None, installed_costs: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_03__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_03__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        result = []
        installed_costs = 13482.9087908759*hp**(9/20)
        result.append(installed_costs)
        return result
