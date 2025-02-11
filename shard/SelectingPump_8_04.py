from kwasak import kwasak_static

class SelectingPump:

    @kwasak_static
    def eqn_8_04(hp: float = None, installed_costs: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_04__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        hp = -9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        hp = 9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        return result

    @staticmethod
    def eqn_8_04__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        installed_costs = 10350.7864343909*hp**(2/5)
        result.append(installed_costs)
        return result
