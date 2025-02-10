from ..kwasak import kwasak_staticclass SelectingPump:

    @kwasak_static
    def eqn_8_01(NC: float = None, NS: float = None, SCON: float = None, installation_cost: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_01__NC(NS: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return result

    @staticmethod
    def eqn_8_01__NS(NC: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return result

    @staticmethod
    def eqn_8_01__SCON(NC: float, NS: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_01__installation_cost(NC: float, NS: float, SCON: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        installation_cost = 1426.00150101399*SCON**(7/20)*(2.0*NC + NS)
        result.append(installation_cost)
        return result

    @kwasak_static
    def eqn_8_02(hp: float = None, installed_costs: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_02__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        hp = 9.18273645546364e-9*installed_costs**2
        result.append(hp)
        return result

    @staticmethod
    def eqn_8_02__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        installed_costs = 10435.5162785557*sqrt(hp)
        result.append(installed_costs)
        return result

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

    @kwasak_static
    def eqn_8_05(Eff: float = None, actual_brake_horsepower: float = None, theoretical_adiabatic_horsepower: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_05__Eff(actual_brake_horsepower: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        Eff = theoretical_adiabatic_horsepower/actual_brake_horsepower
        result.append(Eff)
        return result

    @staticmethod
    def eqn_8_05__actual_brake_horsepower(Eff: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        actual_brake_horsepower = theoretical_adiabatic_horsepower/Eff
        result.append(actual_brake_horsepower)
        return result

    @staticmethod
    def eqn_8_05__theoretical_adiabatic_horsepower(Eff: float, actual_brake_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        theoretical_adiabatic_horsepower = Eff*actual_brake_horsepower
        result.append(theoretical_adiabatic_horsepower)
        return result

    @kwasak_static
    def eqn_8_06(M: float = None, P_1: float = None, P_2: float = None, R: float = None, T: float = None, adiabatic_hp: float = None, k: float = None, w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_06__M(P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        M = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*adiabatic_hp*(k - 1))
        result.append(M)
        return result

    @staticmethod
    def eqn_8_06__P_1(M: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_1 = P_2/(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_1)
        return result

    @staticmethod
    def eqn_8_06__P_2(M: float, P_1: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_2 = P_1*(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_2)
        return result

    @staticmethod
    def eqn_8_06__R(M: float, P_1: float, P_2: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        R = 1980000*M*adiabatic_hp*(k - 1)/(T*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(R)
        return result

    @staticmethod
    def eqn_8_06__T(M: float, P_1: float, P_2: float, R: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        T = 1980000*M*adiabatic_hp*(k - 1)/(R*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(T)
        return result

    @staticmethod
    def eqn_8_06__adiabatic_hp(M: float, P_1: float, P_2: float, R: float, T: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        adiabatic_hp = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*M*(k - 1))
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_06__k(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_06__w(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        w = 1980000*M*adiabatic_hp*(k - 1)/(R*T*k*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(w)
        return result

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

    @kwasak_static
    def eqn_8_08(P_1: float = None, P_2: float = None, adiabatic_power_watts: float = None, f: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_08__P_1(P_2: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_08__P_2(P_1: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_08__adiabatic_power_watts(P_1: float, P_2: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_power_watts = 0.0833333333333333*f*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_power_watts)
        return result

    @staticmethod
    def eqn_8_08__f(P_1: float, P_2: float, adiabatic_power_watts: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        f = 12.0*adiabatic_power_watts/((P_2/P_1)**0.286 - 1.0)
        result.append(f)
        return result

    @kwasak_static
    def eqn_8_09(E_j: float = None, E_m: float = None, e: float = None, r: float = None, s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_09__E_j(E_m: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_j = 0.341296928327645*E_m*r*s/e
        result.append(E_j)
        return result

    @staticmethod
    def eqn_8_09__E_m(E_j: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_m = 2.93*E_j*e/(r*s)
        result.append(E_m)
        return result

    @staticmethod
    def eqn_8_09__e(E_j: float, E_m: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        e = 0.341296928327645*E_m*r*s/E_j
        result.append(e)
        return result

    @staticmethod
    def eqn_8_09__r(E_j: float, E_m: float, e: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        r = 2.93*E_j*e/(E_m*s)
        result.append(r)
        return result

    @staticmethod
    def eqn_8_09__s(E_j: float, E_m: float, e: float, r: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        s = 2.93*E_j*e/(E_m*r)
        result.append(s)
        return result


