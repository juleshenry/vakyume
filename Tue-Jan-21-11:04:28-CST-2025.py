from math import log, sqrt, exp
from sympy import I, Piecewise, LambertW, Eq


class SelectingPump:

    @staticmethod
    def eqn_8_1__SCON(NC: float, NS: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        # [FAILED TO PARSE VIA SYMPY]
        result = []
        return [ None]
        denominator = (16000 * (NS + 2*NC) / 1000) ** 0.35 - installation_cost
        return [ None]
        return [ SCON]

    @staticmethod
    def eqn_8_1__installation_cost(NC: float, NS: float, SCON: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        installation_cost = 1426.00150101399*SCON**(7/20)*(2.0*NC + NS)
        result.append(installation_cost)
        return result

    @staticmethod
    def eqn_8_1__NC(NS: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return result

    @staticmethod
    def eqn_8_1__NS(NC: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return result

    @staticmethod
    def eqn_8_2__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        hp = 9.18273645546364e-9*installed_costs**2
        result.append(hp)
        return result

    @staticmethod
    def eqn_8_2__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 33000 * (hp / 10) ** 0.5
        result = []
        installed_costs = 10435.5162785557*sqrt(hp)
        result.append(installed_costs)
        return result

    @staticmethod
    def eqn_8_3__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        # [FAILED TO PARSE VIA SYMPY]
        return [ ((installed_costs / 38000 ** 0.45) ** (1 / 0.45))]

    @staticmethod
    def eqn_8_3__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 38000 * (hp / 10) ** 0.45
        result = []
        installed_costs = 13482.9087908759*hp**(9/20)
        result.append(installed_costs)
        return result

    @staticmethod
    def eqn_8_4__hp(installed_costs: float):
        # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        hp = -9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        hp = 9.1741667595569e-11*installed_costs**(5/2)
        result.append(hp)
        return result

    @staticmethod
    def eqn_8_4__installed_costs(hp: float):
        # [.pyeqn] installed_costs = 26000 * (hp / 10) ** 0.4
        result = []
        installed_costs = 10350.7864343909*hp**(2/5)
        result.append(installed_costs)
        return result

    @staticmethod
    def eqn_8_5__Eff(actual_brake_horsepower: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        Eff = theoretical_adiabatic_horsepower/actual_brake_horsepower
        result.append(Eff)
        return result

    @staticmethod
    def eqn_8_5__actual_brake_horsepower(Eff: float, theoretical_adiabatic_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        actual_brake_horsepower = theoretical_adiabatic_horsepower/Eff
        result.append(actual_brake_horsepower)
        return result

    @staticmethod
    def eqn_8_5__theoretical_adiabatic_horsepower(Eff: float, actual_brake_horsepower: float):
        # [.pyeqn] Eff = theoretical_adiabatic_horsepower / actual_brake_horsepower
        result = []
        theoretical_adiabatic_horsepower = Eff*actual_brake_horsepower
        result.append(theoretical_adiabatic_horsepower)
        return result

    @staticmethod
    def eqn_8_6__w(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        w = 1980000*M*adiabatic_hp*(k - 1)/(R*T*k*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(w)
        return result

    @staticmethod
    def eqn_8_6__adiabatic_hp(M: float, P_1: float, P_2: float, R: float, T: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        adiabatic_hp = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*M*(k - 1))
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_6__P_1(M: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_1 = P_2/(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_1)
        return result

    @staticmethod
    def eqn_8_6__T(M: float, P_1: float, P_2: float, R: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        T = 1980000*M*adiabatic_hp*(k - 1)/(R*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(T)
        return result

    @staticmethod
    def eqn_8_6__R(M: float, P_1: float, P_2: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        R = 1980000*M*adiabatic_hp*(k - 1)/(T*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(R)
        return result

    @staticmethod
    def eqn_8_6__M(P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        M = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*adiabatic_hp*(k - 1))
        result.append(M)
        return result

    @staticmethod
    def eqn_8_6__P_2(M: float, P_1: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_2 = P_1*(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_2)
        return result

    @staticmethod
    def eqn_8_6__k(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        # [FAILED TO PARSE VIA SYMPY]
        term = (w * M * R * T) / ((P_2 / P_1) ** ((k - 1) / k) - 1)
        closed_form_expr = lambda adiabatic_hp: solve(0 == (adiabatic_hp + term * ((M * 550 * 3600) / (w * R * T))), k)
        return [ closed_form_expr]

    @staticmethod
    def eqn_8_7__w(P_1: float, P_2: float, adiabatic_hp: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        w = 20.0*adiabatic_hp/((P_2/P_1)**0.286 - 1.0)
        result.append(w)
        return result

    @staticmethod
    def eqn_8_7__P_1(P_2: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # [FAILED TO PARSE VIA SYMPY]
        zero = 0
        denominator = (w / 20) * ((P_2 / Pe**0.286) - 1)
        isolated_P_1_term = zero + adiabatic_hp / denominator
        P_1 = exp(isolated_P_1_term) * Pe**0.286
        return [ P_1]

    @staticmethod
    def eqn_8_7__P_2(P_1: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        # [FAILED TO PARSE VIA SYMPY]
        return [ ((w / 20) * (pow((P_2 / P_1), 1/0.286) - 1)) + adiabatic_hp == 0]
        return [ exp((adiabatic_hp - (w / 20)) / ((pow(P_1, 1/0.286) - 1))) * P_1]

    @staticmethod
    def eqn_8_7__adiabatic_hp(P_1: float, P_2: float, w: float):
        # [.pyeqn] adiabatic_hp = (w / 20) * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_hp = 0.05*w*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_8__f(P_1: float, P_2: float, adiabatic_power_watts: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        f = 12.0*adiabatic_power_watts/((P_2/P_1)**0.286 - 1.0)
        result.append(f)
        return result

    @staticmethod
    def eqn_8_8__adiabatic_power_watts(P_1: float, P_2: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        result = []
        adiabatic_power_watts = 0.0833333333333333*f*((P_2/P_1)**(143/500) - 1.0)
        result.append(adiabatic_power_watts)
        return result

    @staticmethod
    def eqn_8_8__P_1(P_2: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [FAILED TO PARSE VIA SYMPY]
        P_1 = symbols('P_1')
        equation = Eq(0, (f / 12) * ((P_2 / P_1) ** 0.286 - 1) - adiabatic_power_watts)
        solution = solve(equation, P_1)
        return solution[0]

    @staticmethod
    def eqn_8_8__P_2(P_1: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [FAILED TO PARSE VIA SYMPY]
        P_2 = (adiabatic_power_watts / (f / 12 * ((P_1 / P_2) ** 0.286 - 1))) + P_1
        return [ P_2]

    @staticmethod
    def eqn_8_9__e(E_j: float, E_m: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        e = 0.341296928327645*E_m*r*s/E_j
        result.append(e)
        return result

    @staticmethod
    def eqn_8_9__E_j(E_m: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_j = 0.341296928327645*E_m*r*s/e
        result.append(E_j)
        return result

    @staticmethod
    def eqn_8_9__r(E_j: float, E_m: float, e: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        r = 2.93*E_j*e/(E_m*s)
        result.append(r)
        return result

    @staticmethod
    def eqn_8_9__E_m(E_j: float, e: float, r: float, s: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        E_m = 2.93*E_j*e/(r*s)
        result.append(E_m)
        return result

    @staticmethod
    def eqn_8_9__s(E_j: float, E_m: float, e: float, r: float):
        # [.pyeqn] r = 2.93 * (E_j * e) / (E_m * s)
        result = []
        s = 2.93*E_j*e/(E_m*r)
        result.append(s)
        return result


class SteamJetInjectors:

    @staticmethod
    def eqn_9_1__A(rho_s: float, v: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        A = w_s/(rho_s*v)
        result.append(A)
        return result

    @staticmethod
    def eqn_9_1__v(A: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        v = w_s/(A*rho_s)
        result.append(v)
        return result

    @staticmethod
    def eqn_9_1__rho_s(A: float, v: float, w_s: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        rho_s = w_s/(A*v)
        result.append(rho_s)
        return result

    @staticmethod
    def eqn_9_1__w_s(A: float, rho_s: float, v: float):
        # [.pyeqn] w_s = v * A * rho_s
        result = []
        w_s = A*rho_s*v
        result.append(w_s)
        return result

    @staticmethod
    def eqn_9_2__d_n(P_m: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        d_n = -0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        d_n = 0.0339853079285911*sqrt(w_s/(P_m*rho_s)**0.5)
        result.append(d_n)
        return result

    @staticmethod
    def eqn_9_2__P_m(d_n: float, rho_s: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        P_m = 1.334027668054e-6*w_s**2/(d_n**4*rho_s)
        result.append(P_m)
        return result

    @staticmethod
    def eqn_9_2__rho_s(P_m: float, d_n: float, w_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        rho_s = 1.334027668054e-6*w_s**2/(P_m*d_n**4)
        result.append(rho_s)
        return result

    @staticmethod
    def eqn_9_2__w_s(P_m: float, d_n: float, rho_s: float):
        # [.pyeqn] w_s = 865.8 * d_n**2 * (P_m * rho_s) ** 0.5
        result = []
        w_s = 865.8*d_n**2*sqrt(P_m*rho_s)
        result.append(w_s)
        return result

    @staticmethod
    def eqn_9_3__t_e(P_s: float, V: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        t_e = 0.001*V*(2300.0 - 3.0*P_s)/w_j
        result.append(t_e)
        return result

    @staticmethod
    def eqn_9_3__w_j(P_s: float, V: float, t_e: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        w_j = 0.001*V*(2300.0 - 3.0*P_s)/t_e
        result.append(w_j)
        return result

    @staticmethod
    def eqn_9_3__P_s(V: float, t_e: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        P_s = 33.3333333333333*(23.0*V - 10.0*t_e*w_j)/V
        result.append(P_s)
        return result

    @staticmethod
    def eqn_9_3__V(P_s: float, t_e: float, w_j: float):
        # [.pyeqn] t_e = (2.3 - 0.003 * P_s) * V / w_j
        result = []
        V = -1000.0*t_e*w_j/(3.0*P_s - 2300.0)
        result.append(V)
        return result

    @staticmethod
    def eqn_9_4__w_s(AEL: float, SC: float, r: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        w_s = AEL*SC*r
        result.append(w_s)
        return result

    @staticmethod
    def eqn_9_4__r(AEL: float, SC: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        r = w_s/(AEL*SC)
        result.append(r)
        return result

    @staticmethod
    def eqn_9_4__SC(AEL: float, r: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        SC = w_s/(AEL*r)
        result.append(SC)
        return result

    @staticmethod
    def eqn_9_4__AEL(SC: float, r: float, w_s: float):
        # [.pyeqn] w_s = AEL * r * SC
        result = []
        AEL = w_s/(SC*r)
        result.append(AEL)
        return result

    @staticmethod
    def eqn_9_5__t_h(V: float, r_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        t_h = V*r_h/w_h
        result.append(t_h)
        return result

    @staticmethod
    def eqn_9_5__r_h(V: float, t_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        r_h = t_h*w_h/V
        result.append(r_h)
        return result

    @staticmethod
    def eqn_9_5__w_h(V: float, r_h: float, t_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        w_h = V*r_h/t_h
        result.append(w_h)
        return result

    @staticmethod
    def eqn_9_5__V(r_h: float, t_h: float, w_h: float):
        # [.pyeqn] w_h = r_h * V / t_h
        result = []
        V = t_h*w_h/r_h
        result.append(V)
        return result


class LiquidRing:

    @staticmethod
    def eqn_10_1__D_r(sig_R: float, w: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        D_r = 229.357798165138*sig_R/w
        result.append(D_r)
        return result

    @staticmethod
    def eqn_10_1__w(D_r: float, sig_R: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        w = 229.357798165138*sig_R/D_r
        result.append(w)
        return result

    @staticmethod
    def eqn_10_1__sig_R(D_r: float, w: float):
        # [.pyeqn] sig_R = 0.00436 * D_r * w
        result = []
        sig_R = 0.00436*D_r*w
        result.append(sig_R)
        return result

    @staticmethod
    def eqn_10_2__PS(Q_gas: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        PS = Q_gas - V*dP/dt
        result.append(PS)
        return result

    @staticmethod
    def eqn_10_2__dP(PS: float, Q_gas: float, V: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dP = dt*(-PS + Q_gas)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_10_2__dt(PS: float, Q_gas: float, V: float, dP: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        dt = -V*dP/(PS - Q_gas)
        result.append(dt)
        return result

    @staticmethod
    def eqn_10_2__V(PS: float, Q_gas: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        V = dt*(-PS + Q_gas)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_10_2__Q_gas(PS: float, V: float, dP: float, dt: float):
        # [.pyeqn] PS = - V * dP / dt + Q_gas
        result = []
        Q_gas = PS + V*dP/dt
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_3__T(N_mfw: float, Q_gas: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        T = 0.108108108108108*Q_gas/N_mfw
        result.append(T)
        return result

    @staticmethod
    def eqn_10_3__N_mfw(Q_gas: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        N_mfw = 0.108108108108108*Q_gas/T
        result.append(N_mfw)
        return result

    @staticmethod
    def eqn_10_3__Q_gas(N_mfw: float, T: float):
        # [.pyeqn] Q_gas = 9.25 * N_mfw * T
        result = []
        Q_gas = 9.25*N_mfw*T
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_4__SP_2(Q_gas: float, SP_1: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_2 = (Q_gas*exp(S_p*t/V) - Q_gas + SP_1)*exp(-S_p*t/V)
        result.append(SP_2)
        return result

    @staticmethod
    def eqn_10_4__S_p(Q_gas: float, SP_1: float, SP_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        S_p = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_4__t(Q_gas: float, SP_1: float, SP_2: float, S_p: float, V: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        t = V*log((Q_gas - SP_1)/(Q_gas - SP_2))/S_p
        result.append(t)
        return result

    @staticmethod
    def eqn_10_4__SP_1(Q_gas: float, SP_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        SP_1 = Q_gas + (-Q_gas + SP_2)*exp(S_p*t/V)
        result.append(SP_1)
        return result

    @staticmethod
    def eqn_10_4__V(Q_gas: float, SP_1: float, SP_2: float, S_p: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        V = S_p*t/log((Q_gas - SP_1)/(Q_gas - SP_2))
        result.append(V)
        return result

    @staticmethod
    def eqn_10_4__Q_gas(SP_1: float, SP_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log((SP_1 - Q_gas) / (SP_2 - Q_gas))
        result = []
        Q_gas = -(SP_1 - SP_2*exp(S_p*t/V))/(exp(S_p*t/V) - 1)
        result.append(Q_gas)
        return result

    @staticmethod
    def eqn_10_5__S_p(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        S_p = V*log(P_1/P_2)/t
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_5__t(P_1: float, P_2: float, S_p: float, V: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_p
        result.append(t)
        return result

    @staticmethod
    def eqn_10_5__P_1(P_2: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_p*t/V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_5__P_2(P_1: float, S_p: float, V: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_p*t/V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_5__V(P_1: float, P_2: float, S_p: float, t: float):
        # [.pyeqn] t = V / S_p * log(P_1 / P_2)
        result = []
        V = S_p*t/log(P_1/P_2)
        result.append(V)
        return result

    @staticmethod
    def eqn_10_6__t(P_1: float, P_2: float, S_a: float, V: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        t = V*log(P_1/P_2)/S_a
        result.append(t)
        return result

    @staticmethod
    def eqn_10_6__P_1(P_2: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_1 = P_2*exp(S_a*t/V)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_10_6__S_a(P_1: float, P_2: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        S_a = V*log(P_1/P_2)/t
        result.append(S_a)
        return result

    @staticmethod
    def eqn_10_6__P_2(P_1: float, S_a: float, V: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        P_2 = P_1*exp(-S_a*t/V)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_10_6__V(P_1: float, P_2: float, S_a: float, t: float):
        # [.pyeqn] S_a = V / t * log(P_1 / P_2)
        result = []
        V = S_a*t/log(P_1/P_2)
        result.append(V)
        return result

    @staticmethod
    def eqn_10_8__c_p(bhp: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        c_p = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(delta_T*f_a*rho)
        result.append(c_p)
        return result

    @staticmethod
    def eqn_10_8__f_a(bhp: float, c_p: float, delta_T: float, delta_h_i: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        f_a = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*rho)
        result.append(f_a)
        return result

    @staticmethod
    def eqn_10_8__delta_h_i(bhp: float, c_p: float, delta_T: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_h_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/w_i
        result.append(delta_h_i)
        return result

    @staticmethod
    def eqn_10_8__rho(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        rho = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*delta_T*f_a)
        result.append(rho)
        return result

    @staticmethod
    def eqn_10_8__delta_T(bhp: float, c_p: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        delta_T = 0.124688279301746*(2545.0*bhp + delta_h_i*w_i)/(c_p*f_a*rho)
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_10_8__w_i(bhp: float, c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        w_i = 0.02*(-127250.0*bhp + 401.0*c_p*delta_T*f_a*rho)/delta_h_i
        result.append(w_i)
        return result

    @staticmethod
    def eqn_10_8__bhp(c_p: float, delta_T: float, delta_h_i: float, f_a: float, rho: float, w_i: float):
        # [.pyeqn] delta_T = (2545 * bhp + w_i * delta_h_i) / ( 8.02 * f_a * rho * c_p )
        result = []
        bhp = 0.00315127701375246*c_p*delta_T*f_a*rho - 0.000392927308447937*delta_h_i*w_i
        result.append(bhp)
        return result

    @staticmethod
    def eqn_10_9__delta_T(T_c: float, T_s: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        delta_T = T_c - T_s
        result.append(delta_T)
        return result

    @staticmethod
    def eqn_10_9__T_s(T_c: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_s = T_c - delta_T
        result.append(T_s)
        return result

    @staticmethod
    def eqn_10_9__T_c(T_s: float, delta_T: float):
        # [.pyeqn] T_c = T_s + delta_T
        result = []
        T_c = T_s + delta_T
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_10__rho(bhp: float, bhp_0: float, mu: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        # >>> from sympy import symbols as s
        # >>> bhp_0,rho,mu,bhp=[s('bhp_0'),s('rho'),s('mu'),s('bhp')]
        # >>> eee = sympy.Eq(0,  bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16) - bhp)
        # >>> sympy.solve(eee, rho) 
        # [FAILED TO PARSE VIA SYMPY]
        # obv, soln is ((((bhp / bhp_0) - 0.5) / 0.0155) / mu ** 0.16) ** (1/0.84))
        # [clearly, LLM is needed...]
        # 1. solve ;
        # 2. py instruct
        # 3. pipe-in python code and de-flower imports, library calls in discord
        rho = symbols('rho')
        equation = Eq(0, bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16) - bhp)
        closed_form_solution = solve(equation, rho)
        return closed_form_solution[0] if closed_form_solution else None

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
    def eqn_10_10__bhp_0(bhp: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp_0 = 2000.0*bhp/(31.0*mu**0.16*rho**0.84 + 1000.0)
        result.append(bhp_0)
        return result

    @staticmethod
    def eqn_10_10__bhp(bhp_0: float, mu: float, rho: float):
        # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
        result = []
        bhp = 0.0005*bhp_0*(31.0*mu**(4/25)*rho**(21/25) + 1000.0)
        result.append(bhp)
        return result

    @staticmethod
    def eqn_10_11__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_s = T_c - 10
        result.append(T_s)
        return result

    @staticmethod
    def eqn_10_11__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 10
        result = []
        T_c = T_s + 10
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_12__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_s = T_c - 5
        result.append(T_s)
        return result

    @staticmethod
    def eqn_10_12__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 5
        result = []
        T_c = T_s + 5
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_13__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_s = T_c - 25
        result.append(T_s)
        return result

    @staticmethod
    def eqn_10_13__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 25
        result = []
        T_c = T_s + 25
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_14__T_s(T_c: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_s = T_c - 12
        result.append(T_s)
        return result

    @staticmethod
    def eqn_10_14__T_c(T_s: float):
        # [.pyeqn] T_c = T_s + 12
        result = []
        T_c = T_s + 12
        result.append(T_c)
        return result

    @staticmethod
    def eqn_10_15__P(S_Th: float, S_p: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        P = S_Th*p_s/(S_Th - S_p)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_15__p_s(P: float, S_Th: float, S_p: float):
        # [.pyeqn] S_p = S_Th * (P - p_s) / P
        result = []
        p_s = P*(S_Th - S_p)/S_Th
        result.append(p_s)
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
    def eqn_10_17__S_Th(P: float, S_0: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_Th = S_0*((P - p_s)/(P - p_0))**(3/5)
        result.append(S_Th)
        return result

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

    @staticmethod
    def eqn_10_17__S_0(P: float, S_Th: float, p_0: float, p_s: float):
        # [.pyeqn] S_Th = S_0 * ((P-p_s) / (P - p_0)) ** 0.6
        result = []
        S_0 = S_Th/((P - p_s)/(P - p_0))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_18__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_p = S_Th*(P*T_i + 460*P - T_i*p_s - 460*p_s)/(P*T_e + 460*P - T_e*p_c - 460*p_c)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_18__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        S_Th = S_p*(P*T_e + 460*P - T_e*p_c - 460*p_c)/(P*T_i + 460*P - T_i*p_s - 460*p_s)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_18__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_e = (P*S_Th*T_i + 460*P*S_Th - 460*P*S_p - S_Th*T_i*p_s - 460*S_Th*p_s + 460*S_p*p_c)/(S_p*(P - p_c))
        result.append(T_e)
        return result

    @staticmethod
    def eqn_10_18__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        P = (S_Th*T_i*p_s + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*T_i + 460*S_Th - S_p*T_e - 460*S_p)
        result.append(P)
        return result

    @staticmethod
    def eqn_10_18__p_s(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_s = (P*S_Th*T_i + 460*P*S_Th - P*S_p*T_e - 460*P*S_p + S_p*T_e*p_c + 460*S_p*p_c)/(S_Th*(T_i + 460))
        result.append(p_s)
        return result

    @staticmethod
    def eqn_10_18__p_c(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        p_c = (-P*S_Th*T_i - 460*P*S_Th + P*S_p*T_e + 460*P*S_p + S_Th*T_i*p_s + 460*S_Th*p_s)/(S_p*(T_e + 460))
        result.append(p_c)
        return result

    @staticmethod
    def eqn_10_18__T_i(P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * (P - p_s)*(460 + T_i) / ((P - p_c)*(460 + T_e) )
        result = []
        T_i = (-460*P*S_Th + P*S_p*T_e + 460*P*S_p + 460*S_Th*p_s - S_p*T_e*p_c - 460*S_p*p_c)/(S_Th*(P - p_s))
        result.append(T_i)
        return result

    @staticmethod
    def eqn_10_19__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_p = S_Th*((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_19__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_Th = S_p/((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_Th)
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
    def eqn_10_19__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        P = (T_e*p_c*(S_p/S_Th)**(5/3) - T_i*p_s + 460.0*p_c*(S_p/S_Th)**(5/3) - 460.0*p_s)/(T_e*(S_p/S_Th)**1.66666666666667 - T_i + 460.0*(S_p/S_Th)**1.66666666666667 - 460.0)
        result.append(P)
        P = (0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/(0.487139289628747*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - T_i + 224.084073229223*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0)
        result.append(P)
        P = (0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/(0.487139289628747*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - T_i + 224.084073229223*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0)
        result.append(P)
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
    def eqn_10_20__S_p(P: float, S_0: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_p = S_0/((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_20__p_0(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        p_0 = P*(-P*T_e*(S_0/S_p)**(5/3) + P*T_i - 460.0*P*(S_0/S_p)**(5/3) + 460.0*P + T_e*p_s*(S_0/S_p)**(5/3) - T_i*p_c - 460.0*p_c + 460.0*p_s*(S_0/S_p)**(5/3))/(P*T_i + 460.0*P - T_i*p_c - 460.0*p_c)
        result.append(p_0)
        p_0 = P*(-0.487139289628747*P*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 - T_i*p_c - 460.0*p_c + 224.084073229223*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5)/(P*T_i + 460.0*P - T_i*p_c - 460.0*p_c)
        result.append(p_0)
        p_0 = P*(-0.487139289628747*P*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 - T_i*p_c - 460.0*p_c + 224.084073229223*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5)/(P*T_i + 460.0*P - T_i*p_c - 460.0*p_c)
        result.append(p_0)
        return result

    @staticmethod
    def eqn_10_20__T_e(P: float, S_0: float, S_p: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        T_e = (P**2*T_i - 460.0*P**2*(S_0/S_p)**(5/3) + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + 460.0*P*p_s*(S_0/S_p)**(5/3) + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(S_0/S_p)**(5/3)*(P - p_s))
        result.append(T_e)
        T_e = 2.05280095711867*(P**2*T_i - 224.084073229223*P**2*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + 224.084073229223*P*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P - p_s)*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5)
        result.append(T_e)
        T_e = 2.05280095711867*(P**2*T_i - 224.084073229223*P**2*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + 224.084073229223*P*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P - p_s)*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5)
        result.append(T_e)
        return result

    @staticmethod
    def eqn_10_20__P(S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        # [FAILED TO PARSE VIA SYMPY]
        """
        Solve the given equation for P.
        :param S_0: Float (base sensitivity level).
        :param S_p: Float (price sensitivity level).
        :param T_e: Float (external temperature in Kelvin).
        :param T_i: Float (internal temperature in Kelvin).
        :param p_0: Float (reference price level).
        :param p_c: Float (commodity price level).
        :param p_s: Float (special price level).
        """
        return [ S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / ]
        (P * (P - p_s)*(460 + T_e))**0.6 - S_0)
        P = 1.0
        tolerance = 1e-6
        while True:
        f_new = RHS(P)
        df_dp = lambda P: -S_p * ((460 + T_i) / (P * (P - p_s)*(460 + T_e))**0.6 - 1/6)
        delta_P = f_new / df_dp(P)
        if abs(delta_P) < tolerance:
        return [ P]
        else:
        P -= delta_P
        P_solution = find_P()
        return [ P_solution]

    @staticmethod
    def eqn_10_20__p_s(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        p_s = (P**2*T_e*(S_0/S_p)**(5/3) - P**2*T_i + 460.0*P**2*(S_0/S_p)**(5/3) - 460.0*P**2 + P*T_i*p_0 + P*T_i*p_c + 460.0*P*p_0 + 460.0*P*p_c - T_i*p_0*p_c - 460.0*p_0*p_c)/(P*(S_0/S_p)**(5/3)*(T_e + 460.0))
        result.append(p_s)
        p_s = 2.05280095711867*(0.487139289628747*P**2*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 - P**2*T_i + 224.084073229223*P**2*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 - 460.0*P**2 + P*T_i*p_0 + P*T_i*p_c + 460.0*P*p_0 + 460.0*P*p_c - T_i*p_0*p_c - 460.0*p_0*p_c)/(P*(T_e + 460.0)*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5)
        result.append(p_s)
        p_s = 2.05280095711867*(0.487139289628747*P**2*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 - P**2*T_i + 224.084073229223*P**2*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 - 460.0*P**2 + P*T_i*p_0 + P*T_i*p_c + 460.0*P*p_0 + 460.0*P*p_c - T_i*p_0*p_c - 460.0*p_0*p_c)/(P*(T_e + 460.0)*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5)
        result.append(p_s)
        return result

    @staticmethod
    def eqn_10_20__p_c(P: float, S_0: float, S_p: float, T_e: float, T_i: float, p_0: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        p_c = P*(-P*T_e*(S_0/S_p)**(5/3) + P*T_i - 460.0*P*(S_0/S_p)**(5/3) + 460.0*P + T_e*p_s*(S_0/S_p)**(5/3) - T_i*p_0 - 460.0*p_0 + 460.0*p_s*(S_0/S_p)**(5/3))/(P*T_i + 460.0*P - T_i*p_0 - 460.0*p_0)
        result.append(p_c)
        p_c = P*(-0.487139289628747*P*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 - T_i*p_0 - 460.0*p_0 + 224.084073229223*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5)/(P*T_i + 460.0*P - T_i*p_0 - 460.0*p_0)
        result.append(p_c)
        p_c = P*(-0.487139289628747*P*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 - T_i*p_0 - 460.0*p_0 + 224.084073229223*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5)/(P*T_i + 460.0*P - T_i*p_0 - 460.0*p_0)
        result.append(p_c)
        return result

    @staticmethod
    def eqn_10_20__S_0(P: float, S_p: float, T_e: float, T_i: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        S_0 = S_p*((P**2*T_i + 460.0*P**2 - P*T_i*p_0 - P*T_i*p_c - 460.0*P*p_0 - 460.0*P*p_c + T_i*p_0*p_c + 460.0*p_0*p_c)/(P*(P*T_e + 460.0*P - T_e*p_s - 460.0*p_s)))**(3/5)
        result.append(S_0)
        return result

    @staticmethod
    def eqn_10_20__T_i(P: float, S_0: float, S_p: float, T_e: float, p_0: float, p_c: float, p_s: float):
        # [.pyeqn] S_0 = S_p * ((P - p_0)*(460 + T_i) * (P - p_c) / (P * (P - p_s)*(460 + T_e) ) )**0.6
        result = []
        T_i = (P**2*T_e*(S_0/S_p)**(5/3) + 460.0*P**2*(S_0/S_p)**(5/3) - 460.0*P**2 - P*T_e*p_s*(S_0/S_p)**(5/3) + 460.0*P*p_0 + 460.0*P*p_c - 460.0*P*p_s*(S_0/S_p)**(5/3) - 460.0*p_0*p_c)/(P**2 - P*p_0 - P*p_c + p_0*p_c)
        result.append(T_i)
        T_i = (0.487139289628747*P**2*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + 224.084073229223*P**2*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 - 460.0*P**2 - 0.487139289628747*P*T_e*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P*p_0 + 460.0*P*p_c - 224.084073229223*P*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 - I*(S_0/S_p)**0.333333333333333)**5 - 460.0*p_0*p_c)/(P**2 - P*p_0 - P*p_c + p_0*p_c)
        result.append(T_i)
        T_i = (0.487139289628747*P**2*T_e*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + 224.084073229223*P**2*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 - 460.0*P**2 - 0.487139289628747*P*T_e*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 + 460.0*P*p_0 + 460.0*P*p_c - 224.084073229223*P*p_s*(-0.577350269189626*(S_0/S_p)**0.333333333333333 + I*(S_0/S_p)**0.333333333333333)**5 - 460.0*p_0*p_c)/(P**2 - P*p_0 - P*p_c + p_0*p_c)
        result.append(T_i)
        return result

    @staticmethod
    def eqn_10_21__P(P_d: float, P_prime: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P = P_d*P_prime/760
        result.append(P)
        return result

    @staticmethod
    def eqn_10_21__P_prime(P: float, P_d: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_prime = 760*P/P_d
        result.append(P_prime)
        return result

    @staticmethod
    def eqn_10_21__P_d(P: float, P_prime: float):
        # [.pyeqn] P_prime = P / P_d * 760
        result = []
        P_d = 760*P/P_prime
        result.append(P_d)
        return result


class RotaryPistonVane:

    @staticmethod
    def eqn_11_1__PS(Q_0: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        PS = Q_0 + Q_external_gas_throughput - V*dP/dT
        result.append(PS)
        return result

    @staticmethod
    def eqn_11_1__Q_external_gas_throughput(PS: float, Q_0: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_external_gas_throughput = PS - Q_0 + V*dP/dT
        result.append(Q_external_gas_throughput)
        return result

    @staticmethod
    def eqn_11_1__dP(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dP = dT*(-PS + Q_0 + Q_external_gas_throughput)/V
        result.append(dP)
        return result

    @staticmethod
    def eqn_11_1__Q_0(PS: float, Q_external_gas_throughput: float, V: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        Q_0 = PS - Q_external_gas_throughput + V*dP/dT
        result.append(Q_0)
        return result

    @staticmethod
    def eqn_11_1__V(PS: float, Q_0: float, Q_external_gas_throughput: float, dP: float, dT: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        V = dT*(-PS + Q_0 + Q_external_gas_throughput)/dP
        result.append(V)
        return result

    @staticmethod
    def eqn_11_1__dT(PS: float, Q_0: float, Q_external_gas_throughput: float, V: float, dP: float):
        # [.pyeqn] PS = -V * dP / dT + Q_external_gas_throughput + Q_0
        result = []
        dT = V*dP/(-PS + Q_0 + Q_external_gas_throughput)
        result.append(dT)
        return result

    @staticmethod
    def eqn_11_2__t(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        left_side = 0 - (V / S_vol_pump_speed * log((SP_1 - (Q + Q_0))/ (SP_2 - (Q + Q_0)))) - (t)
        return [ left_side - t]
        tolerance = 1e-6
        max_iterations = 1000
        initial_guess = 1.0
        for _ in range(max_iterations):
        t_new = f(initial_guess)
        if abs(t_new) < tolerance:
        return [ t_new]
        initial_guess -= t_new / derivative_f(initial_guess)
        raise ValueError("Failed to converge after {} iterations".format(max_iterations))
        return [ V / S_vol_pump_speed * log((SP_1 - (Q + Q_0))/ (SP_2 - (Q + Q_0)))]

    @staticmethod
    def eqn_11_2__S_vol_pump_speed(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        return [ (V / S_vol_pump_speed) * log((SP_1 - (Q_external_gas_throughput + Q_0))/(SP_2 - (Q + Q_0))) - t]
        S_vol_pump_speed = symbols('S_vol_pump_speed')
        equation = eqn_11_2__S_vol_pump_speed(Q=0, Q_0=0, Q_external_gas_throughput=0, SP_1=0, SP_2=0, V=0, t=0) - 0
        solution = solve(equation, S_vol_pump_speed)
        print(f"Closed-form expression for S_vol_pump_speed: {solution}")

    @staticmethod
    def eqn_11_2__Q_0(Q: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        Q_0 = V / (S_vol_pump_speed * ln((SP_1 - Q_external_gas_throughput - SP_2 + Q) / (SP_1 - SP_2))) - t
        return [ Q_0]

    @staticmethod
    def eqn_11_2__SP_1(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        return [ 0 - (V / S_vol_pump_speed * log((SP_1 - (Q + Q_0))/(SP_2 - (Q + Q_0))) - t)]
        SP_1_solution = fsolve(equation, SP_1_initial_guess)[0]
        return [ SP_1_solution]

    @staticmethod
    def eqn_11_2__Q_external_gas_throughput(Q: float, Q_0: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        Q_ext_gas = symbols('Q_ext_gas')
        rearranged_eq = V / S_vol_pump_speed * log((SP_1 - (Q + Q_0)) / (SP_2 - (Q_ext_gas + Q_0))) - t
        closed_form = solve(rearranged_eq, Q_ext_gas)

    @staticmethod
    def eqn_11_2__V(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        return [ V = (SP_1 - (Q + Q_0)) / (SP_2 - (Q + Q_0)) * ln((SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - S_vol_pump_speed*t))]

    @staticmethod
    def eqn_11_2__Q(Q_0: float, Q_external_gas_throughput: float, SP_1: float, SP_2: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        lhs = (V / S_vol_pump_speed * log((SP_1 - Q_external_gas_throughput - Q_0) / (SP_2 - (Q + Q_0)))) - t
        return [exp(lhs)]

    @staticmethod
    def eqn_11_2__SP_2(Q: float, Q_0: float, Q_external_gas_throughput: float, SP_1: float, S_vol_pump_speed: float, V: float, t: float):
        # [.pyeqn] t = V / S_vol_pump_speed * ln( (SP_1 - (Q_external_gas_throughput + Q_0))/ (SP_2 - (Q + Q_0)))
        # [FAILED TO PARSE VIA SYMPY]
        lhs = V / S_vol_pump_speed * log((SP_1 - (Q + Q_0))) - t
        return [ SP_2 = SP_1 - ((V * lhs) / (S_vol_pump_speed * log(SP_1 - (Q_external_gas_throughput + Q_0))))]

    @staticmethod
    def eqn_11_3__t_c(F_s: float, t: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t_c = t/F_s
        result.append(t_c)
        return result

    @staticmethod
    def eqn_11_3__F_s(t: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        F_s = t/t_c
        result.append(F_s)
        return result

    @staticmethod
    def eqn_11_3__t(F_s: float, t_c: float):
        # [.pyeqn] t = t_c * F_s
        result = []
        t = F_s*t_c
        result.append(t)
        return result

    @staticmethod
    def eqn_11_4__p_g(p_s: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_g = p_s - p_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_4__p_v(p_g: float, p_s: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_v = 0
        result.append(p_v)
        p_v = -p_g + p_s
        result.append(p_v)
        return result

    @staticmethod
    def eqn_11_4__p_s(p_g: float, p_v: float):
        # [.pyeqn] p_v / (p_v + p_g) = p_v / p_s
        result = []
        p_s = p_g + p_v
        result.append(p_s)
        return result

    @staticmethod
    def eqn_11_5__p_v_max(P_0_v: float, P_D: float, p_g: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_v_max = -P_0_v*p_g/(P_0_v - P_D)
        result.append(p_v_max)
        return result

    @staticmethod
    def eqn_11_5__p_g(P_0_v: float, P_D: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        p_g = p_v_max*(-P_0_v + P_D)/P_0_v
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_5__P_0_v(P_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_0_v = P_D*p_v_max/(p_g + p_v_max)
        result.append(P_0_v)
        return result

    @staticmethod
    def eqn_11_5__P_D(P_0_v: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = P_0_v * p_g / (P_D - P_0_v)
        result = []
        P_D = P_0_v*(p_g + p_v_max)/p_v_max
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_6__S_B(P_0_V: float, P_D: float, P_v_0: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = S_D*(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)/(P_D*(P_0_V - p_b))
        result.append(S_B)
        return result

    @staticmethod
    def eqn_11_6__S_D(P_0_V: float, P_D: float, P_v_0: float, S_B: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = P_D*S_B*(P_0_V - p_b)/(P_D*p_v_max - P_v_0*p_g - P_v_0*p_v_max)
        result.append(S_D)
        return result

    @staticmethod
    def eqn_11_6__p_g(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (-P_0_V*P_D*S_B + P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_v_max)/(P_v_0*S_D)
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_6__P_v_0(P_0_V: float, P_D: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = P_D*(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)/(S_D*(p_g + p_v_max))
        result.append(P_v_0)
        return result

    @staticmethod
    def eqn_11_6__p_b(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (P_0_V*P_D*S_B - P_D*S_D*p_v_max + P_v_0*S_D*p_g + P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(p_b)
        return result

    @staticmethod
    def eqn_11_6__P_D(P_0_V: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0*S_D*(p_g + p_v_max)/(-P_0_V*S_B + S_B*p_b + S_D*p_v_max)
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_6__P_0_V(P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float, p_v_max: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (P_D*S_B*p_b + P_D*S_D*p_v_max - P_v_0*S_D*p_g - P_v_0*S_D*p_v_max)/(P_D*S_B)
        result.append(P_0_V)
        return result

    @staticmethod
    def eqn_11_6__p_v_max(P_0_V: float, P_D: float, P_v_0: float, S_B: float, S_D: float, p_b: float, p_g: float):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V*P_D*S_B - P_D*S_B*p_b + P_v_0*S_D*p_g)/(S_D*(P_D - P_v_0))
        result.append(p_v_max)
        return result
