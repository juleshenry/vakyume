from math import log, sqrt, exp


class VacuumTheory:

    @staticmethod
    def eqn_1_3__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return T

    @staticmethod
    def eqn_1_3__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return k

    @staticmethod
    def eqn_1_3__v(T: float, k: float, m: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        v = -1.73205080756888*sqrt(T*k/m)
        result.append(v)
        v = 1.73205080756888*sqrt(T*k/m)
        result.append(v)
        return v

    @staticmethod
    def eqn_1_3__m(T: float, k: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        m = 3.0*T*k/v**2
        result.append(m)
        return m

    @staticmethod
    def eqn_1_7__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return p

    @staticmethod
    def eqn_1_7__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return V

    @staticmethod
    def eqn_1_7__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return n

    @staticmethod
    def eqn_1_7__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return R

    @staticmethod
    def eqn_1_7__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return T

    @staticmethod
    def eqn_1_8__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return V

    @staticmethod
    def eqn_1_8__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return P

    @staticmethod
    def eqn_1_8__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return R

    @staticmethod
    def eqn_1_8__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return m

    @staticmethod
    def eqn_1_8__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return T

    @staticmethod
    def eqn_1_8__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return M

    @staticmethod
    def eqn_1_9__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return P

    @staticmethod
    def eqn_1_9__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return R

    @staticmethod
    def eqn_1_9__T(M: float, P: float, R: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return T

    @staticmethod
    def eqn_1_9__rho(M: float, P: float, R: float, T: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_1_9__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return M

    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return P_2

    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return V_1

    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return V_2

    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return P_1

    @staticmethod
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return T_2

    @staticmethod
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return T_1

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return q

    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return P

    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return W

    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return T

    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60) ft^3/min
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return M

    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return Total_P

    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return sum_partial_pressures

    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return n_a

    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return n

    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return y_a

    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return P

    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return y_a

    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return p_a


class FluidFlowVacuumLines:

    @staticmethod
    def eqn_2_1__D(Re: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        D = Re*mu/(rho*v)
        result.append(D)
        return D

    @staticmethod
    def eqn_2_1__Re(D: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        Re = D*rho*v/mu
        result.append(Re)
        return Re

    @staticmethod
    def eqn_2_1__v(D: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        v = Re*mu/(D*rho)
        result.append(v)
        return v

    @staticmethod
    def eqn_2_1__rho(D: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        rho = Re*mu/(D*v)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_2_1__mu(D: float, Re: float, rho: float, v: float):
        # [.pyeqn] Re = rho * D * v / mu
        result = []
        mu = D*rho*v/Re
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_2__delta(lambd: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        delta = -0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        delta = 0.474424998328794*sqrt(lambd/psi)
        result.append(delta)
        return delta

    @staticmethod
    def eqn_2_2__lambd(delta: float, psi: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        lambd = 4.44288293815837*delta**2*psi
        result.append(lambd)
        return lambd

    @staticmethod
    def eqn_2_2__psi(delta: float, lambd: float):
        # [.pyeqn] lambd = 3.141592653589793 * delta ** 2 * psi * 2 ** 0.5
        result = []
        psi = 0.225079079039277*lambd/delta**2
        result.append(psi)
        return psi

    @staticmethod
    def eqn_2_3__D(kn: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        D = lambd/kn
        result.append(D)
        return D

    @staticmethod
    def eqn_2_3__lambd(D: float, kn: float):
        # [.pyeqn] kn = lambd / D
        result = []
        lambd = D*kn
        result.append(lambd)
        return lambd

    @staticmethod
    def eqn_2_3__kn(D: float, lambd: float):
        # [.pyeqn] kn = lambd / D
        result = []
        kn = lambd/D
        result.append(kn)
        return kn

    @staticmethod
    def eqn_2_4___beta(mu: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        _beta = mu*vel_grad
        result.append(_beta)
        return _beta

    @staticmethod
    def eqn_2_4__vel_grad(_beta: float, mu: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        vel_grad = _beta/mu
        result.append(vel_grad)
        return vel_grad

    @staticmethod
    def eqn_2_4__mu(_beta: float, vel_grad: float):
        # [.pyeqn] _beta = mu * vel_grad
        result = []
        mu = _beta/vel_grad
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_5__D(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        D = -2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = -2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        D = 2.52647511098426*(L*mu*q/delta_P)**(1/4)
        result.append(D)
        return D

    @staticmethod
    def eqn_2_5__q(D: float, L: float, delta_P: float, mu: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        q = 0.0245436926061703*D**4*delta_P/(L*mu)
        result.append(q)
        return q

    @staticmethod
    def eqn_2_5__L(D: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        L = 0.0245436926061703*D**4*delta_P/(mu*q)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_5__delta_P(D: float, L: float, mu: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        delta_P = 40.7436654315252*L*mu*q/D**4
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_2_5__mu(D: float, L: float, delta_P: float, q: float):
        # [.pyeqn] q = 3.141592653589793 * (D ** 4) * delta_P / (128 * L * mu)
        result = []
        mu = 0.0245436926061703*D**4*delta_P/(L*q)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_6__rho(lambd: float, mu: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        rho = 2.85714285714286*mu/(lambd*v_a)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_2_6__v_a(lambd: float, mu: float, rho: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        v_a = 2.85714285714286*mu/(lambd*rho)
        result.append(v_a)
        return v_a

    @staticmethod
    def eqn_2_6__lambd(mu: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        lambd = 2.85714285714286*mu/(rho*v_a)
        result.append(lambd)
        return lambd

    @staticmethod
    def eqn_2_6__mu(lambd: float, rho: float, v_a: float):
        # [.pyeqn] mu = 0.35 * rho * lambd * v_a
        result = []
        mu = 0.35*lambd*rho*v_a
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_7__T(k: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        T = 0.392699081698724*m*v_a**2/k
        result.append(T)
        return T

    @staticmethod
    def eqn_2_7__k(T: float, m: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        k = 0.392699081698724*m*v_a**2/T
        result.append(k)
        return k

    @staticmethod
    def eqn_2_7__v_a(T: float, k: float, m: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        v_a = 1.59576912160573*sqrt(T*k/m)
        result.append(v_a)
        return v_a

    @staticmethod
    def eqn_2_7__m(T: float, k: float, v_a: float):
        # [.pyeqn] v_a = ((8 * k * T) / (3.141592653589793 * m)) ** 0.5
        result = []
        m = 2.54647908947033*T*k/v_a**2
        result.append(m)
        return m

    @staticmethod
    def eqn_2_8__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return mu_c

    @staticmethod
    def eqn_2_8__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return P_c

    @staticmethod
    def eqn_2_8__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return T_c

    @staticmethod
    def eqn_2_8__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return M

    @staticmethod
    def eqn_2_8__mu_c(M: float, P_c: float, T_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        mu_c = 7.7*sqrt(M)*P_c**(2/3)/T_c**(1/6)
        result.append(mu_c)
        return mu_c

    @staticmethod
    def eqn_2_8__P_c(M: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        P_c = -0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        P_c = 0.046801946114055*(T_c**0.166666666666667*mu_c/M**0.5)**(3/2)
        result.append(P_c)
        return P_c

    @staticmethod
    def eqn_2_8__T_c(M: float, P_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        T_c = 208422.380089*M**3*P_c**4/mu_c**6
        result.append(T_c)
        return T_c

    @staticmethod
    def eqn_2_8__M(P_c: float, T_c: float, mu_c: float):
        # [.pyeqn] mu_c = (7.7 * (M ** 0.5) * P_c ** (2 / 3)) / T_c ** (1 / 6)
        result = []
        M = 0.0168662506324844*T_c**(1/3)*mu_c**2/P_c**(4/3)
        result.append(M)
        return M

    @staticmethod
    def eqn_2_10__delta_P(Suc_Pres: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        delta_P = -Suc_Pres + oper_press
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_2_10__oper_press(Suc_Pres: float, delta_P: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        oper_press = Suc_Pres + delta_P
        result.append(oper_press)
        return oper_press

    @staticmethod
    def eqn_2_10__Suc_Pres(delta_P: float, oper_press: float):
        # [.pyeqn] Suc_Pres = oper_press - delta_P
        result = []
        Suc_Pres = -delta_P + oper_press
        result.append(Suc_Pres)
        return Suc_Pres

    @staticmethod
    def eqn_2_11__D(L: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        D = L*f*v**2/(2*g_c*h_r)
        result.append(D)
        return D

    @staticmethod
    def eqn_2_11__L(D: float, f: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        L = 2*D*g_c*h_r/(f*v**2)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_11__f(D: float, L: float, g_c: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        f = 2*D*g_c*h_r/(L*v**2)
        result.append(f)
        return f

    @staticmethod
    def eqn_2_11__h_r(D: float, L: float, f: float, g_c: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        h_r = L*f*v**2/(2*D*g_c)
        result.append(h_r)
        return h_r

    @staticmethod
    def eqn_2_11__g_c(D: float, L: float, f: float, h_r: float, v: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        g_c = L*f*v**2/(2*D*h_r)
        result.append(g_c)
        return g_c

    @staticmethod
    def eqn_2_11__v(D: float, L: float, f: float, g_c: float, h_r: float):
        # [.pyeqn] h_r = f * L * v ** 2 / (D * 2 * g_c)
        result = []
        v = -sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        v = sqrt(2)*sqrt(D*g_c*h_r/(L*f))
        result.append(v)
        return v

    @staticmethod
    def eqn_2_12__L(d: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        L = 0.464037122969838*d*delta_P*g/(f*rho*v**2)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_12__f(L: float, d: float, delta_P: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        f = 0.464037122969838*d*delta_P*g/(L*rho*v**2)
        result.append(f)
        return f

    @staticmethod
    def eqn_2_12__delta_P(L: float, d: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        delta_P = 2.155*L*f*rho*v**2/(d*g)
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_2_12__v(L: float, d: float, delta_P: float, f: float, g: float, rho: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        v = -0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        v = 0.681202703290172*sqrt(d*delta_P*g/(L*f*rho))
        result.append(v)
        return v

    @staticmethod
    def eqn_2_12__g(L: float, d: float, delta_P: float, f: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        g = 2.155*L*f*rho*v**2/(d*delta_P)
        result.append(g)
        return g

    @staticmethod
    def eqn_2_12__d(L: float, delta_P: float, f: float, g: float, rho: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        d = 2.155*L*f*rho*v**2/(delta_P*g)
        result.append(d)
        return d

    @staticmethod
    def eqn_2_12__rho(L: float, d: float, delta_P: float, f: float, g: float, v: float):
        # [.pyeqn] delta_P = 4.31 * rho * f * L * v ** 2 / (2 * d * g)
        result = []
        rho = 0.464037122969838*d*delta_P*g/(L*f*v**2)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_2_13__L(d: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        L = 0.465116279069767*d**5*delta_P/(f*q**2*rho)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_13__f(L: float, d: float, delta_P: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        f = 0.465116279069767*d**5*delta_P/(L*q**2*rho)
        result.append(f)
        return f

    @staticmethod
    def eqn_2_13__q(L: float, d: float, delta_P: float, f: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        q = -0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        q = 0.681994339470473*sqrt(d**5*delta_P/(L*f*rho))
        result.append(q)
        return q

    @staticmethod
    def eqn_2_13__delta_P(L: float, d: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        delta_P = 2.15*L*f*q**2*rho/d**5
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_2_13__d(L: float, delta_P: float, f: float, q: float, rho: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        d = 1.16543402167043*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) - 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = -0.942855929354115*(L*f*q**2*rho/delta_P)**(1/5) + 0.685024930457783*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) - 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        d = 0.360138918518902*(L*f*q**2*rho/delta_P)**(1/5) + 1.10839362062173*I*(L*f*q**2*rho/delta_P)**(1/5)
        result.append(d)
        return d

    @staticmethod
    def eqn_2_13__rho(L: float, d: float, delta_P: float, f: float, q: float):
        # [.pyeqn] delta_P = 2.15 * rho * f * L * q ** 2 / (d ** 5)
        result = []
        rho = 0.465116279069767*d**5*delta_P/(L*f*q**2)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_2_14__v_s(M: float, R: float, T: float, g_c: float, k: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        v_s = sqrt(R*T*g_c*k/M)
        result.append(v_s)
        return v_s

    @staticmethod
    def eqn_2_14__R(M: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        R = M*v_s**2/(T*g_c*k)
        result.append(R)
        return R

    @staticmethod
    def eqn_2_14__g_c(M: float, R: float, T: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        g_c = M*v_s**2/(R*T*k)
        result.append(g_c)
        return g_c

    @staticmethod
    def eqn_2_14__T(M: float, R: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        T = M*v_s**2/(R*g_c*k)
        result.append(T)
        return T

    @staticmethod
    def eqn_2_14__k(M: float, R: float, T: float, g_c: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        k = M*v_s**2/(R*T*g_c)
        result.append(k)
        return k

    @staticmethod
    def eqn_2_14__M(R: float, T: float, g_c: float, k: float, v_s: float):
        # [.pyeqn] v_s = (k * g_c * R / M * T) ** 0.5
        result = []
        M = R*T*g_c*k/v_s**2
        result.append(M)
        return M

    @staticmethod
    def eqn_2_15__Re(f: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        Re = 0.009971220736/f**4
        result.append(Re)
        return Re

    @staticmethod
    def eqn_2_15__f(Re: float):
        # [.pyeqn] f = 0.316 / Re ** (0.25)
        result = []
        f = 0.316/Re**(1/4)
        result.append(f)
        return f

    @staticmethod
    def eqn_2_16__Re(f: float):
        # [.pyeqn] f = 64 / Re
        result = []
        Re = 64/f
        result.append(Re)
        return Re

    @staticmethod
    def eqn_2_16__f(Re: float):
        # [.pyeqn] f = 64 / Re
        result = []
        f = 64/Re
        result.append(f)
        return f

    @staticmethod
    def eqn_2_16__L(d: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        L = 28.9855072463768*d**2*delta_P/(mu*v)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_16__delta_P(L: float, d: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        delta_P = 0.0345*L*mu*v/d**2
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_2_16__d(L: float, delta_P: float, mu: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        d = -0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        d = 0.185741756210067*sqrt(L*mu*v/delta_P)
        result.append(d)
        return d

    @staticmethod
    def eqn_2_16__v(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        v = 28.9855072463768*d**2*delta_P/(L*mu)
        result.append(v)
        return v

    @staticmethod
    def eqn_2_16__mu(L: float, d: float, delta_P: float, v: float):
        # [.pyeqn] delta_P = 0.0345* mu * L * v / d**2
        result = []
        mu = 28.9855072463768*d**2*delta_P/(L*v)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_16__L(d: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        L = 9.52380952380952*d**4*delta_P/(mu*q)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_16__q(L: float, d: float, delta_P: float, mu: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        q = 9.52380952380952*d**4*delta_P/(L*mu)
        result.append(q)
        return q

    @staticmethod
    def eqn_2_16__delta_P(L: float, d: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        delta_P = 0.105*L*mu*q/d**4
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_2_16__d(L: float, delta_P: float, mu: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        d = -0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*I*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = -0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        d = 0.569242509762222*(L*mu*q/delta_P)**(1/4)
        result.append(d)
        return d

    @staticmethod
    def eqn_2_16__mu(L: float, d: float, delta_P: float, q: float):
        # [.pyeqn] delta_P = 0.105 * mu * L * q / d**4
        result = []
        mu = 9.52380952380952*d**4*delta_P/(L*q)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_18a__R_ll(D_eq: float):
        # [.pyeqn] D_eq = 4 * R_ll
        result = []
        R_ll = D_eq/4
        result.append(R_ll)
        return R_ll

    @staticmethod
    def eqn_2_18a__D_eq(R_ll: float):
        # [.pyeqn] D_eq = 4 * R_ll
        result = []
        D_eq = 4*R_ll
        result.append(D_eq)
        return D_eq

    @staticmethod
    def eqn_2_18b__h(R_ll: float, w: float):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        h = 2*R_ll*w/(-2*R_ll + w)
        result.append(h)
        return h

    @staticmethod
    def eqn_2_18b__w(R_ll: float, h: float):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        w = 2*R_ll*h/(-2*R_ll + h)
        result.append(w)
        return w

    @staticmethod
    def eqn_2_18b__R_ll(h: float, w: float):
        # [.pyeqn] R_ll = w * h / (2 * (w + h))
        result = []
        R_ll = h*w/(2*(h + w))
        result.append(R_ll)
        return R_ll

    @staticmethod
    def eqn_2_19a__Re(R_ll: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        Re = 4*R_ll*rho*v/mu
        result.append(Re)
        return Re

    @staticmethod
    def eqn_2_19a__v(R_ll: float, Re: float, mu: float, rho: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        v = Re*mu/(4*R_ll*rho)
        result.append(v)
        return v

    @staticmethod
    def eqn_2_19a__R_ll(Re: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        R_ll = Re*mu/(4*rho*v)
        result.append(R_ll)
        return R_ll

    @staticmethod
    def eqn_2_19a__rho(R_ll: float, Re: float, mu: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        rho = Re*mu/(4*R_ll*v)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_2_19a__mu(R_ll: float, Re: float, rho: float, v: float):
        # [.pyeqn] Re = 4 * R_ll * rho * v / mu
        result = []
        mu = 4*R_ll*rho*v/Re
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_19b__w(Re: float, h: float, mu: float, rho: float, v: float):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        w = Re*h*mu/(-Re*mu + 2*h*rho*v)
        result.append(w)
        return w

    @staticmethod
    def eqn_2_19b__Re(h: float, mu: float, rho: float, v: float, w: float):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        Re = 2*h*rho*v*w/(mu*(h + w))
        result.append(Re)
        return Re

    @staticmethod
    def eqn_2_19b__v(Re: float, h: float, mu: float, rho: float, w: float):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        v = Re*mu*(h + w)/(2*h*rho*w)
        result.append(v)
        return v

    @staticmethod
    def eqn_2_19b__rho(Re: float, h: float, mu: float, v: float, w: float):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        rho = Re*mu*(h + w)/(2*h*v*w)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_2_19b__h(Re: float, mu: float, rho: float, v: float, w: float):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        h = Re*mu*w/(-Re*mu + 2*rho*v*w)
        result.append(h)
        return h

    @staticmethod
    def eqn_2_19b__mu(Re: float, h: float, rho: float, v: float, w: float):
        # [.pyeqn] Re = (2 * w * h * rho * v) / ((w + h) * mu)
        result = []
        mu = 2*h*rho*v*w/(Re*(h + w))
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_20__sum_pipe(L: float, sum_equivalent_length: float):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        sum_pipe = L - sum_equivalent_length
        result.append(sum_pipe)
        return sum_pipe

    @staticmethod
    def eqn_2_20__L(sum_equivalent_length: float, sum_pipe: float):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        L = sum_equivalent_length + sum_pipe
        result.append(L)
        return L

    @staticmethod
    def eqn_2_20__sum_equivalent_length(L: float, sum_pipe: float):
        # [.pyeqn] L = sum_pipe + sum_equivalent_length
        result = []
        sum_equivalent_length = L - sum_pipe
        result.append(sum_equivalent_length)
        return sum_equivalent_length

    @staticmethod
    def eqn_2_22__P_s(Q_throughput: float, S_p: float):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        P_s = Q_throughput/S_p
        result.append(P_s)
        return P_s

    @staticmethod
    def eqn_2_22__Q_throughput(P_s: float, S_p: float):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        Q_throughput = P_s*S_p
        result.append(Q_throughput)
        return Q_throughput

    @staticmethod
    def eqn_2_22__S_p(P_s: float, Q_throughput: float):
        # [.pyeqn] Q_throughput = S_p * P_s
        result = []
        S_p = Q_throughput/P_s
        result.append(S_p)
        return S_p

    @staticmethod
    def eqn_2_25__P_2(C: float, P_1: float, Q_throughput: float):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        P_2 = P_1 - Q_throughput/C
        result.append(P_2)
        return P_2

    @staticmethod
    def eqn_2_25__Q_throughput(C: float, P_1: float, P_2: float):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        Q_throughput = C*(P_1 - P_2)
        result.append(Q_throughput)
        return Q_throughput

    @staticmethod
    def eqn_2_25__P_1(C: float, P_2: float, Q_throughput: float):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        P_1 = P_2 + Q_throughput/C
        result.append(P_1)
        return P_1

    @staticmethod
    def eqn_2_25__C(P_1: float, P_2: float, Q_throughput: float):
        # [.pyeqn] C = Q_throughput / (P_1 - P_2)
        result = []
        C = Q_throughput/(P_1 - P_2)
        result.append(C)
        return C

    @staticmethod
    def eqn_2_26__D(L: float, P_downstream: float, P_p: float, P_upstream: float, mu: float, q: float):
        # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        D = -14953.4878122122*I*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = 14953.4878122122*I*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = -14953.4878122122*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        D = 14953.4878122122*(-L*mu*q/(1.22718463030851e+15*P_downstream - 1.22718463030851e+15*P_upstream))**(1/4)
        result.append(D)
        return D

    @staticmethod
    def eqn_2_26__q(D: float, L: float, P_downstream: float, P_p: float, P_upstream: float, mu: float):
        # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        q = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(L*mu)
        result.append(q)
        return q

    @staticmethod
    def eqn_2_26__L(D: float, P_downstream: float, P_p: float, P_upstream: float, mu: float, q: float):
        # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        L = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(mu*q)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_26__P_downstream(D: float, L: float, P_p: float, P_upstream: float, mu: float, q: float):
        # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_downstream = P_upstream - 40.7436654315252*L*mu*q/D**4
        result.append(P_downstream)
        return P_downstream

    @staticmethod
    def eqn_2_26__P_upstream(D: float, L: float, P_downstream: float, P_p: float, mu: float, q: float):
        # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_upstream = P_downstream + 40.7436654315252*L*mu*q/D**4
        result.append(P_upstream)
        return P_upstream

    @staticmethod
    def eqn_2_26__P_p(D: float, L: float, P_downstream: float, P_upstream: float, mu: float, q: float):
        # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        P_p = 0.0
        result.append(P_p)
        return P_p

    @staticmethod
    def eqn_2_26__mu(D: float, L: float, P_downstream: float, P_p: float, P_upstream: float, q: float):
        # [.pyeqn] q * P_p = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p * (P_upstream - P_downstream)
        result = []
        mu = 0.0245436926061703*D**4*(-P_downstream + P_upstream)/(L*q)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_28__D(C: float, L: float, P_p: float, mu: float):
        # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        D = -2.52647511098426*I*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = 2.52647511098426*I*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = -2.52647511098426*(C*L*mu/P_p)**(1/4)
        result.append(D)
        D = 2.52647511098426*(C*L*mu/P_p)**(1/4)
        result.append(D)
        return D

    @staticmethod
    def eqn_2_28__L(C: float, D: float, P_p: float, mu: float):
        # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        L = 0.0245436926061703*D**4*P_p/(C*mu)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_28__C(D: float, L: float, P_p: float, mu: float):
        # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        C = 0.0245436926061703*D**4*P_p/(L*mu)
        result.append(C)
        return C

    @staticmethod
    def eqn_2_28__P_p(C: float, D: float, L: float, mu: float):
        # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        P_p = 40.7436654315252*C*L*mu/D**4
        result.append(P_p)
        return P_p

    @staticmethod
    def eqn_2_28__mu(C: float, D: float, L: float, P_p: float):
        # [.pyeqn] C = 3.141592653589793 * D ** 4 / (128 * mu * L) * P_p
        result = []
        mu = 0.0245436926061703*D**4*P_p/(C*L)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_29__S_2(C: float, S_1: float):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_2 = C*S_1/(C - S_1)
        result.append(S_2)
        return S_2

    @staticmethod
    def eqn_2_29__S_1(C: float, S_2: float):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        S_1 = C*S_2/(C + S_2)
        result.append(S_1)
        return S_1

    @staticmethod
    def eqn_2_29__C(S_1: float, S_2: float):
        # [.pyeqn] S_1 ** -1 = S_2 ** -1 + 1 / C
        result = []
        C = -S_1*S_2/(S_1 - S_2)
        result.append(C)
        return C

    @staticmethod
    def eqn_2_31__S_pump_speed(C: float, S_p: float):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_pump_speed = C*S_p/(C + S_p)
        result.append(S_pump_speed)
        return S_pump_speed

    @staticmethod
    def eqn_2_31__S_p(C: float, S_pump_speed: float):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        S_p = C*S_pump_speed/(C - S_pump_speed)
        result.append(S_p)
        return S_p

    @staticmethod
    def eqn_2_31__C(S_p: float, S_pump_speed: float):
        # [.pyeqn] S_pump_speed = (S_p * C) / (S_p + C)
        result = []
        C = S_p*S_pump_speed/(S_p - S_pump_speed)
        result.append(C)
        return C

    @staticmethod
    def eqn_2_32__C_series(geometric_sum_C: float):
        # [.pyeqn] 1 / C_series = geometric_sum_C
        result = []
        C_series = 1/geometric_sum_C
        result.append(C_series)
        return C_series

    @staticmethod
    def eqn_2_32__geometric_sum_C(C_series: float):
        # [.pyeqn] 1 / C_series = geometric_sum_C
        result = []
        geometric_sum_C = 1/C_series
        result.append(geometric_sum_C)
        return geometric_sum_C

    @staticmethod
    def eqn_2_33__C_paralell(arithmetic_sum_C: float):
        # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
        result = []
        C_paralell = 1/arithmetic_sum_C
        result.append(C_paralell)
        return C_paralell

    @staticmethod
    def eqn_2_33__arithmetic_sum_C(C_paralell: float):
        # [.pyeqn] 1 / C_paralell = arithmetic_sum_C
        result = []
        arithmetic_sum_C = 1/C_paralell
        result.append(arithmetic_sum_C)
        return arithmetic_sum_C

    @staticmethod
    def eqn_2_34__L(C: float, C_1: float, C_2: float, D: float, P_p: float, mu: float):
        # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        L = D**3*(C_1*D*P_p + C_2*mu)/(C*mu)
        result.append(L)
        return L

    @staticmethod
    def eqn_2_34__D(C: float, C_1: float, C_2: float, L: float, P_p: float, mu: float):
        # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        D = Piecewise((-sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (-sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((-sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (-sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) + C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 - sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        D = Piecewise((sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*(-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), Eq(C*L*mu/(C_1*P_p), 0)), (sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))/2 + sqrt(2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) - 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(2*C_1**2*P_p**2) - C_2**3*mu**3/(4*C_1**3*P_p**3*sqrt(-2*C*L*mu/(3*C_1*P_p*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3)) + 2*(sqrt(C**3*L**3*mu**3/(27*C_1**3*P_p**3) + (-C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(8*C_1**2*P_p**2) - 3*C_2**6*mu**6/(2048*C_1**6*P_p**6))**2/4) + C_2**2*mu**2*(-C*L*mu/(C_1*P_p) - 3*C_2**4*mu**4/(256*C_1**4*P_p**4))/(16*C_1**2*P_p**2) + 3*C_2**6*mu**6/(4096*C_1**6*P_p**6))**(1/3) + C_2**2*mu**2/(4*C_1**2*P_p**2))))/2 - C_2*mu/(4*C_1*P_p), True))
        result.append(D)
        return D

    @staticmethod
    def eqn_2_34__C_1(C: float, C_2: float, D: float, L: float, P_p: float, mu: float):
        # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C_1 = mu*(C*L - C_2*D**3)/(D**4*P_p)
        result.append(C_1)
        return C_1

    @staticmethod
    def eqn_2_34__C_2(C: float, C_1: float, D: float, L: float, P_p: float, mu: float):
        # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C_2 = C*L/D**3 - C_1*D*P_p/mu
        result.append(C_2)
        return C_2

    @staticmethod
    def eqn_2_34__C(C_1: float, C_2: float, D: float, L: float, P_p: float, mu: float):
        # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        C = D**3*(C_1*D*P_p + C_2*mu)/(L*mu)
        result.append(C)
        return C

    @staticmethod
    def eqn_2_34__P_p(C: float, C_1: float, C_2: float, D: float, L: float, mu: float):
        # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        P_p = mu*(C*L - C_2*D**3)/(C_1*D**4)
        result.append(P_p)
        return P_p

    @staticmethod
    def eqn_2_34__mu(C: float, C_1: float, C_2: float, D: float, L: float, P_p: float):
        # [.pyeqn] C = C_1 * (D ** 4 / (mu * L)) * P_p + C_2 * (D ** 3 / L)
        result = []
        mu = C_1*D**4*P_p/(C*L - C_2*D**3)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_2_11__F_p(C_L: float, C_T: float):
        # [.pyeqn] C_T = C_L * F_p
        result = []
        F_p = C_T/C_L
        result.append(F_p)
        return F_p

    @staticmethod
    def eqn_2_11__C_T(C_L: float, F_p: float):
        # [.pyeqn] C_T = C_L * F_p
        result = []
        C_T = C_L*F_p
        result.append(C_T)
        return C_T

    @staticmethod
    def eqn_2_11__C_L(C_T: float, F_p: float):
        # [.pyeqn] C_T = C_L * F_p
        result = []
        C_L = C_T/F_p
        result.append(C_L)
        return C_L

    @staticmethod
    def eqn_2_36__F_t(C: float, C_0: float):
        # [.pyeqn] C = C_0 * F_t
        result = []
        F_t = C/C_0
        result.append(F_t)
        return F_t

    @staticmethod
    def eqn_2_36__C_0(C: float, F_t: float):
        # [.pyeqn] C = C_0 * F_t
        result = []
        C_0 = C/F_t
        result.append(C_0)
        return C_0

    @staticmethod
    def eqn_2_36__C(C_0: float, F_t: float):
        # [.pyeqn] C = C_0 * F_t
        result = []
        C = C_0*F_t
        result.append(C)
        return C

    @staticmethod
    def eqn_2_37__A(C: float, F_t: float, M: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        A = 0.000681714375311032*C**2*M/(F_t*T)
        result.append(A)
        return A

    @staticmethod
    def eqn_2_37__C(A: float, F_t: float, M: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        C = 38.3*sqrt(A*F_t*T/M)
        result.append(C)
        return C

    @staticmethod
    def eqn_2_37__F_t(A: float, C: float, M: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        F_t = 0.000681714375311032*C**2*M/(A*T)
        result.append(F_t)
        return F_t

    @staticmethod
    def eqn_2_37__T(A: float, C: float, F_t: float, M: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        T = 0.000681714375311032*C**2*M/(A*F_t)
        result.append(T)
        return T

    @staticmethod
    def eqn_2_37__M(A: float, C: float, F_t: float, T: float):
        # [.pyeqn] C = 38.3 * (T * A * F_t / M) ** 0.5
        result = []
        M = 1466.89*A*F_t*T/C**2
        result.append(M)
        return M


class PressMgmt:

    @staticmethod
    def eqn_3_1__Abs_Pressure(BarometricPressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Abs_Pressure = BarometricPressure - Vacuum
        result.append(Abs_Pressure)
        return Abs_Pressure

    @staticmethod
    def eqn_3_1__Vacuum(Abs_Pressure: float, BarometricPressure: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        Vacuum = -Abs_Pressure + BarometricPressure
        result.append(Vacuum)
        return Vacuum

    @staticmethod
    def eqn_3_1__BarometricPressure(Abs_Pressure: float, Vacuum: float):
        # [.pyeqn] Abs_Pressure = BarometricPressure - Vacuum
        result = []
        BarometricPressure = Abs_Pressure + Vacuum
        result.append(BarometricPressure)
        return BarometricPressure

    @staticmethod
    def eqn_3_2__H(G: float, G_C: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        H = G/(G_C*P*rho)
        result.append(H)
        return H

    @staticmethod
    def eqn_3_2__G_C(G: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G_C = G/(H*P*rho)
        result.append(G_C)
        return G_C

    @staticmethod
    def eqn_3_2__P(G: float, G_C: float, H: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        P = G/(G_C*H*rho)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_2__G(G_C: float, H: float, P: float, rho: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        G = G_C*H*P*rho
        result.append(G)
        return G

    @staticmethod
    def eqn_3_2__rho(G: float, G_C: float, H: float, P: float):
        # [.pyeqn] P = G / (G_C * rho * H)
        result = []
        rho = G/(G_C*H*P)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_3_3__H_2(H_1: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_2 = H_1 - P + P_P
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_3__P(H_1: float, H_2: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P = H_1 - H_2 + P_P
        result.append(P)
        return P

    @staticmethod
    def eqn_3_3__H_1(H_2: float, P: float, P_P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        H_1 = H_2 + P - P_P
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_3__P_P(H_1: float, H_2: float, P: float):
        # [.pyeqn] P_P - P = H_2 - H_1
        result = []
        P_P = -H_1 + H_2 + P
        result.append(P_P)
        return P_P

    @staticmethod
    def eqn_3_4__KAPPA(P: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        KAPPA = P*V
        result.append(KAPPA)
        return KAPPA

    @staticmethod
    def eqn_3_4__P(KAPPA: float, V: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        P = KAPPA/V
        result.append(P)
        return P

    @staticmethod
    def eqn_3_4__V(KAPPA: float, P: float):
        # [.pyeqn] P * V = KAPPA
        result = []
        V = KAPPA/P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_5__P(P_P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P = P_P*V_P/V
        result.append(P)
        return P

    @staticmethod
    def eqn_3_5__V_P(P: float, P_P: float, V: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V_P = P*V/P_P
        result.append(V_P)
        return V_P

    @staticmethod
    def eqn_3_5__V(P: float, P_P: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        V = P_P*V_P/P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_5__P_P(P: float, V: float, V_P: float):
        # [.pyeqn] P_P = P * (V / V_P)
        result = []
        P_P = P*V/V_P
        result.append(P_P)
        return P_P

    @staticmethod
    def eqn_3_6__H_2(H_1: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_2 = H_1 + P*V/V_P - P
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_6__V_P(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V_P = P*V/(-H_1 + H_2 + P)
        result.append(V_P)
        return V_P

    @staticmethod
    def eqn_3_6__V(H_1: float, H_2: float, P: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        V = V_P*(-H_1 + H_2 + P)/P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_6__P(H_1: float, H_2: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        P = V_P*(-H_1 + H_2)/(V - V_P)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_6__H_1(H_2: float, P: float, V: float, V_P: float):
        # [.pyeqn] P = V_P * (H_2 - H_1) / (V - V_P)
        result = []
        H_1 = H_2 - P*V/V_P + P
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_8__A_C(H_2: float, V_P: float):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        A_C = V_P/H_2
        result.append(A_C)
        return A_C

    @staticmethod
    def eqn_3_8__H_2(A_C: float, V_P: float):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        H_2 = V_P/A_C
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_8__V_P(A_C: float, H_2: float):
        # [.pyeqn] V_P = A_C * H_2
        result = []
        V_P = A_C*H_2
        result.append(V_P)
        return V_P

    @staticmethod
    def eqn_3_9__H_2(A_C: float, H_1: float, P: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_2 = (A_C*(H_1 - P) - sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        H_2 = (A_C*(H_1 - P) + sqrt(A_C*(A_C*H_1**2 - 2*A_C*H_1*P + A_C*P**2 + 4*P*V)))/(2*A_C)
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_9__V(A_C: float, H_1: float, H_2: float, P: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        V = A_C*H_2*(-H_1 + H_2 + P)/P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_9__P(A_C: float, H_1: float, H_2: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        P = A_C*H_2*(H_1 - H_2)/(A_C*H_2 - V)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_9__H_1(A_C: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        H_1 = H_2 + P - P*V/(A_C*H_2)
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_9__A_C(H_1: float, H_2: float, P: float, V: float):
        # [.pyeqn] P = A_C * H_2 * (H_2 - H_1) / (V - A_C * H_2)
        result = []
        A_C = P*V/(H_2*(-H_1 + H_2 + P))
        result.append(A_C)
        return A_C

    @staticmethod
    def eqn_3_11__A_C(H_2: float, P: float, V: float):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        A_C = P*V/H_2**2
        result.append(A_C)
        return A_C

    @staticmethod
    def eqn_3_11__H_2(A_C: float, P: float, V: float):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        H_2 = -sqrt(P*V/A_C)
        result.append(H_2)
        H_2 = sqrt(P*V/A_C)
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_11__P(A_C: float, H_2: float, V: float):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        P = A_C*H_2**2/V
        result.append(P)
        return P

    @staticmethod
    def eqn_3_11__V(A_C: float, H_2: float, P: float):
        # [.pyeqn] P = A_C / V * (H_2) ** 2
        result = []
        V = A_C*H_2**2/P
        result.append(V)
        return V

    @staticmethod
    def eqn_3_12__H_2(KAPPA_1: float, P: float):
        # [.pyeqn] P = KAPPA_1 * H_2 ** 2
        result = []
        H_2 = -sqrt(P/KAPPA_1)
        result.append(H_2)
        H_2 = sqrt(P/KAPPA_1)
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_12__P(H_2: float, KAPPA_1: float):
        # [.pyeqn] P = KAPPA_1 * H_2 ** 2
        result = []
        P = H_2**2*KAPPA_1
        result.append(P)
        return P

    @staticmethod
    def eqn_3_12__KAPPA_1(H_2: float, P: float):
        # [.pyeqn] P = KAPPA_1 * H_2 ** 2
        result = []
        KAPPA_1 = P/H_2**2
        result.append(KAPPA_1)
        return KAPPA_1

    @staticmethod
    def eqn_3_13__KAPPA_2(H_1: float, H_2: float, P: float):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        KAPPA_2 = -P/(H_1 - H_2)
        result.append(KAPPA_2)
        return KAPPA_2

    @staticmethod
    def eqn_3_13__P(H_1: float, H_2: float, KAPPA_2: float):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        P = KAPPA_2*(-H_1 + H_2)
        result.append(P)
        return P

    @staticmethod
    def eqn_3_13__H_1(H_2: float, KAPPA_2: float, P: float):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_1 = H_2 - P/KAPPA_2
        result.append(H_1)
        return H_1

    @staticmethod
    def eqn_3_13__H_2(H_1: float, KAPPA_2: float, P: float):
        # [.pyeqn] P = KAPPA_2 * (H_2 - H_1)
        result = []
        H_2 = H_1 + P/KAPPA_2
        result.append(H_2)
        return H_2

    @staticmethod
    def eqn_3_15__V_PMIN():
        # [.pyeqn] V_PMIN = 3.141592653589793 / 4
        result = []
        V_PMIN = 0.785398163397448
        result.append(V_PMIN)
        return V_PMIN

    @staticmethod
    def eqn_3_16__V_div_V_P_MAX():
        # [.pyeqn] V_div_V_P_MAX = 200000 / (3.141592653589793 / 4)
        result = []
        V_div_V_P_MAX = 254647.908947033
        result.append(V_div_V_P_MAX)
        return V_div_V_P_MAX

    @staticmethod
    def eqn_3_17__P_MIN():
        # [.pyeqn] P_MIN = (3.141592653589793 / 4) / (200000)
        result = []
        P_MIN = 0.00000392699081698724
        result.append(P_MIN)
        return P_MIN


class AirLeak:

    @staticmethod
    def eqn_4_7__W(W_T: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W = W_T - sum_individual_leak_rates
        result.append(W)
        return W

    @staticmethod
    def eqn_4_7__W_T(W: float, sum_individual_leak_rates: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        W_T = W + sum_individual_leak_rates
        result.append(W_T)
        return W_T

    @staticmethod
    def eqn_4_7__sum_individual_leak_rates(W: float, W_T: float):
        # [.pyeqn] W_T = W + sum_individual_leak_rates
        result = []
        sum_individual_leak_rates = -W + W_T
        result.append(sum_individual_leak_rates)
        return sum_individual_leak_rates

    @staticmethod
    def eqn_4_10__del_P(T: float, V: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        del_P = 0.319795330988168*T*leakage*t/V
        result.append(del_P)
        return del_P

    @staticmethod
    def eqn_4_10__V(T: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        V = 0.319795330988168*T*leakage*t/del_P
        result.append(V)
        return V

    @staticmethod
    def eqn_4_10__leakage(T: float, V: float, del_P: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        leakage = 3.127*V*del_P/(T*t)
        result.append(leakage)
        return leakage

    @staticmethod
    def eqn_4_10__t(T: float, V: float, del_P: float, leakage: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        t = 3.127*V*del_P/(T*leakage)
        result.append(t)
        return t

    @staticmethod
    def eqn_4_10__T(V: float, del_P: float, leakage: float, t: float):
        # [.pyeqn] leakage = 0.0059 * V * del_P / t * 530 / T   lb/hr
        result = []
        T = 3.127*V*del_P/(leakage*t)
        result.append(T)
        return T


class ProcessAppI:

    @staticmethod
    def eqn_5_1__K_i(x_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        K_i = y_i/x_i
        result.append(K_i)
        return K_i

    @staticmethod
    def eqn_5_1__x_i(K_i: float, y_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        x_i = y_i/K_i
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_5_1__y_i(K_i: float, x_i: float):
        # [.pyeqn] K_i = y_i / x_i
        result = []
        y_i = K_i*x_i
        result.append(y_i)
        return y_i

    @staticmethod
    def eqn_5_2__x_1(K_1: float, K_2: float, alpha_1_2: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        pass # unable to solve

    @staticmethod
    def eqn_5_2__y_2(K_1: float, K_2: float, alpha_1_2: float, x_1: float, x_2: float, y_1: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        pass # unable to solve

    @staticmethod
    def eqn_5_2__alpha_1_2(K_1: float, K_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        alpha_1_2 = K_1/K_2
        result.append(alpha_1_2)
        return alpha_1_2

    @staticmethod
    def eqn_5_2__K_1(K_2: float, alpha_1_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_1 = K_2*alpha_1_2
        result.append(K_1)
        return K_1

    @staticmethod
    def eqn_5_2__x_2(K_1: float, K_2: float, alpha_1_2: float, x_1: float, y_1: float, y_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        pass # unable to solve

    @staticmethod
    def eqn_5_2__y_1(K_1: float, K_2: float, alpha_1_2: float, x_1: float, x_2: float, y_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        pass # unable to solve

    @staticmethod
    def eqn_5_2__K_2(K_1: float, alpha_1_2: float, x_1: float, x_2: float, y_1: float, y_2: float):
        # [.pyeqn] alpha_1_2 = K_1 / K_2 = y_1 * x_2 / (y_2 * x_1)
        result = []
        K_2 = K_1/alpha_1_2
        result.append(K_2)
        return K_2

    @staticmethod
    def eqn_5_3__x_i(P_0_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        x_i = p_i/P_0_i
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_5_3__p_i(P_0_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        p_i = P_0_i*x_i
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_5_3__P_0_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_0_i
        result = []
        P_0_i = p_i/x_i
        result.append(P_0_i)
        return P_0_i

    @staticmethod
    def eqn_5_4__P(P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P = P_0_i*x_i/y_i
        result.append(P)
        return P

    @staticmethod
    def eqn_5_4__x_i(P: float, P_0_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        x_i = P*y_i/P_0_i
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_5_4__y_i(P: float, P_0_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        y_i = P_0_i*x_i/P
        result.append(y_i)
        return y_i

    @staticmethod
    def eqn_5_4__P_0_i(P: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * P_0_i
        result = []
        P_0_i = P*y_i/x_i
        result.append(P_0_i)
        return P_0_i

    @staticmethod
    def eqn_5_5__alpha_12(P_0_1: float, P_0_2: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        alpha_12 = P_0_1/P_0_2
        result.append(alpha_12)
        return alpha_12

    @staticmethod
    def eqn_5_5__P_0_1(P_0_2: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_1 = P_0_2*alpha_12
        result.append(P_0_1)
        return P_0_1

    @staticmethod
    def eqn_5_5__P_0_2(P_0_1: float, alpha_12: float):
        # [.pyeqn] alpha_12 = P_0_1 / P_0_2
        result = []
        P_0_2 = P_0_1/alpha_12
        result.append(P_0_2)
        return P_0_2

    @staticmethod
    def eqn_5_6__gamma_i(P_0_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        gamma_i = p_i/(P_0_i*x_i)
        result.append(gamma_i)
        return gamma_i

    @staticmethod
    def eqn_5_6__x_i(P_0_i: float, gamma_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        x_i = p_i/(P_0_i*gamma_i)
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_5_6__p_i(P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        p_i = P_0_i*gamma_i*x_i
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_5_6__P_0_i(gamma_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * gamma_i * P_0_i
        result = []
        P_0_i = p_i/(gamma_i*x_i)
        result.append(P_0_i)
        return P_0_i

    @staticmethod
    def eqn_5_7__y_i(P: float, P_0_i: float, gamma_i: float, x_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        y_i = P_0_i*gamma_i*x_i/P
        result.append(y_i)
        return y_i

    @staticmethod
    def eqn_5_7__gamma_i(P: float, P_0_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        gamma_i = P*y_i/(P_0_i*x_i)
        result.append(gamma_i)
        return gamma_i

    @staticmethod
    def eqn_5_7__P(P_0_i: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P = P_0_i*gamma_i*x_i/y_i
        result.append(P)
        return P

    @staticmethod
    def eqn_5_7__x_i(P: float, P_0_i: float, gamma_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        x_i = P*y_i/(P_0_i*gamma_i)
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_5_7__P_0_i(P: float, gamma_i: float, x_i: float, y_i: float):
        # [.pyeqn] y_i * P = x_i * gamma_i * P_0_i
        result = []
        P_0_i = P*y_i/(gamma_i*x_i)
        result.append(P_0_i)
        return P_0_i

    @staticmethod
    def eqn_5_8__P_0_1(P_0_2: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_1 = P_0_2*alpha_12*gamma_2/gamma_1
        result.append(P_0_1)
        return P_0_1

    @staticmethod
    def eqn_5_8__gamma_1(P_0_1: float, P_0_2: float, alpha_12: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_1 = P_0_2*alpha_12*gamma_2/P_0_1
        result.append(gamma_1)
        return gamma_1

    @staticmethod
    def eqn_5_8__gamma_2(P_0_1: float, P_0_2: float, alpha_12: float, gamma_1: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        gamma_2 = P_0_1*gamma_1/(P_0_2*alpha_12)
        result.append(gamma_2)
        return gamma_2

    @staticmethod
    def eqn_5_8__P_0_2(P_0_1: float, alpha_12: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        P_0_2 = P_0_1*gamma_1/(alpha_12*gamma_2)
        result.append(P_0_2)
        return P_0_2

    @staticmethod
    def eqn_5_8__alpha_12(P_0_1: float, P_0_2: float, gamma_1: float, gamma_2: float):
        # [.pyeqn] alpha_12 = gamma_1 * P_0_1 / (gamma_2 * P_0_2)
        result = []
        alpha_12 = P_0_1*gamma_1/(P_0_2*gamma_2)
        result.append(alpha_12)
        return alpha_12

    @staticmethod
    def eqn_5_9__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return D

    @staticmethod
    def eqn_5_9__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return L_0

    @staticmethod
    def eqn_5_9__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = L_0 / (L_0 + D)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return V_1

    @staticmethod
    def eqn_5_10a__D(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        D = -L_0 + V_1
        result.append(D)
        return D

    @staticmethod
    def eqn_5_10a__L_0(D: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        L_0 = 0
        result.append(L_0)
        L_0 = -D + V_1
        result.append(L_0)
        return L_0

    @staticmethod
    def eqn_5_10a__V_1(D: float, L_0: float):
        # [.pyeqn] L_0 / V_1 = (L_0 / D) / (L_0 / D + 1)
        result = []
        V_1 = D + L_0
        result.append(V_1)
        return V_1

    @staticmethod
    def eqn_5_10b__L_0(R: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        L_0 = R*V_1/(R + 1)
        result.append(L_0)
        return L_0

    @staticmethod
    def eqn_5_10b__V_1(L_0: float, R: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        V_1 = L_0 + L_0/R
        result.append(V_1)
        return V_1

    @staticmethod
    def eqn_5_10b__R(L_0: float, V_1: float):
        # [.pyeqn] L_0 / V_1 = R / (R + 1)
        result = []
        R = -L_0/(L_0 - V_1)
        result.append(R)
        return R

    @staticmethod
    def eqn_5_11__L_N(B: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        L_N = B + V_0
        result.append(L_N)
        return L_N

    @staticmethod
    def eqn_5_11__B(L_N: float, V_0: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        B = L_N - V_0
        result.append(B)
        return B

    @staticmethod
    def eqn_5_11__V_0(B: float, L_N: float):
        # [.pyeqn] L_N / V_0 = (V_0 + B) / V_0
        result = []
        V_0 = -B + L_N
        result.append(V_0)
        return V_0

    @staticmethod
    def eqn_5_12__T(E: float, N_ES: float, N_t: float):
        # [.pyeqn] N_t = N_ES / E ** T
        result = []
        T = log(N_ES/N_t)
        result.append(T)
        return T

    @staticmethod
    def eqn_5_12__N_ES(E: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / E ** T
        result = []
        N_ES = N_t*exp(T)
        result.append(N_ES)
        return N_ES

    @staticmethod
    def eqn_5_12__N_t(E: float, N_ES: float, T: float):
        # [.pyeqn] N_t = N_ES / E ** T
        result = []
        N_t = N_ES*exp(-T)
        result.append(N_t)
        return N_t

    @staticmethod
    def eqn_5_12__E(N_ES: float, N_t: float, T: float):
        # [.pyeqn] N_t = N_ES / E ** T
        pass # unable to solve

    @staticmethod
    def eqn_5_13__HETP(H_p: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        HETP = H_p/N_ES
        result.append(HETP)
        return HETP

    @staticmethod
    def eqn_5_13__N_ES(HETP: float, H_p: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        N_ES = H_p/HETP
        result.append(N_ES)
        return N_ES

    @staticmethod
    def eqn_5_13__H_p(HETP: float, N_ES: float):
        # [.pyeqn] H_p = N_ES * HETP
        result = []
        H_p = HETP*N_ES
        result.append(H_p)
        return H_p

    @staticmethod
    def eqn_5_14__M(P_0: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        M = 294.213699178261*T*W_E**2/P_0**2
        result.append(M)
        return M

    @staticmethod
    def eqn_5_14__W_E(M: float, P_0: float, T: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        W_E = 0.0583*P_0*sqrt(M/T)
        result.append(W_E)
        return W_E

    @staticmethod
    def eqn_5_14__P_0(M: float, T: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        P_0 = 17.1526586620926*W_E/sqrt(M/T)
        result.append(P_0)
        return P_0

    @staticmethod
    def eqn_5_14__T(M: float, P_0: float, W_E: float):
        # [.pyeqn] W_E = 0.0583 * P_0 * (M / T) ** 0.5
        result = []
        T = 0.00339889*M*P_0**2/W_E**2
        result.append(T)
        return T

    @staticmethod
    def eqn_5_15__P_0_1(M_1: float, M_2: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_1 = P_0_2*a_M_12/(M_2/M_1)**(2/5)
        result.append(P_0_1)
        return P_0_1

    @staticmethod
    def eqn_5_15__M_1(M_2: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_1 = -M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        M_1 = M_2/(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_1)
        return M_1

    @staticmethod
    def eqn_5_15__M_2(M_1: float, P_0_1: float, P_0_2: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        M_2 = -M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        M_2 = M_1*(P_0_2*a_M_12/P_0_1)**(5/2)
        result.append(M_2)
        return M_2

    @staticmethod
    def eqn_5_15__a_M_12(M_1: float, M_2: float, P_0_1: float, P_0_2: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        a_M_12 = P_0_1*(M_2/M_1)**(2/5)/P_0_2
        result.append(a_M_12)
        return a_M_12

    @staticmethod
    def eqn_5_15__P_0_2(M_1: float, M_2: float, P_0_1: float, a_M_12: float):
        # [.pyeqn] a_M_12 = (P_0_1) / (P_0_2) * (M_2 / M_1) ** 0.4
        result = []
        P_0_2 = P_0_1*(M_2/M_1)**(2/5)/a_M_12
        result.append(P_0_2)
        return P_0_2

    @staticmethod
    def eqn_5_16__x_i(H_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        x_i = p_i/H_i
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_5_16__p_i(H_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        p_i = H_i*x_i
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_5_16__H_i(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * H_i
        result = []
        H_i = p_i/x_i
        result.append(H_i)
        return H_i

    @staticmethod
    def eqn_5_17__x_1(lnH_2_1: float, lnH_2_3: float, lnH_2_mi: float, x_3: float):
        # [.pyeqn] ln(H_2_mi) = x_1 * ln(H_2_1) + x_3 * ln(H_2_3)
        result = []
        x_1 = (-x_3*log(H_2_3) + log(H_2_mi))/log(H_2_1)
        result.append(x_1)
        return x_1

    @staticmethod
    def eqn_5_17__x_3(lnH_2_1: float, lnH_2_3: float, lnH_2_mi: float, x_1: float):
        # [.pyeqn] ln(H_2_mi) = x_1 * ln(H_2_1) + x_3 * ln(H_2_3)
        result = []
        x_3 = (-x_1*log(H_2_1) + log(H_2_mi))/log(H_2_3)
        result.append(x_3)
        return x_3

    @staticmethod
    def eqn_5_17__lnH_2_mi(lnH_2_1: float, lnH_2_3: float, x_1: float, x_3: float):
        # [.pyeqn] ln(H_2_mi) = x_1 * ln(H_2_1) + x_3 * ln(H_2_3)
        pass # unable to solve

    @staticmethod
    def eqn_5_17__lnH_2_3(lnH_2_1: float, lnH_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] ln(H_2_mi) = x_1 * ln(H_2_1) + x_3 * ln(H_2_3)
        pass # unable to solve

    @staticmethod
    def eqn_5_17__lnH_2_1(lnH_2_3: float, lnH_2_mi: float, x_1: float, x_3: float):
        # [.pyeqn] ln(H_2_mi) = x_1 * ln(H_2_1) + x_3 * ln(H_2_3)
        pass # unable to solve


class ProcessAppIi:

    @staticmethod
    def eqn_6_1__del_h_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        result = []
        del_h_v = (T_1*c_p*w_1 - T_R*c_p*w_1 + w_2*c_p(T_2 - T_R))/w_v
        result.append(del_h_v)
        return del_h_v

    @staticmethod
    def eqn_6_1__T_2(T_1: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        pass # unable to solve

    @staticmethod
    def eqn_6_1__w_1(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        result = []
        w_1 = (del_h_v*w_v - w_2*c_p(T_2 - T_R))/(c_p*(T_1 - T_R))
        result.append(w_1)
        return w_1

    @staticmethod
    def eqn_6_1__w_2(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        result = []
        w_2 = (-T_1*c_p*w_1 + T_R*c_p*w_1 + del_h_v*w_v)/c_p(T_2 - T_R)
        result.append(w_2)
        return w_2

    @staticmethod
    def eqn_6_1__c_p(T_1: float, T_2: float, T_R: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        result = []
        c_p = (del_h_v*w_v - w_2*c_p(T_2 - T_R))/(w_1*(T_1 - T_R))
        result.append(c_p)
        return c_p

    @staticmethod
    def eqn_6_1__T_1(T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        result = []
        T_1 = (T_R*c_p*w_1 + del_h_v*w_v - w_2*c_p(T_2 - T_R))/(c_p*w_1)
        result.append(T_1)
        return T_1

    @staticmethod
    def eqn_6_1__T_R(T_1: float, T_2: float, c_p: float, del_h_v: float, w_1: float, w_2: float, w_v: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        pass # unable to solve

    @staticmethod
    def eqn_6_1__w_v(T_1: float, T_2: float, T_R: float, c_p: float, del_h_v: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = w_v * del_h_v
        result = []
        w_v = (T_1*c_p*w_1 - T_R*c_p*w_1 + w_2*c_p(T_2 - T_R))/del_h_v
        result.append(w_v)
        return w_v

    @staticmethod
    def eqn_6_2__T_2(Q_v: float, T_1: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = 12000 * Q_v
        pass # unable to solve

    @staticmethod
    def eqn_6_2__w_1(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = 12000 * Q_v
        result = []
        w_1 = (12000*Q_v - w_2*c_p(T_2 - T_R))/(c_p*(T_1 - T_R))
        result.append(w_1)
        return w_1

    @staticmethod
    def eqn_6_2__w_2(Q_v: float, T_1: float, T_2: float, T_R: float, c_p: float, w_1: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = 12000 * Q_v
        result = []
        w_2 = (12000*Q_v - T_1*c_p*w_1 + T_R*c_p*w_1)/c_p(T_2 - T_R)
        result.append(w_2)
        return w_2

    @staticmethod
    def eqn_6_2__c_p(Q_v: float, T_1: float, T_2: float, T_R: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = 12000 * Q_v
        result = []
        c_p = (12000*Q_v - w_2*c_p(T_2 - T_R))/(w_1*(T_1 - T_R))
        result.append(c_p)
        return c_p

    @staticmethod
    def eqn_6_2__T_1(Q_v: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = 12000 * Q_v
        result = []
        T_1 = (12000*Q_v + T_R*c_p*w_1 - w_2*c_p(T_2 - T_R))/(c_p*w_1)
        result.append(T_1)
        return T_1

    @staticmethod
    def eqn_6_2__Q_v(T_1: float, T_2: float, T_R: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = 12000 * Q_v
        result = []
        Q_v = T_1*c_p*w_1/12000 - T_R*c_p*w_1/12000 + w_2*c_p(T_2 - T_R)/12000
        result.append(Q_v)
        return Q_v

    @staticmethod
    def eqn_6_2__T_R(Q_v: float, T_1: float, T_2: float, c_p: float, w_1: float, w_2: float):
        # [.pyeqn] w_1 * c_p * (T_1 - T_R) + w_2 * c_p (T_2 - T_R) = 12000 * Q_v
        pass # unable to solve

    @staticmethod
    def eqn_6_4__delta_h_v(Q_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        delta_h_v = 12000*Q_v/w_v
        result.append(delta_h_v)
        return delta_h_v

    @staticmethod
    def eqn_6_4__Q_v(delta_h_v: float, w_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        Q_v = delta_h_v*w_v/12000
        result.append(Q_v)
        return Q_v

    @staticmethod
    def eqn_6_4__w_v(Q_v: float, delta_h_v: float):
        # [.pyeqn] w_v = 12000 * Q_v / delta_h_v
        result = []
        w_v = 12000*Q_v/delta_h_v
        result.append(w_v)
        return w_v

    @staticmethod
    def eqn_6_6__delta_h_v(Q_r: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        delta_h_v = 24*Q_r/f_m
        result.append(delta_h_v)
        return delta_h_v

    @staticmethod
    def eqn_6_6__f_m(Q_r: float, delta_h_v: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        f_m = 24*Q_r/delta_h_v
        result.append(f_m)
        return f_m

    @staticmethod
    def eqn_6_6__Q_r(delta_h_v: float, f_m: float):
        # [.pyeqn] f_m = 24 * Q_r / delta_h_v = 0.0266 * Q_r
        result = []
        Q_r = delta_h_v*f_m/24
        result.append(Q_r)
        return Q_r

    @staticmethod
    def eqn_6_7__C_1(C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_1)
        return C_1

    @staticmethod
    def eqn_6_7__delta_h_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/m_v
        result.append(delta_h_v)
        return delta_h_v

    @staticmethod
    def eqn_6_7__T_2(C_1: float, C_2: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*m_v)/(c_p*m_b)
        result.append(T_2)
        return T_2

    @staticmethod
    def eqn_6_7__C_2(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        C_2 = (C_1*delta_h_c*m_b + c_p*m_b*(T_1 - T_2) - delta_h_v*m_v)/(delta_h_c*m_b)
        result.append(C_2)
        return C_2

    @staticmethod
    def eqn_6_7__m_b(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_b = delta_h_v*m_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return m_b

    @staticmethod
    def eqn_6_7__m_v(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        m_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/delta_h_v
        result.append(m_v)
        return m_v

    @staticmethod
    def eqn_6_7__delta_h_c(C_1: float, C_2: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*m_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return delta_h_c

    @staticmethod
    def eqn_6_7__c_p(C_1: float, C_2: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*m_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return c_p

    @staticmethod
    def eqn_6_7__T_1(C_1: float, C_2: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, m_v: float):
        # [.pyeqn] m_v * delta_h_v = m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*m_v)/(c_p*m_b)
        result.append(T_1)
        return T_1

    @staticmethod
    def eqn_6_8__C_1(T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        C_1 = (C_2*delta_h_c*m_b + c_p*m_b*(-T_1 + T_2) + delta_h_v*delta_t*w_v)/(delta_h_c*m_b)
        result.append(C_1)
        return C_1

    @staticmethod
    def eqn_6_8__delta_h_v(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_t*w_v)
        result.append(delta_h_v)
        return delta_h_v

    @staticmethod
    def eqn_6_8__T_2(C_1: float, T_1: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_2 = (T_1*c_p*m_b + delta_h_c*m_b*(C_1 - C_2) - delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_2)
        return T_2

    @staticmethod
    def eqn_6_8__m_b(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        m_b = delta_h_v*delta_t*w_v/(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)
        result.append(m_b)
        return m_b

    @staticmethod
    def eqn_6_8__delta_h_c(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        delta_h_c = (-T_1*c_p*m_b + T_2*c_p*m_b + delta_h_v*delta_t*w_v)/(m_b*(C_1 - C_2))
        result.append(delta_h_c)
        return delta_h_c

    @staticmethod
    def eqn_6_8__c_p(C_1: float, T_1: float, T_2: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        c_p = (-C_1*delta_h_c*m_b + C_2*delta_h_c*m_b + delta_h_v*delta_t*w_v)/(m_b*(T_1 - T_2))
        result.append(c_p)
        return c_p

    @staticmethod
    def eqn_6_8__T_1(C_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float, w_v: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        T_1 = (T_2*c_p*m_b + delta_h_c*m_b*(-C_1 + C_2) + delta_h_v*delta_t*w_v)/(c_p*m_b)
        result.append(T_1)
        return T_1

    @staticmethod
    def eqn_6_8__w_v(C_1: float, T_1: float, T_2: float, c_p: float, delta_h_c: float, delta_h_v: float, m_b: float):
        # [.pyeqn] w_v  = (m_b * c_p * (T_1 - T_2) + m_b * delta_h_c * (C_1 - C_2))/(delta_t * delta_h_v)
        result = []
        w_v = m_b*(C_1*delta_h_c - C_2*delta_h_c + T_1*c_p - T_2*c_p)/(delta_h_v*delta_t)
        result.append(w_v)
        return w_v

    @staticmethod
    def eqn_6_9__delta_P(A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        delta_P = A*dV_dt*r_M/(A**2 - dV_dt*m*mu*r)
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_6_9__dV_dt(A: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
        result.append(dV_dt)
        return dV_dt

    @staticmethod
    def eqn_6_9__r_M(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
        result.append(r_M)
        return r_M

    @staticmethod
    def eqn_6_9__r(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
        result.append(r)
        return r

    @staticmethod
    def eqn_6_9__m(A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
        result.append(m)
        return m

    @staticmethod
    def eqn_6_9__A(dV_dt: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        A = (dV_dt*r_M - sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        A = (dV_dt*r_M + sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        return A

    @staticmethod
    def eqn_6_9__mu(A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_6_10__s(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        s = log(A*delta_P/(dV_dt*mu*r_c*tau))/log(delta_P)
        result.append(s)
        return s

    @staticmethod
    def eqn_6_10__dV_dt(A: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        dV_dt = A*delta_P**(1 - s)/(mu*r_c*tau)
        result.append(dV_dt)
        return dV_dt

    @staticmethod
    def eqn_6_10__tau(A: float, dV_dt: float, delta_P: float, mu: float, r_c: float, s: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        tau = A*delta_P**(1 - s)/(dV_dt*mu*r_c)
        result.append(tau)
        return tau

    @staticmethod
    def eqn_6_10__delta_P(A: float, dV_dt: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        delta_P = (dV_dt*mu*r_c*tau/A)**(-1/(s - 1))
        result.append(delta_P)
        return delta_P

    @staticmethod
    def eqn_6_10__r_c(A: float, dV_dt: float, delta_P: float, mu: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        r_c = A*delta_P**(1 - s)/(dV_dt*mu*tau)
        result.append(r_c)
        return r_c

    @staticmethod
    def eqn_6_10__A(dV_dt: float, delta_P: float, mu: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        A = dV_dt*delta_P**(s - 1)*mu*r_c*tau
        result.append(A)
        return A

    @staticmethod
    def eqn_6_10__mu(A: float, dV_dt: float, delta_P: float, r_c: float, s: float, tau: float):
        # [.pyeqn] dV_dt = (A * delta_P**(1 - s) ) / ( mu * tau * r_c )
        result = []
        mu = A*delta_P**(1 - s)/(dV_dt*r_c*tau)
        result.append(mu)
        return mu

    @staticmethod
    def eqn_6_11a__h_d(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        h_d = delta_h_i*delta_m*m_b/(A_d*delta_T*t_R)
        result.append(h_d)
        return h_d

    @staticmethod
    def eqn_6_11a__delta_h_i(A_d: float, delta_T: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_h_i = A_d*delta_T*h_d*t_R/(delta_m*m_b)
        result.append(delta_h_i)
        return delta_h_i

    @staticmethod
    def eqn_6_11a__delta_T(A_d: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_T = delta_h_i*delta_m*m_b/(A_d*h_d*t_R)
        result.append(delta_T)
        return delta_T

    @staticmethod
    def eqn_6_11a__A_d(delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        A_d = delta_h_i*delta_m*m_b/(delta_T*h_d*t_R)
        result.append(A_d)
        return A_d

    @staticmethod
    def eqn_6_11a__m_b(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        m_b = A_d*delta_T*h_d*t_R/(delta_h_i*delta_m)
        result.append(m_b)
        return m_b

    @staticmethod
    def eqn_6_11a__delta_m(A_d: float, delta_T: float, delta_h_i: float, h_d: float, m_b: float, t_R: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        delta_m = A_d*delta_T*h_d*t_R/(delta_h_i*m_b)
        result.append(delta_m)
        return delta_m

    @staticmethod
    def eqn_6_11a__t_R(A_d: float, delta_T: float, delta_h_i: float, delta_m: float, h_d: float, m_b: float):
        # [.pyeqn] t_R = delta_h_i * m_b * delta_m / (A_d * h_d * delta_T)
        result = []
        t_R = delta_h_i*delta_m*m_b/(A_d*delta_T*h_d)
        result.append(t_R)
        return t_R


class Precondensors:

    @staticmethod
    def eqn_7_1__P(p_i: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        P = p_i/y_i
        result.append(P)
        return P

    @staticmethod
    def eqn_7_1__p_i(P: float, y_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        p_i = P*y_i
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_7_1__y_i(P: float, p_i: float):
        # [.pyeqn] y_i = p_i / P
        result = []
        y_i = p_i/P
        result.append(y_i)
        return y_i

    @staticmethod
    def eqn_7_2__P_i_0(p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        P_i_0 = p_i/x_i
        result.append(P_i_0)
        return P_i_0

    @staticmethod
    def eqn_7_2__x_i(P_i_0: float, p_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        x_i = p_i/P_i_0
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_7_2__p_i(P_i_0: float, x_i: float):
        # [.pyeqn] p_i = x_i * P_i_0
        result = []
        p_i = P_i_0*x_i
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_7_3__P_i_0(epsilon_i: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        P_i_0 = p_i/(epsilon_i*x_i)
        result.append(P_i_0)
        return P_i_0

    @staticmethod
    def eqn_7_3__x_i(P_i_0: float, epsilon_i: float, p_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        x_i = p_i/(P_i_0*epsilon_i)
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_7_3__p_i(P_i_0: float, epsilon_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        p_i = P_i_0*epsilon_i*x_i
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_7_3__epsilon_i(P_i_0: float, p_i: float, x_i: float):
        # [.pyeqn] p_i = x_i * epsilon_i * P_i_0
        result = []
        epsilon_i = p_i/(P_i_0*x_i)
        result.append(epsilon_i)
        return epsilon_i

    @staticmethod
    def eqn_7_4__p_c(P: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_c = P - p_nc
        result.append(p_c)
        return p_c

    @staticmethod
    def eqn_7_4__P(p_c: float, p_nc: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        P = p_c + p_nc
        result.append(P)
        return P

    @staticmethod
    def eqn_7_4__p_nc(P: float, p_c: float):
        # [.pyeqn] p_nc = P - p_c
        result = []
        p_nc = P - p_c
        result.append(p_nc)
        return p_nc

    @staticmethod
    def eqn_7_4__p(P_c: float, n_i: float, n_nc: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc = p_i / (p - P_c)
        pass # unable to solve

    @staticmethod
    def eqn_7_4__p_i(P_c: float, n_i: float, n_nc: float, p: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc = p_i / (p - P_c)
        result = []
        p_i = n_i*p_nc/n_nc
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_7_4__p_nc(P_c: float, n_i: float, n_nc: float, p: float, p_i: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc = p_i / (p - P_c)
        result = []
        p_nc = n_nc*p_i/n_i
        result.append(p_nc)
        return p_nc

    @staticmethod
    def eqn_7_4__P_c(n_i: float, n_nc: float, p: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc = p_i / (p - P_c)
        pass # unable to solve

    @staticmethod
    def eqn_7_4__n_i(P_c: float, n_nc: float, p: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc = p_i / (p - P_c)
        result = []
        n_i = n_nc*p_i/p_nc
        result.append(n_i)
        return n_i

    @staticmethod
    def eqn_7_4__n_nc(P_c: float, n_i: float, p: float, p_i: float, p_nc: float):
        # [.pyeqn] n_i / n_nc = p_i / p_nc = p_i / (p - P_c)
        result = []
        n_nc = n_i*p_nc/p_i
        result.append(n_nc)
        return n_nc

    @staticmethod
    def eqn_7_5__P(N_i: float, N_nc: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P = P_c + N_nc*p_i/N_i
        result.append(P)
        return P

    @staticmethod
    def eqn_7_5__p_i(N_i: float, N_nc: float, P: float, P_c: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        p_i = N_i*(P - P_c)/N_nc
        result.append(p_i)
        return p_i

    @staticmethod
    def eqn_7_5__N_i(N_nc: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_i = N_nc*p_i/(P - P_c)
        result.append(N_i)
        return N_i

    @staticmethod
    def eqn_7_5__P_c(N_i: float, N_nc: float, P: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        P_c = P - N_nc*p_i/N_i
        result.append(P_c)
        return P_c

    @staticmethod
    def eqn_7_5__N_nc(N_i: float, P: float, P_c: float, p_i: float):
        # [.pyeqn] N_i = N_nc * (p_i) / (P - P_c)
        result = []
        N_nc = N_i*(P - P_c)/p_i
        result.append(N_nc)
        return N_nc

    @staticmethod
    def eqn_7_6__P(M: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*x_i/(29*W_i) + p_c
        result.append(P)
        return P

    @staticmethod
    def eqn_7_6__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air)
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_7_6__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*x_i/(29*W_i) + P
        result.append(p_c)
        return p_c

    @staticmethod
    def eqn_7_6__P_i_0(M: float, P: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*x_i)
        result.append(P_i_0)
        return P_i_0

    @staticmethod
    def eqn_7_6__W_air(M: float, P: float, P_i_0: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*x_i)
        result.append(W_air)
        return W_air

    @staticmethod
    def eqn_7_6__W_i(M: float, P: float, P_i_0: float, W_air: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*x_i/(29*(P - p_c))
        result.append(W_i)
        return W_i

    @staticmethod
    def eqn_7_6__M(P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*x_i)
        result.append(M)
        return M

    @staticmethod
    def eqn_7_7__P(M: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P = M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + p_c
        result.append(P)
        return P

    @staticmethod
    def eqn_7_7__x_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        x_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*epsilon_i)
        result.append(x_i)
        return x_i

    @staticmethod
    def eqn_7_7__epsilon_i(M: float, P: float, P_i_0: float, W_air: float, W_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        epsilon_i = 29*W_i*(P - p_c)/(M*P_i_0*W_air*x_i)
        result.append(epsilon_i)
        return epsilon_i

    @staticmethod
    def eqn_7_7__p_c(M: float, P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        p_c = -M*P_i_0*W_air*epsilon_i*x_i/(29*W_i) + P
        result.append(p_c)
        return p_c

    @staticmethod
    def eqn_7_7__P_i_0(M: float, P: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        P_i_0 = 29*W_i*(P - p_c)/(M*W_air*epsilon_i*x_i)
        result.append(P_i_0)
        return P_i_0

    @staticmethod
    def eqn_7_7__W_air(M: float, P: float, P_i_0: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_air = 29*W_i*(P - p_c)/(M*P_i_0*epsilon_i*x_i)
        result.append(W_air)
        return W_air

    @staticmethod
    def eqn_7_7__W_i(M: float, P: float, P_i_0: float, W_air: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        W_i = M*P_i_0*W_air*epsilon_i*x_i/(29*(P - p_c))
        result.append(W_i)
        return W_i

    @staticmethod
    def eqn_7_7__M(P: float, P_i_0: float, W_air: float, W_i: float, epsilon_i: float, p_c: float, x_i: float):
        # [.pyeqn] W_i = W_air * (M * x_i * epsilon_i * P_i_0) / (29 * (P - p_c))
        result = []
        M = 29*W_i*(P - p_c)/(P_i_0*W_air*epsilon_i*x_i)
        result.append(M)
        return M

    @staticmethod
    def eqn_7_8__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        c_p = Q_condensor_heat_duty/(L_c*del_T)
        result.append(c_p)
        return c_p

    @staticmethod
    def eqn_7_8__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        L_c = Q_condensor_heat_duty/(c_p*del_T)
        result.append(L_c)
        return L_c

    @staticmethod
    def eqn_7_8__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(L_c*c_p)
        result.append(del_T)
        return del_T

    @staticmethod
    def eqn_7_8__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T)
        result = []
        Q_condensor_heat_duty = L_c*c_p*del_T
        result.append(Q_condensor_heat_duty)
        return Q_condensor_heat_duty

    @staticmethod
    def eqn_7_9__L_c(Q_condensor_heat_duty: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        L_c = 0.124688279301746*Q_condensor_heat_duty/(c_p*del_T*rho)
        result.append(L_c)
        return L_c

    @staticmethod
    def eqn_7_9__rho(L_c: float, Q_condensor_heat_duty: float, c_p: float, del_T: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        rho = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*del_T)
        result.append(rho)
        return rho

    @staticmethod
    def eqn_7_9__c_p(L_c: float, Q_condensor_heat_duty: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        c_p = 0.124688279301746*Q_condensor_heat_duty/(L_c*del_T*rho)
        result.append(c_p)
        return c_p

    @staticmethod
    def eqn_7_9__del_T(L_c: float, Q_condensor_heat_duty: float, c_p: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        del_T = 0.124688279301746*Q_condensor_heat_duty/(L_c*c_p*rho)
        result.append(del_T)
        return del_T

    @staticmethod
    def eqn_7_9__Q_condensor_heat_duty(L_c: float, c_p: float, del_T: float, rho: float):
        # [.pyeqn] L_c = Q_condensor_heat_duty / (c_p * del_T * rho * 8.02)
        result = []
        Q_condensor_heat_duty = 8.02*L_c*c_p*del_T*rho
        result.append(Q_condensor_heat_duty)
        return Q_condensor_heat_duty

    @staticmethod
    def eqn_7_10__del_T(L_c_P: float, Q_condensor_heat_duty: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        del_T = Q_condensor_heat_duty/(500*L_c_P)
        result.append(del_T)
        return del_T

    @staticmethod
    def eqn_7_10__Q_condensor_heat_duty(L_c_P: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        Q_condensor_heat_duty = 500*L_c_P*del_T
        result.append(Q_condensor_heat_duty)
        return Q_condensor_heat_duty

    @staticmethod
    def eqn_7_10__L_c_P(Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] L_c_P = Q_condensor_heat_duty / (500 * del_T)
        result = []
        L_c_P = Q_condensor_heat_duty/(500*del_T)
        result.append(L_c_P)
        return L_c_P

    @staticmethod
    def eqn_7_11__V_c(Q_condensor_heat_duty: float, U_v: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        V_c = Q_condensor_heat_duty/(U_v*del_T_LM)
        result.append(V_c)
        return V_c

    @staticmethod
    def eqn_7_11__U_v(Q_condensor_heat_duty: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        U_v = Q_condensor_heat_duty/(V_c*del_T_LM)
        result.append(U_v)
        return U_v

    @staticmethod
    def eqn_7_11__del_T_LM(Q_condensor_heat_duty: float, U_v: float, V_c: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(U_v*V_c)
        result.append(del_T_LM)
        return del_T_LM

    @staticmethod
    def eqn_7_11__Q_condensor_heat_duty(U_v: float, V_c: float, del_T_LM: float):
        # [.pyeqn] V_c = Q_condensor_heat_duty / (U_v * del_T_LM)
        result = []
        Q_condensor_heat_duty = U_v*V_c*del_T_LM
        result.append(Q_condensor_heat_duty)
        return Q_condensor_heat_duty

    @staticmethod
    def eqn_7_12__U(A: float, Q_condensor_heat_duty: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        U = Q_condensor_heat_duty/(A*del_T)
        result.append(U)
        return U

    @staticmethod
    def eqn_7_12__del_T(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        del_T = Q_condensor_heat_duty/(A*U)
        result.append(del_T)
        return del_T

    @staticmethod
    def eqn_7_12__A(Q_condensor_heat_duty: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        A = Q_condensor_heat_duty/(U*del_T)
        result.append(A)
        return A

    @staticmethod
    def eqn_7_12__Q_condensor_heat_duty(A: float, U: float, del_T: float):
        # [.pyeqn] Q_condensor_heat_duty = U * A * del_T
        result = []
        Q_condensor_heat_duty = A*U*del_T
        result.append(Q_condensor_heat_duty)
        return Q_condensor_heat_duty

    @staticmethod
    def eqn_7_14a__U(A: float, Q_condensor_heat_duty: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        U = Q_condensor_heat_duty/(A*del_T_LM)
        result.append(U)
        return U

    @staticmethod
    def eqn_7_14a__del_T_LM(A: float, Q_condensor_heat_duty: float, U: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        del_T_LM = Q_condensor_heat_duty/(A*U)
        result.append(del_T_LM)
        return del_T_LM

    @staticmethod
    def eqn_7_14a__A(Q_condensor_heat_duty: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        A = Q_condensor_heat_duty/(U*del_T_LM)
        result.append(A)
        return A

    @staticmethod
    def eqn_7_14a__Q_condensor_heat_duty(A: float, U: float, del_T_LM: float):
        # [.pyeqn] A = Q_condensor_heat_duty / (U * del_T_LM)
        result = []
        Q_condensor_heat_duty = A*U*del_T_LM
        result.append(Q_condensor_heat_duty)
        return Q_condensor_heat_duty

    @staticmethod
    def eqn_7_14b__del_T_1(A: float, Q_condensor_heat_duty: float, U: float, del_T_2: float, lndel_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_1 = del_T_2 + exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_1)
        return del_T_1

    @staticmethod
    def eqn_7_14b__U(A: float, Q_condensor_heat_duty: float, del_T_1: float, del_T_2: float, lndel_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        U = Q_condensor_heat_duty/(A*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(U)
        return U

    @staticmethod
    def eqn_7_14b__lndel_T_1(A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        pass # unable to solve

    @staticmethod
    def eqn_7_14b__del_T_2(A: float, Q_condensor_heat_duty: float, U: float, del_T_1: float, lndel_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        del_T_2 = del_T_1 - exp(LambertW(Q_condensor_heat_duty/(A*U)))
        result.append(del_T_2)
        return del_T_2

    @staticmethod
    def eqn_7_14b__A(Q_condensor_heat_duty: float, U: float, del_T_1: float, del_T_2: float, lndel_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        A = Q_condensor_heat_duty/(U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2))
        result.append(A)
        return A

    @staticmethod
    def eqn_7_14b__Q_condensor_heat_duty(A: float, U: float, del_T_1: float, del_T_2: float, lndel_T_1: float):
        # [.pyeqn] A = (Q_condensor_heat_duty / (U * (del_T_1 - del_T_2))) / ln(del_T_1 - del_T_2)
        result = []
        Q_condensor_heat_duty = A*U*(del_T_1 - del_T_2)*log(del_T_1 - del_T_2)
        result.append(Q_condensor_heat_duty)
        return Q_condensor_heat_duty

    @staticmethod
    def eqn_7_15__U(sum_R: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        U = 1/sum_R
        result.append(U)
        return U

    @staticmethod
    def eqn_7_15__sum_R(U: float):
        # [.pyeqn] 1 / U = sum_R
        result = []
        sum_R = 1/U
        result.append(sum_R)
        return sum_R

    @staticmethod
    def eqn_7_16__h_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_0 = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_f_0*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_0)
        return h_0

    @staticmethod
    def eqn_7_16__k_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_0*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(k_w)
        return k_w

    @staticmethod
    def eqn_7_16__D_0(D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_f_0*U_0*h_0 - U_0 + h_0)/(U_0*h_0*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return D_0

    @staticmethod
    def eqn_7_16__U_0(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_0*h_i*k_w/(D_0*D_LM*R_fi*h_0*h_i*k_w + D_0*D_LM*h_0*k_w + D_0*D_i*h_0*h_i*x_w + D_LM*D_i*R_f_0*h_0*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return U_0

    @staticmethod
    def eqn_7_16__x_w(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_f_0*k_w/D_0 - D_LM*k_w/(D_0*h_0) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return x_w

    @staticmethod
    def eqn_7_16__D_LM(D_0: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_0*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_0*h_i + D_0*U_0*h_0 + D_i*R_f_0*U_0*h_0*h_i + D_i*U_0*h_i - D_i*h_0*h_i))
        result.append(D_LM)
        return D_LM

    @staticmethod
    def eqn_7_16__R_fi(D_0: float, D_LM: float, D_i: float, R_f_0: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_f_0/D_0 - D_i/(D_0*h_0) + D_i/(D_0*U_0)
        result.append(R_fi)
        return R_fi

    @staticmethod
    def eqn_7_16__D_i(D_0: float, D_LM: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_0*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_0*x_w + D_LM*R_f_0*U_0*h_0*k_w + D_LM*U_0*k_w - D_LM*h_0*k_w))
        result.append(D_i)
        return D_i

    @staticmethod
    def eqn_7_16__h_i(D_0: float, D_LM: float, D_i: float, R_f_0: float, R_fi: float, U_0: float, h_0: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_0*k_w/(D_0*D_LM*R_fi*U_0*h_0*k_w + D_0*D_i*U_0*h_0*x_w + D_LM*D_i*R_f_0*U_0*h_0*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_0*k_w)
        result.append(h_i)
        return h_i

    @staticmethod
    def eqn_7_16__R_f_0(D_0: float, D_LM: float, D_i: float, R_fi: float, U_0: float, h_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 =  1 / h_0 + R_f_0 + x_w * D_0 / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_f_0 = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - 1/h_0 + 1/U_0
        result.append(R_f_0)
        return R_f_0

    @staticmethod
    def eqn_7_17__R_0(R_nc: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_0 = R_nc + 1/h_c
        result.append(R_0)
        return R_0

    @staticmethod
    def eqn_7_17__R_nc(R_0: float, h_c: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        R_nc = R_0 - 1/h_c
        result.append(R_nc)
        return R_nc

    @staticmethod
    def eqn_7_17__h_c(R_0: float, R_nc: float):
        # [.pyeqn] R_0 = R_nc + 1 / h_c
        result = []
        h_c = 1/(R_0 - R_nc)
        result.append(h_c)
        return h_c

    @staticmethod
    def eqn_7_18__k_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        k_w = -D_0*D_i*U_0*h_c*h_i*x_w/(D_LM*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(k_w)
        return k_w

    @staticmethod
    def eqn_7_18__D_0(D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_0 = D_LM*D_i*h_i*k_w*(-R_fo*U_0*h_c - R_nc*U_0*h_c - U_0 + h_c)/(U_0*h_c*(D_LM*R_fi*h_i*k_w + D_LM*k_w + D_i*h_i*x_w))
        result.append(D_0)
        return D_0

    @staticmethod
    def eqn_7_18__U_0(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        U_0 = D_LM*D_i*h_c*h_i*k_w/(D_0*D_LM*R_fi*h_c*h_i*k_w + D_0*D_LM*h_c*k_w + D_0*D_i*h_c*h_i*x_w + D_LM*D_i*R_fo*h_c*h_i*k_w + D_LM*D_i*R_nc*h_c*h_i*k_w + D_LM*D_i*h_i*k_w)
        result.append(U_0)
        return U_0

    @staticmethod
    def eqn_7_18__R_nc(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_nc = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_fo - 1/h_c + 1/U_0
        result.append(R_nc)
        return R_nc

    @staticmethod
    def eqn_7_18__D_LM(D_0: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_LM = -D_0*D_i*U_0*h_c*h_i*x_w/(k_w*(D_0*R_fi*U_0*h_c*h_i + D_0*U_0*h_c + D_i*R_fo*U_0*h_c*h_i + D_i*R_nc*U_0*h_c*h_i + D_i*U_0*h_i - D_i*h_c*h_i))
        result.append(D_LM)
        return D_LM

    @staticmethod
    def eqn_7_18__R_fi(D_0: float, D_LM: float, D_i: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fi = -1/h_i - D_i*x_w/(D_LM*k_w) - D_i*R_fo/D_0 - D_i*R_nc/D_0 - D_i/(D_0*h_c) + D_i/(D_0*U_0)
        result.append(R_fi)
        return R_fi

    @staticmethod
    def eqn_7_18__D_i(D_0: float, D_LM: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        D_i = -D_0*D_LM*U_0*h_c*k_w*(R_fi*h_i + 1)/(h_i*(D_0*U_0*h_c*x_w + D_LM*R_fo*U_0*h_c*k_w + D_LM*R_nc*U_0*h_c*k_w + D_LM*U_0*k_w - D_LM*h_c*k_w))
        result.append(D_i)
        return D_i

    @staticmethod
    def eqn_7_18__h_i(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_i = -D_0*D_LM*U_0*h_c*k_w/(D_0*D_LM*R_fi*U_0*h_c*k_w + D_0*D_i*U_0*h_c*x_w + D_LM*D_i*R_fo*U_0*h_c*k_w + D_LM*D_i*R_nc*U_0*h_c*k_w + D_LM*D_i*U_0*k_w - D_LM*D_i*h_c*k_w)
        result.append(h_i)
        return h_i

    @staticmethod
    def eqn_7_18__x_w(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        x_w = -D_LM*R_fi*k_w/D_i - D_LM*k_w/(D_i*h_i) - D_LM*R_fo*k_w/D_0 - D_LM*R_nc*k_w/D_0 - D_LM*k_w/(D_0*h_c) + D_LM*k_w/(D_0*U_0)
        result.append(x_w)
        return x_w

    @staticmethod
    def eqn_7_18__R_fo(D_0: float, D_LM: float, D_i: float, R_fi: float, R_nc: float, U_0: float, h_c: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        R_fo = -D_0*R_fi/D_i - D_0/(D_i*h_i) - D_0*x_w/(D_LM*k_w) - R_nc - 1/h_c + 1/U_0
        result.append(R_fo)
        return R_fo

    @staticmethod
    def eqn_7_18__h_c(D_0: float, D_LM: float, D_i: float, R_fi: float, R_fo: float, R_nc: float, U_0: float, h_i: float, k_w: float, x_w: float):
        # [.pyeqn] 1 / U_0 = R_nc + 1 / h_c + R_fo + (x_w * D_0) / (k_w * D_LM) + R_fi * D_0 / D_i + D_0 / (h_i * D_i)
        result = []
        h_c = -D_LM*D_i*U_0*h_i*k_w/(D_0*D_LM*R_fi*U_0*h_i*k_w + D_0*D_LM*U_0*k_w + D_0*D_i*U_0*h_i*x_w + D_LM*D_i*R_fo*U_0*h_i*k_w + D_LM*D_i*R_nc*U_0*h_i*k_w - D_LM*D_i*h_i*k_w)
        result.append(h_c)
        return h_c


class SelectingPump:

    @staticmethod
    def eqn_8_1__NC(NS: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return NC

    @staticmethod
    def eqn_8_1__NS(NC: float, SCON: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return NS

    @staticmethod
    def eqn_8_1__SCON(NC: float, NS: float, installation_cost: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
