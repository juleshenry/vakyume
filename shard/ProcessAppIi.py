from kwasak import kwasak_static

class ProcessAppIi:

    @kwasak_static
    def eqn_6_09(A: float = None, dV_dt: float = None, delta_P: float = None, m: float = None, mu: float = None, r: float = None, r_M: float = None,**kwargs):
        return


    @staticmethod
    def eqn_6_09__A(dV_dt: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        A = (dV_dt*r_M - sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        A = (dV_dt*r_M + sqrt(dV_dt*(dV_dt*r_M**2 + 4*delta_P**2*m*mu*r)))/(2*delta_P)
        result.append(A)
        return result

    @staticmethod
    def eqn_6_09__dV_dt(A: float, delta_P: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        dV_dt = A**2*delta_P/(A*r_M + delta_P*m*mu*r)
        result.append(dV_dt)
        return result

    @staticmethod
    def eqn_6_09__delta_P(A: float, dV_dt: float, m: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        delta_P = A*dV_dt*r_M/(A**2 - dV_dt*m*mu*r)
        result.append(delta_P)
        return result

    @staticmethod
    def eqn_6_09__m(A: float, dV_dt: float, delta_P: float, mu: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        m = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*mu*r)
        result.append(m)
        return result

    @staticmethod
    def eqn_6_09__mu(A: float, dV_dt: float, delta_P: float, m: float, r: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        mu = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*r)
        result.append(mu)
        return result

    @staticmethod
    def eqn_6_09__r(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r_M: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r = A*(A*delta_P - dV_dt*r_M)/(dV_dt*delta_P*m*mu)
        result.append(r)
        return result

    @staticmethod
    def eqn_6_09__r_M(A: float, dV_dt: float, delta_P: float, m: float, mu: float, r: float):
        # [.pyeqn] dV_dt = (A * delta_P) / (mu * (m / A) * r * delta_P + r_M)
        result = []
        r_M = A*delta_P/dV_dt - delta_P*m*mu*r/A
        result.append(r_M)
        return result
