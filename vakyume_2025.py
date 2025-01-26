from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static

from suck_consts import * 

class VacuumTheory:

    @kwasak_static
    def eqn_1_3(m: float = None,T: float = None,k: float = None,v: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_3__m(T: float, k: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        m = 3.0*T*k/v**2
        result.append(m)
        return result

    @staticmethod
    def eqn_1_3__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_1_3__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        k = 0.333333333333333*m*v**2/T
        result.append(k)
        return result

    @staticmethod
    def eqn_1_3__v(T: float, k: float, m: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        v = -1.73205080756888*sqrt(T*k/m)
        result.append(v)
        v = 1.73205080756888*sqrt(T*k/m)
        result.append(v)
        return result

    @kwasak_static
    def eqn_1_7(T: float = None,n: float = None,V: float = None,p: float = None,R: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_7__T(R: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        T = V*p/(R*n)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_7__n(R: float, T: float, V: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        n = V*p/(R*T)
        result.append(n)
        return result

    @staticmethod
    def eqn_1_7__V(R: float, T: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        V = R*T*n/p
        result.append(V)
        return result

    @staticmethod
    def eqn_1_7__p(R: float, T: float, V: float, n: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        p = R*T*n/V
        result.append(p)
        return result

    @staticmethod
    def eqn_1_7__R(T: float, V: float, n: float, p: float):
        # [.pyeqn] p * V = n * R * T
        result = []
        R = V*p/(T*n)
        result.append(R)
        return result

    @kwasak_static
    def eqn_1_8(M: float = None,T: float = None,P: float = None,V: float = None,m: float = None,R: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_8__M(P: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        M = R*T*m/(P*V)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_8__T(M: float, P: float, R: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        T = M*P*V/(R*m)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_8__P(M: float, R: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        P = R*T*m/(M*V)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_8__V(M: float, P: float, R: float, T: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        V = R*T*m/(M*P)
        result.append(V)
        return result

    @staticmethod
    def eqn_1_8__m(M: float, P: float, R: float, T: float, V: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        m = M*P*V/(R*T)
        result.append(m)
        return result

    @staticmethod
    def eqn_1_8__R(M: float, P: float, T: float, V: float, m: float):
        # [.pyeqn] P * V = m / M * R * T
        result = []
        R = M*P*V/(T*m)
        result.append(R)
        return result

    @kwasak_static
    def eqn_1_9(T: float = None,M: float = None,rho: float = None,P: float = None,R: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_9__T(M: float, P: float, R: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        T = M*P/(R*rho)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_9__M(P: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        M = R*T*rho/P
        result.append(M)
        return result

    @staticmethod
    def eqn_1_9__rho(M: float, P: float, R: float, T: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        rho = M*P/(R*T)
        result.append(rho)
        return result

    @staticmethod
    def eqn_1_9__P(M: float, R: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        P = R*T*rho/M
        result.append(P)
        return result

    @staticmethod
    def eqn_1_9__R(M: float, P: float, T: float, rho: float):
        # [.pyeqn] rho = P * M / (R * T)
        result = []
        R = M*P/(T*rho)
        result.append(R)
        return result

    @kwasak_static
    def eqn_1_10(T_1: float = None,P_2: float = None,V_1: float = None,T_2: float = None,P_1: float = None,V_2: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_1 = P_1*T_2*V_1/(P_2*V_2)
        result.append(T_1)
        return result

    @staticmethod
    def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        P_2 = P_1*T_2*V_1/(T_1*V_2)
        result.append(P_2)
        return result

    @staticmethod
    def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_1 = P_2*T_1*V_2/(P_1*T_2)
        result.append(V_1)
        return result

    @staticmethod
    def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        T_2 = P_2*T_1*V_2/(P_1*V_1)
        result.append(T_2)
        return result

    @staticmethod
    def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        co = 'P_2, T_1, T_2, V_1, V_2'.split(',')
        print(*[(a,b,) for a,b in zip(co,[P_2, T_1, T_2, V_1, V_2])],sep='\n')
        P_1 = P_2*T_1*V_2/(T_2*V_1)
        result.append(P_1)
        return result

    @staticmethod
    def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
        # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
        result = []
        V_2 = P_1*T_2*V_1/(P_2*T_1)
        result.append(V_2)
        return result

    @kwasak_static
    def eqn_1_11(T: float = None,M: float = None,W: float = None,P: float = None,q: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_11__T(M: float, P: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        T = 738*M*P*q/(6821*W)
        result.append(T)
        return result

    @staticmethod
    def eqn_1_11__M(P: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        M = 6821*T*W/(738*P*q)
        result.append(M)
        return result

    @staticmethod
    def eqn_1_11__W(M: float, P: float, T: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        W = 738*M*P*q/(6821*T)
        result.append(W)
        return result

    @staticmethod
    def eqn_1_11__P(M: float, T: float, W: float, q: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        P = 6821*T*W/(738*M*q)
        result.append(P)
        return result

    @staticmethod
    def eqn_1_11__q(M: float, P: float, T: float, W: float):
        # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
        result = []
        q = 6821*T*W/(738*M*P)
        result.append(q)
        return result

    @kwasak_static
    def eqn_1_12(sum_partial_pressures: float = None,Total_P: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_12__sum_partial_pressures(Total_P: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        sum_partial_pressures = Total_P
        result.append(sum_partial_pressures)
        return result

    @staticmethod
    def eqn_1_12__Total_P(sum_partial_pressures: float):
        # [.pyeqn] Total_P = sum_partial_pressures
        result = []
        Total_P = sum_partial_pressures
        result.append(Total_P)
        return result

    @kwasak_static
    def eqn_1_13a(y_a: float = None,n: float = None,n_a: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_13a__y_a(n: float, n_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        y_a = n_a/n
        result.append(y_a)
        return result

    @staticmethod
    def eqn_1_13a__n(n_a: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n = n_a/y_a
        result.append(n)
        return result

    @staticmethod
    def eqn_1_13a__n_a(n: float, y_a: float):
        # [.pyeqn] y_a = n_a / n
        result = []
        n_a = n*y_a
        result.append(n_a)
        return result

    @kwasak_static
    def eqn_1_13b(p_a: float = None,y_a: float = None,P: float = None, **kwargs):
        return


    @staticmethod
    def eqn_1_13b__p_a(P: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        p_a = P*y_a
        result.append(p_a)
        return result

    @staticmethod
    def eqn_1_13b__y_a(P: float, p_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        y_a = p_a/P
        result.append(y_a)
        return result

    @staticmethod
    def eqn_1_13b__P(p_a: float, y_a: float):
        # [.pyeqn] y_a = p_a / P
        result = []
        P = p_a/y_a
        result.append(P)
        return result

print(sak_funx:=list(filter(lambda a:a.startswith('eqn') and '__'not in a,dir(VacuumTheory))))
print(funx:=list(filter(lambda a:a.startswith('eqn') and '__' in a,dir(VacuumTheory))))

# iterate all methods and fill with dummy values.
import inspect


def construir(s):
    candz= list(filter(lambda a:a.startswith(s),funx))[:2] # first two'll do
    a,b = map(lambda a:inspect.signature(getattr(VacuumTheory,a)),candz)
    woa = lambda d:str(d).replace(')','').replace('(','').replace(', ','').split(TYPE)
    tokes = set([u for u in woa(a)+woa(b)if u])
    return sorted(list(tokes))
        
sak_funx_di = {s: construir(s) for s in sak_funx}
import random 
#KWARGS {'P_1': 1.66261, 'P_2': 0.57988, 'T_1': 0.0068, 'T_2': 9.97603, 'V_1': 1.15482, 'V_2': 2.25261}
fr = { 'P_2': 0.57988, 'T_1': 0.0068, 'T_2': 9.97603, 'V_1': 1.15482, 'V_2': 2.25261}
s=VacuumTheory.eqn_1_10__P_1(**{'P_2': 0.57988, 'T_1': 0.0068, 'T_2': 9.97603, 'V_1': 1.15482, 'V_2': 2.25261}) #[0.0007710117693077729]
print(s)
1/0
for o,d in sak_funx_di.items():
    print('EQN',o)
    preset_dummy_args = {d:round(random.random() * 10,5) for d in d}
    print('KWARGS',preset_dummy_args)
    for dd in d:
        print(dd)
        nao = list(filter(lambda a:a!=dd, list(d)))
        kew = {n:preset_dummy_args[n] for n in nao}
        print('calling .. ',o+'_'+dd,)
        print('w/',kew)

        print(getattr(VacuumTheory(),o)(**kew))
    1/0
        