from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class RotaryPistonVane:
    @kwasak_static
    def eqn_11_06(
        P_0_V: float = None,
        P_D: float = None,
        P_v_0: float = None,
        S_B: float = None,
        S_D: float = None,
        p_b: float = None,
        p_g: float = None,
        p_v_max: float = None,
        **kwargs
    ):
        return

    @staticmethod
    def eqn_11_06__P_0_V(
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_0_V = (
            P_D * S_B * p_b
            + P_D * S_D * p_v_max
            - P_v_0 * S_D * p_g
            - P_v_0 * S_D * p_v_max
        ) / (P_D * S_B)
        result.append(P_0_V)
        return result

    @staticmethod
    def eqn_11_06__P_D(
        P_0_V: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_D = P_v_0 * S_D * (p_g + p_v_max) / (-P_0_V * S_B + S_B * p_b + S_D * p_v_max)
        result.append(P_D)
        return result

    @staticmethod
    def eqn_11_06__P_v_0(
        P_0_V: float,
        P_D: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        P_v_0 = (
            P_D * (-P_0_V * S_B + S_B * p_b + S_D * p_v_max) / (S_D * (p_g + p_v_max))
        )
        result.append(P_v_0)
        return result

    @staticmethod
    def eqn_11_06__S_B(
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_D: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_B = (
            S_D
            * (P_D * p_v_max - P_v_0 * p_g - P_v_0 * p_v_max)
            / (P_D * (P_0_V - p_b))
        )
        result.append(S_B)
        return result

    @staticmethod
    def eqn_11_06__S_D(
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        p_b: float,
        p_g: float,
        p_v_max: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        S_D = (
            P_D * S_B * (P_0_V - p_b) / (P_D * p_v_max - P_v_0 * p_g - P_v_0 * p_v_max)
        )
        result.append(S_D)
        return result

    @staticmethod
    def eqn_11_06__p_b(
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_g: float,
        p_v_max: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_b = (
            P_0_V * P_D * S_B
            - P_D * S_D * p_v_max
            + P_v_0 * S_D * p_g
            + P_v_0 * S_D * p_v_max
        ) / (P_D * S_B)
        result.append(p_b)
        return result

    @staticmethod
    def eqn_11_06__p_g(
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_v_max: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_g = (
            -P_0_V * P_D * S_B
            + P_D * S_B * p_b
            + P_D * S_D * p_v_max
            - P_v_0 * S_D * p_v_max
        ) / (P_v_0 * S_D)
        result.append(p_g)
        return result

    @staticmethod
    def eqn_11_06__p_v_max(
        P_0_V: float,
        P_D: float,
        P_v_0: float,
        S_B: float,
        S_D: float,
        p_b: float,
        p_g: float,
    ):
        # [.pyeqn] p_v_max = S_B / S_D * P_D * (P_0_V - p_b) / (P_D - P_v_0) + P_v_0 / (P_D - P_v_0) * p_g
        result = []
        p_v_max = (P_0_V * P_D * S_B - P_D * S_B * p_b + P_v_0 * S_D * p_g) / (
            S_D * (P_D - P_v_0)
        )
        result.append(p_v_max)
        return result

if __name__ == "__main__":
    import tru
    y = {}
    for u,o in enumerate(filter(lambda o:str(o)[0].isalpha() and str(o)[0].capitalize()==str(o)[0] and str(o) not in map(lambda a:a.strip(),'I, Piecewise, LambertW, Eq, symbols'.split(',')),dir())):
        print(f'@@@{u+1}.',o, type(o))
        # try:
        truth = False
        for tempt in range(budget:=5):
            try:
                truth = truth or tru.Verify(vars()[o]).verify() 
            except ValueError as ve:
                if (m:="math domain error") in str(ve):pass
                # elif(m:=)
                print("[ERROR]"+":"*99,m)
                # print(str(ve));1/0
        print("+"*8*8,*((truth,) if (b:=isinstance(truth,bool)) else (truth.items())),sep=('\n\t'if not b else ''))
        y[o] = truth
    print(*[yo for yo in y.items()],sep=('\n'))
    
def export_unfinished():
    # This might fail if run as a module without the __main__ block having run
    return locals().get('y', {})

