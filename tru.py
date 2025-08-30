from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import inspect
from suck_consts import *
import random


def test_a():
    class VacuumTheory:

        @kwasak_static
        def eqn_1_3(
            m: float = None, T: float = None, k: float = None, v: float = None, **kwargs
        ):
            return

        @staticmethod
        def eqn_1_3__m(T: float, k: float, v: float):
            # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
            result = []
            m = 3.0 * T * k / v**2
            result.append(m)
            return result

        @staticmethod
        def eqn_1_3__T(k: float, m: float, v: float):
            # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
            result = []
            T = 0.333333333333333 * m * v**2 / k
            result.append(T)
            return result

        @staticmethod
        def eqn_1_3__k(T: float, m: float, v: float):
            # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
            result = []
            k = 0.333333333333333 * m * v**2 / T
            result.append(k)
            return result

        @staticmethod
        def eqn_1_3__v(T: float, k: float, m: float):
            # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
            result = []
            v = -1.73205080756888 * sqrt(T * k / m)
            result.append(v)
            v = 1.73205080756888 * sqrt(T * k / m)
            result.append(v)
            return result

        @kwasak_static
        def eqn_1_7(
            T: float = None,
            n: float = None,
            V: float = None,
            p: float = None,
            R: float = None,
            **kwargs,
        ):
            return

        @staticmethod
        def eqn_1_7__T(R: float, V: float, n: float, p: float):
            # [.pyeqn] p * V = n * R * T
            result = []
            T = V * p / (R * n)
            result.append(T)
            return result

        @staticmethod
        def eqn_1_7__n(R: float, T: float, V: float, p: float):
            # [.pyeqn] p * V = n * R * T
            result = []
            n = V * p / (R * T)
            result.append(n)
            return result

        @staticmethod
        def eqn_1_7__V(R: float, T: float, n: float, p: float):
            # [.pyeqn] p * V = n * R * T
            result = []
            V = R * T * n / p
            result.append(V)
            return result

        @staticmethod
        def eqn_1_7__p(R: float, T: float, V: float, n: float):
            # [.pyeqn] p * V = n * R * T
            result = []
            p = R * T * n / V
            result.append(p)
            return result

        @staticmethod
        def eqn_1_7__R(T: float, V: float, n: float, p: float):
            # [.pyeqn] p * V = n * R * T
            result = []
            R = V * p / (T * n)
            result.append(R)
            return result

        @kwasak_static
        def eqn_1_8(
            M: float = None,
            T: float = None,
            P: float = None,
            V: float = None,
            m: float = None,
            R: float = None,
            **kwargs,
        ):
            return

        @staticmethod
        def eqn_1_8__M(P: float, R: float, T: float, V: float, m: float):
            # [.pyeqn] P * V = m / M * R * T
            result = []
            M = R * T * m / (P * V)
            result.append(M)
            return result

        @staticmethod
        def eqn_1_8__T(M: float, P: float, R: float, V: float, m: float):
            # [.pyeqn] P * V = m / M * R * T
            result = []
            T = M * P * V / (R * m)
            result.append(T)
            return result

        @staticmethod
        def eqn_1_8__P(M: float, R: float, T: float, V: float, m: float):
            # [.pyeqn] P * V = m / M * R * T
            result = []
            P = R * T * m / (M * V)
            result.append(P)
            return result

        @staticmethod
        def eqn_1_8__V(M: float, P: float, R: float, T: float, m: float):
            # [.pyeqn] P * V = m / M * R * T
            result = []
            V = R * T * m / (M * P)
            result.append(V)
            return result

        @staticmethod
        def eqn_1_8__m(M: float, P: float, R: float, T: float, V: float):
            # [.pyeqn] P * V = m / M * R * T
            result = []
            m = M * P * V / (R * T)
            result.append(m)
            return result

        @staticmethod
        def eqn_1_8__R(M: float, P: float, T: float, V: float, m: float):
            # [.pyeqn] P * V = m / M * R * T
            result = []
            R = M * P * V / (T * m)
            result.append(R)
            return result

        @kwasak_static
        def eqn_1_9(
            T: float = None,
            M: float = None,
            rho: float = None,
            P: float = None,
            R: float = None,
            **kwargs,
        ):
            return

        @staticmethod
        def eqn_1_9__T(M: float, P: float, R: float, rho: float):
            # [.pyeqn] rho = P * M / (R * T)
            result = []
            T = M * P / (R * rho)
            result.append(T)
            return result

        @staticmethod
        def eqn_1_9__M(P: float, R: float, T: float, rho: float):
            # [.pyeqn] rho = P * M / (R * T)
            result = []
            M = R * T * rho / P
            result.append(M)
            return result

        @staticmethod
        def eqn_1_9__rho(M: float, P: float, R: float, T: float):
            # [.pyeqn] rho = P * M / (R * T)
            result = []
            rho = M * P / (R * T)
            result.append(rho)
            return result

        @staticmethod
        def eqn_1_9__P(M: float, R: float, T: float, rho: float):
            # [.pyeqn] rho = P * M / (R * T)
            result = []
            P = R * T * rho / M
            result.append(P)
            return result

        @staticmethod
        def eqn_1_9__R(M: float, P: float, T: float, rho: float):
            # [.pyeqn] rho = P * M / (R * T)
            result = []
            R = M * P / (T * rho)
            result.append(R)
            return result

        @kwasak_static
        def eqn_1_10(
            T_1: float = None,
            P_2: float = None,
            V_1: float = None,
            T_2: float = None,
            P_1: float = None,
            V_2: float = None,
            **kwargs,
        ):
            return

        @staticmethod
        def eqn_1_10__T_1(P_1: float, P_2: float, T_2: float, V_1: float, V_2: float):
            # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
            result = []
            T_1 = P_1 * T_2 * V_1 / (P_2 * V_2)
            result.append(T_1)
            return result

        @staticmethod
        def eqn_1_10__P_2(P_1: float, T_1: float, T_2: float, V_1: float, V_2: float):
            # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
            result = []
            P_2 = P_1 * T_2 * V_1 / (T_1 * V_2)
            result.append(P_2)
            return result

        @staticmethod
        def eqn_1_10__V_1(P_1: float, P_2: float, T_1: float, T_2: float, V_2: float):
            # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
            result = []
            V_1 = P_2 * T_1 * V_2 / (P_1 * T_2)
            result.append(V_1)
            return result

        @staticmethod
        def eqn_1_10__T_2(P_1: float, P_2: float, T_1: float, V_1: float, V_2: float):
            # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
            result = []
            T_2 = P_2 * T_1 * V_2 / (P_1 * V_1)
            result.append(T_2)
            return result

        @staticmethod
        def eqn_1_10__P_1(P_2: float, T_1: float, T_2: float, V_1: float, V_2: float):
            # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
            result = []
            P_1 = P_2 * T_1 * V_2 / (T_2 * V_1)
            result.append(P_1)
            return result

        @staticmethod
        def eqn_1_10__V_2(P_1: float, P_2: float, T_1: float, T_2: float, V_1: float):
            # [.pyeqn] P_1 * V_1 / T_1 = P_2 * V_2 / T_2
            result = []
            V_2 = P_1 * T_2 * V_1 / (P_2 * T_1)
            result.append(V_2)
            return result

        @kwasak_static
        def eqn_1_11(
            T: float = None,
            M: float = None,
            W: float = None,
            P: float = None,
            q: float = None,
            **kwargs,
        ):
            return

        @staticmethod
        def eqn_1_11__T(M: float, P: float, W: float, q: float):
            # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
            result = []
            T = 738 * M * P * q / (6821 * W)
            result.append(T)
            return result

        @staticmethod
        def eqn_1_11__M(P: float, T: float, W: float, q: float):
            # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
            result = []
            M = 6821 * T * W / (738 * P * q)
            result.append(M)
            return result

        @staticmethod
        def eqn_1_11__W(M: float, P: float, T: float, q: float):
            # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
            result = []
            W = 738 * M * P * q / (6821 * T)
            result.append(W)
            return result

        @staticmethod
        def eqn_1_11__P(M: float, T: float, W: float, q: float):
            # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
            result = []
            P = 6821 * T * W / (738 * M * q)
            result.append(P)
            return result

        @staticmethod
        def eqn_1_11__q(M: float, P: float, T: float, W: float):
            # [.pyeqn] q = W * (359 / M) * (760 / P) * (T / 492) * (1/60)
            result = []
            q = 6821 * T * W / (738 * M * P)
            result.append(q)
            return result

        @kwasak_static
        def eqn_1_12(
            sum_partial_pressures: float = None, Total_P: float = None, **kwargs
        ):
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
        def eqn_1_13a(y_a: float = None, n: float = None, n_a: float = None, **kwargs):
            return

        @staticmethod
        def eqn_1_13a__y_a(n: float, n_a: float):
            # [.pyeqn] y_a = n_a / n
            result = []
            y_a = n_a / n
            result.append(y_a)
            return result

        @staticmethod
        def eqn_1_13a__n(n_a: float, y_a: float):
            # [.pyeqn] y_a = n_a / n
            result = []
            n = n_a / y_a
            result.append(n)
            return result

        @staticmethod
        def eqn_1_13a__n_a(n: float, y_a: float):
            # [.pyeqn] y_a = n_a / n
            result = []
            n_a = n * y_a
            result.append(n_a)
            return result

        @kwasak_static
        def eqn_1_13b(p_a: float = None, y_a: float = None, P: float = None, **kwargs):
            return

        @staticmethod
        def eqn_1_13b__p_a(P: float, y_a: float):
            # [.pyeqn] y_a = p_a / P
            result = []
            p_a = P * y_a
            result.append(p_a)
            return result

        @staticmethod
        def eqn_1_13b__y_a(P: float, p_a: float):
            # [.pyeqn] y_a = p_a / P
            result = []
            y_a = p_a / P
            result.append(y_a)
            return result

        @staticmethod
        def eqn_1_13b__P(p_a: float, y_a: float):
            # [.pyeqn] y_a = p_a / P
            result = []
            P = p_a / y_a
            result.append(P)
            return result

    assert Verify(VacuumTheory).verify()
    print("exitoso")


def test_b():
    class C:
        @kwasak_static
        def eqn_10_10(
            rho: float = None,
            mu: float = None,
            bhp: float = None,
            bhp_0: float = None,
            **kwargs,
        ):
            return

        @staticmethod
        def eqn_10_10__rho(bhp: float, bhp_0: float, mu: float):
            # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
            # [Sympy Failover]
            pass  # Ollama offline

        @staticmethod
        def eqn_10_10__mu(bhp: float, bhp_0: float, rho: float):
            # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
            result = []
            mu = (
                -204374584201.104
                * I
                * (bhp / (bhp_0 * rho**0.84) - 0.5 / rho**0.84) ** (25 / 4)
            )
            result.append(mu)
            mu = (
                204374584201.104
                * I
                * (bhp / (bhp_0 * rho**0.84) - 0.5 / rho**0.84) ** (25 / 4)
            )
            result.append(mu)
            mu = -204374584201.104 * (
                bhp / (bhp_0 * rho**0.84) - 0.5 / rho**0.84
            ) ** (25 / 4)
            result.append(mu)
            mu = 204374584201.104 * (
                bhp / (bhp_0 * rho**0.84) - 0.5 / rho**0.84
            ) ** (25 / 4)
            result.append(mu)
            return result

        @staticmethod
        def eqn_10_10__bhp(bhp_0: float, mu: float, rho: float):
            # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
            result = []
            bhp = 0.0005 * bhp_0 * (31.0 * mu ** (4 / 25) * rho ** (21 / 25) + 1000.0)
            result.append(bhp)
            return result

        @staticmethod
        def eqn_10_10__bhp_0(bhp: float, mu: float, rho: float):
            # [.pyeqn] bhp = bhp_0 * (0.5 + 0.0155 * rho ** 0.84 * mu ** 0.16)
            result = []
            bhp_0 = 2000.0 * bhp / (31.0 * mu**0.16 * rho**0.84 + 1000.0)
            result.append(bhp_0)
            return result

    assert Verify(C).verify()
    print("exitosoB")


# print = lambda *a:a
class Verify:
    # iterate all methods and fill with dummy values.

    def __init__(self, lib_class):
        self.lib_class = lib_class

        # Get base equation names and their variants
        self.base_equations = self._get_base_equations(lib_class)
        self.equation_variants = self._get_equation_variants(lib_class)

        # Build dictionary of equations and their parameters
        self.equation_params = {
            base_eq: self.construir(base_eq) for base_eq in self.base_equations
        }

        # Initialize empty dict for dummy arguments
        self.proposed_dummy_args = {}

    def _get_base_equations(self, cls):
        """Get list of base equation names (those starting with 'eqn' without '__')"""
        return [
            name for name in dir(cls) if name.startswith("eqn") and "__" not in name
        ]

    def _get_equation_variants(self, cls):
        """Get list of equation variant names (those starting with 'eqn' with '__')"""
        return [name for name in dir(cls) if name.startswith("eqn") and "__" in name]

    def construir(self, base_equation_name):
        # Find the first two equation variants that start with the base equation name
        candidates = [
            variant
            for variant in self.equation_variants
            if variant.startswith(base_equation_name)
        ][:2]

        # Get the signatures of these candidate methods
        signatures = [
            inspect.signature(getattr(self.lib_class, candidate))
            for candidate in candidates
        ]

        # Function to extract parameter names from a signature string
        def extract_parameters(signature_str):
            # Clean the string: remove parentheses and commas with spaces
            cleaned = (
                str(signature_str).replace("(", "").replace(")", "").replace(", ", "")
            )
            # Split on TYPE (assuming TYPE is a defined separator, e.g., "float = None")
            parts = cleaned.split(TYPE)
            # Filter out empty parts
            return [part for part in parts if part]

        # Collect all unique parameter names from both signatures
        all_parameters = set()
        for sig in signatures:
            all_parameters.update(extract_parameters(sig))

        # Return sorted list of unique parameters
        return sorted(list(all_parameters))

    @staticmethod
    def make_rand():
        return round(random.random() * 4, 5)

    def verify(self):
        failing_dict = {}
        for equation_name, all_params in self.equation_params.items():
            self.proposed_dummy_args = {p: self.make_rand() for p in all_params}
            print("PROPOSEDUMBYARGZ" + str(self.proposed_dummy_args))
            for param_index, target_param in enumerate(all_params):
                if not self.todo_suave(
                    all_params,
                    target_param,
                    param_index,
                    self.proposed_dummy_args,
                    equation_name,
                ):
                    failing_dict[target_param] = f"{equation_name} remains elusive"
            self.proposed_dummy_args = {}
        return True if not failing_dict else failing_dict

    def non_none_soln(self, equation_name, equation_method, input_values):
        # Call the equation method with input values
        try:
            return getattr(self.lib_class(), equation_name)(**input_values)
        except OllamaOffline as oo:
            print(f"Ollama is offline for {equation_method}")
            return False

    def todo_suave(
        self, all_params, target_param, param_index, param_values, equation_name
    ) -> bool:
        # Get list of parameters except the target one
        other_params = [p for p in all_params if p != target_param]

        # Build dict of values for other parameters
        input_values = {param: param_values[param] for param in other_params}

        # Log the equation being called
        equation_method = f"{equation_name}_{target_param}"
        print(f"calling {equation_method} with {input_values}", end="")

        # Call the equation method with input values
        result = self.non_none_soln(equation_name, equation_method, input_values)
        if not param_index:  # assumes first is correct. not always true! TODO
            if not result:
                return False
            self.proposed_dummy_args[target_param] = (
                self.make_rand()
                if not result
                else result  ####TODO:!!!!!!!!!!!!!!!!!!!!!!!!!!!
            )  # should fix if first eqn is None ?
            print("\nAssuming below will be the golden n-tuple...")
            print(self.proposed_dummy_args)

        if isinstance(result, list):
            param_values[target_param] = result
        # Log the result
        print("-" * 9, ">>>", target_param, "=", result)

        def are_similar(a, b) -> bool:
            if (
                isinstance(c := a - b, float)
                or isinstance(c, int)
                or isinstance(c, complex)
            ):
                return abs(c) < 1e-10
            else:
                evaluated = result.subs({target_param: param_values[target_param]})
                if evaluated.is_real:
                    return evaluated < 1e-10
                else:
                    print(evaluated)
                    # Complex Yardstick
                    real, imag = map(lambda o: float(abs(o)), evaluated.as_real_imag())
                    print("RIIII" * 2, real, imag)
                    return real < 1e-10 and imag < 1e-10

        # ultimate copout: if indeed rez is None, it may be the unsolvable cases by IA
        if not result:
            return False
        for r in result:
            pseudo_gold = param_values[target_param]
            if isinstance(pseudo_gold, list):
                for pg in pseudo_gold:
                    if are_similar(r, pg):
                        param_values[target_param] = pg
                        return True
            elif are_similar(r, pseudo_gold):
                return True
        return False
        # return any(ev(r, param_values[target_param]) for r in result) if result else result


if __name__ == "__main__":
    test_b()
"""

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
    return y

"""
