from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class VacuumTheory:
    @kwasak_static
    def eqn_1_3(T=None, k=None, m=None, v=None, **kwargs):
        return

    @staticmethod
    def eqn_1_3__T(k: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        T = 0.333333333333333*m*v**2/k
        result.append(T)
        return result

    @staticmethod
    def eqn_1_3__k(T: float, m: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T => k = (m*v^2) / (3*T)
        result = []
        if T > 0 and m > 0 and sqrt(v).is_real:
            k = (m * pow(v, 2)) / (3.0 * T)
            result.append(k)
        else:
            # Handle cases where input values might lead to errors or non-physical results
            return None
        return [k] if len(result) == 1 else []


    @staticmethod
    def eqn_1_3__m(T: float, k: float, v: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T => m = (k*T) / ((v^2)/2)
        result = []
        if T > 0 and k > 0 and pow(v, 2):
            # Handle cases where input values might lead to errors or non-physical results
            return None

        m = (k * T) / ((pow(v, 2)) / 2.0) if len(result) == 1 else []
        result.append(m) if not np.isnan(m).any() and not np.isinf(m).any() else None
        return [m] if m is not None else []

    # Assuming eqn_10__k was meant to be corrected here, as it's referenced but undefined in the provided code snippet:

    @staticmethod
    def eqn_1_3__v(T: float, k: float, m: float):
        # [.pyeqn] .5 * m * v**2 = 1.5 * k * T
        result = []
        v = -1.73205080756888*sqrt(T*k/m)
        result.append(v)
        v = 1.73205080756888*sqrt(T*k/m)
        result.append(v)
        return result

