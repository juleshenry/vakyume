from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import numpy as np

class SelectingPump:
    @kwasak_static
    def eqn_8_1(NC=None, NS=None, SCON=None, installation_cost=None, **kwargs):
        return

    @staticmethod
    def eqn_8_1__NC(NS: float, SCON: float, installation_cost: float, **kwargs):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NC = -0.5*NS + 0.000350630766969363*installation_cost/SCON**(7/20)
        result.append(NC)
        return result

    @staticmethod
    def eqn_8_1__NS(NC: float, SCON: float, installation_cost: float, **kwargs):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        NS = -2.0*NC + 0.000701261533938727*installation_cost/SCON**(7/20)
        result.append(NS)
        return result

    @staticmethod
    def eqn_8_1__SCON(NC: float, NS: float, installation_cost: float, **kwargs):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        # Error during Sympy solve: Sympy solve failed
        def func(SCON):
            # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
            return eval("(16000 * (NS + 2 * NC) * (x / 1000) ** 0.35) - (installation_cost)".replace('x', str(SCON)))
        raise UnsolvedException("Pending LLM/Manual Repair")

    @staticmethod
    def eqn_8_1__installation_cost(NC: float, NS: float, SCON: float, **kwargs):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        installation_cost = 1426.00150101399*SCON**(7/20)*(2.0*NC + NS)
        result.append(installation_cost)
        return result

