from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np
from kwasak import kwasak_static
from suck_consts import *
class SelectingPump:

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
        OllamaOffline('Ollama is offline')

    @staticmethod
    def eqn_8_01__installation_cost(NC: float, NS: float, SCON: float):
        # [.pyeqn] installation_cost = 16000 * (NS + 2 * NC) * (SCON / 1000) ** 0.35
        result = []
        installation_cost = 1426.00150101399*SCON**(7/20)*(2.0*NC + NS)
        result.append(installation_cost)
        return result
