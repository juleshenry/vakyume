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
    def eqn_8_08(P_1: float = None, P_2: float = None, adiabatic_power_watts: float = None, f: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_08__P_1(P_2: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        OllamaOffline('Ollama is offline')

    @staticmethod
    def eqn_8_08__P_2(P_1: float, adiabatic_power_watts: float, f: float):
        # [.pyeqn] adiabatic_power_watts = f / 12 * ((P_2 / P_1) ** 0.286 - 1)
        # [Sympy Failover]
        OllamaOffline('Ollama is offline')

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