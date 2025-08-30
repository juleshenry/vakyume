from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np
from kwasak import kwasak_static

class SelectingPump:

    @kwasak_static
    def eqn_8_06(M: float = None, P_1: float = None, P_2: float = None, R: float = None, T: float = None, adiabatic_hp: float = None, k: float = None, w: float = None,**kwargs):
        return


    @staticmethod
    def eqn_8_06__M(P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        M = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*adiabatic_hp*(k - 1))
        result.append(M)
        return result

    @staticmethod
    def eqn_8_06__P_1(M: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_1 = P_2/(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_1)
        return result

    @staticmethod
    def eqn_8_06__P_2(M: float, P_1: float, R: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        P_2 = P_1*(1980000*M*adiabatic_hp/(R*T*w) - 1980000*M*adiabatic_hp/(R*T*k*w) + 1)**(k/(k - 1))
        result.append(P_2)
        return result

    @staticmethod
    def eqn_8_06__R(M: float, P_1: float, P_2: float, T: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        R = 1980000*M*adiabatic_hp*(k - 1)/(T*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(R)
        return result

    @staticmethod
    def eqn_8_06__T(M: float, P_1: float, P_2: float, R: float, adiabatic_hp: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        T = 1980000*M*adiabatic_hp*(k - 1)/(R*k*w*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(T)
        return result

    @staticmethod
    def eqn_8_06__adiabatic_hp(M: float, P_1: float, P_2: float, R: float, T: float, k: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        adiabatic_hp = R*T*k*w*((P_2/P_1)**((k - 1)/k) - 1)/(1980000*M*(k - 1))
        result.append(adiabatic_hp)
        return result

    @staticmethod
    def eqn_8_06__k(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, w: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_8_06__w(M: float, P_1: float, P_2: float, R: float, T: float, adiabatic_hp: float, k: float):
        # [.pyeqn] adiabatic_hp = (k / (k - 1) * (w * R * T) / (M * 550 * 3600) * ((P_2 / P_1) ** ((k - 1) / k) - 1))
        result = []
        w = 1980000*M*adiabatic_hp*(k - 1)/(R*T*k*((P_2/P_1)**((k - 1)/k) - 1))
        result.append(w)
        return result
