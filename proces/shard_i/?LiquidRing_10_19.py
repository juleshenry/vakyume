from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve
from scipy.optimize import newton
from kwasak import kwasak_static
import pandas as pd
import numpy as np
from kwasak import kwasak_static

class LiquidRing:

    @kwasak_static
    def eqn_10_19(P: float = None, S_Th: float = None, S_p: float = None, T_e: float = None, T_i: float = None, p_c: float = None, p_s: float = None,**kwargs):
        return


    @staticmethod
    def eqn_10_19__P(S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        # [Sympy Failover]
        pass # Ollama offline

    @staticmethod
    def eqn_10_19__S_Th(P: float, S_p: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_Th = S_p/((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_Th)
        return result

    @staticmethod
    def eqn_10_19__S_p(P: float, S_Th: float, T_e: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        S_p = S_Th*((P*T_i + 460.0*P - T_i*p_s - 460.0*p_s)/(P*T_e + 460.0*P - T_e*p_c - 460.0*p_c))**(3/5)
        result.append(S_p)
        return result

    @staticmethod
    def eqn_10_19__T_e(P: float, S_Th: float, S_p: float, T_i: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_e = (P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P - T_i*p_s + 460.0*p_c*(S_p/S_Th)**(5/3) - 460.0*p_s)/((S_p/S_Th)**(5/3)*(P - p_c))
        result.append(T_e)
        T_e = 2.05280095711867*(P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        T_e = 2.05280095711867*(P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P - T_i*p_s + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*p_s)/((P - p_c)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(T_e)
        return result

    @staticmethod
    def eqn_10_19__T_i(P: float, S_Th: float, S_p: float, T_e: float, p_c: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        T_i = (P*T_e*(S_p/S_Th)**(5/3) + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P - T_e*p_c*(S_p/S_Th)**(5/3) - 460.0*p_c*(S_p/S_Th)**(5/3) + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        T_i = (0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P - 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*p_s)/(P - p_s)
        result.append(T_i)
        return result

    @staticmethod
    def eqn_10_19__p_c(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_s: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_c = (P*T_e*(S_p/S_Th)**(5/3) - P*T_i + 460.0*P*(S_p/S_Th)**(5/3) - 460.0*P + T_i*p_s + 460.0*p_s)/((S_p/S_Th)**(5/3)*(T_e + 460.0))
        result.append(p_c)
        p_c = 2.05280095711867*(0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        p_c = 2.05280095711867*(0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - P*T_i + 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 - 460.0*P + T_i*p_s + 460.0*p_s)/((T_e + 460.0)*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)
        result.append(p_c)
        return result

    @staticmethod
    def eqn_10_19__p_s(P: float, S_Th: float, S_p: float, T_e: float, T_i: float, p_c: float):
        # [.pyeqn] S_p = S_Th * ((P - p_s)*(460 + T_i)  / ( (P - p_c)*(460 + T_e) ))**0.6
        result = []
        p_s = (-P*T_e*(S_p/S_Th)**(5/3) + P*T_i - 460.0*P*(S_p/S_Th)**(5/3) + 460.0*P + T_e*p_c*(S_p/S_Th)**(5/3) + 460.0*p_c*(S_p/S_Th)**(5/3))/(T_i + 460.0)
        result.append(p_s)
        p_s = (-0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 - I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        p_s = (-0.487139289628747*P*T_e*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + P*T_i - 224.084073229223*P*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 460.0*P + 0.487139289628747*T_e*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5 + 224.084073229223*p_c*(-0.577350269189626*(S_p/S_Th)**0.333333333333333 + I*(S_p/S_Th)**0.333333333333333)**5)/(T_i + 460.0)
        result.append(p_s)
        return result
