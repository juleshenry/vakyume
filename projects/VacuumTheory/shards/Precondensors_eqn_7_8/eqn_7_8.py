from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_8__L_c_cap import eqn_7_8__L_c
from .eqn_7_8__Q_condensor_heat_duty_cap import eqn_7_8__Q_condensor_heat_duty
from .eqn_7_8__c_p import eqn_7_8__c_p
from .eqn_7_8__del_T_cap import eqn_7_8__del_T

class Precondensors:
    eqn_7_8__L_c = staticmethod(eqn_7_8__L_c)
    eqn_7_8__Q_condensor_heat_duty = staticmethod(eqn_7_8__Q_condensor_heat_duty)
    eqn_7_8__c_p = staticmethod(eqn_7_8__c_p)
    eqn_7_8__del_T = staticmethod(eqn_7_8__del_T)

    @kwasak_static
    def eqn_7_8(L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, **kwargs):
        return
