from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_8__L_c import eqn_7_8__L_c
from .eqn_7_8__Q_condensor_heat_duty import eqn_7_8__Q_condensor_heat_duty
from .eqn_7_8__c_p import eqn_7_8__c_p
from .eqn_7_8__del_T import eqn_7_8__del_T

class Precondensors:
    eqn_7_8__L_c = eqn_7_8__L_c
    eqn_7_8__Q_condensor_heat_duty = eqn_7_8__Q_condensor_heat_duty
    eqn_7_8__c_p = eqn_7_8__c_p
    eqn_7_8__del_T = eqn_7_8__del_T

    @kwasak
    def eqn_7_8(self, L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None):
        return
