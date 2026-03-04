from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_9__L_c_cap import eqn_7_9__L_c
from .eqn_7_9__Q_condensor_heat_duty_cap import eqn_7_9__Q_condensor_heat_duty
from .eqn_7_9__c_p import eqn_7_9__c_p
from .eqn_7_9__del_T_cap import eqn_7_9__del_T
from .eqn_7_9__rho import eqn_7_9__rho


class Precondensors:
    eqn_7_9__L_c = eqn_7_9__L_c
    eqn_7_9__Q_condensor_heat_duty = eqn_7_9__Q_condensor_heat_duty
    eqn_7_9__c_p = eqn_7_9__c_p
    eqn_7_9__del_T = eqn_7_9__del_T
    eqn_7_9__rho = eqn_7_9__rho

    @kwasak
    def eqn_7_9(
        self, L_c=None, Q_condensor_heat_duty=None, c_p=None, del_T=None, rho=None
    ):
        return
