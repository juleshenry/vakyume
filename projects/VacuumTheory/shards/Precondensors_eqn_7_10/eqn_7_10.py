from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_10__L_c_P_cap import eqn_7_10__L_c_P
from .eqn_7_10__Q_condensor_heat_duty_cap import eqn_7_10__Q_condensor_heat_duty
from .eqn_7_10__del_T_cap import eqn_7_10__del_T

class Precondensors:
    eqn_7_10__L_c_P = eqn_7_10__L_c_P
    eqn_7_10__Q_condensor_heat_duty = eqn_7_10__Q_condensor_heat_duty
    eqn_7_10__del_T = eqn_7_10__del_T

    @kwasak_static
    def eqn_7_10(self, L_c_P=None, Q_condensor_heat_duty=None, del_T=None):
        return
