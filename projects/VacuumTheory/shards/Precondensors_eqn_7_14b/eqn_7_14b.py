from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton, brentq
import numpy as np
from vakyume.config import UnsolvedException, safe_brentq

from vakyume.kwasak import kwasak
from .eqn_7_14b__A_cap import eqn_7_14b__A
from .eqn_7_14b__Q_cap_condensor_heat_duty import eqn_7_14b__Q_condensor_heat_duty
from .eqn_7_14b__U_cap import eqn_7_14b__U
from .eqn_7_14b__del_T_cap_1 import eqn_7_14b__del_T_1
from .eqn_7_14b__del_T_cap_2 import eqn_7_14b__del_T_2

class Precondensors:
    eqn_7_14b__A = eqn_7_14b__A
    eqn_7_14b__Q_condensor_heat_duty = eqn_7_14b__Q_condensor_heat_duty
    eqn_7_14b__U = eqn_7_14b__U
    eqn_7_14b__del_T_1 = eqn_7_14b__del_T_1
    eqn_7_14b__del_T_2 = eqn_7_14b__del_T_2

    @kwasak
    def eqn_7_14b(self, A=None, Q_condensor_heat_duty=None, U=None, del_T_1=None, del_T_2=None):
        return
