from math import log, sqrt, exp, pow, e
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak_static
from .eqn_7_14b__A_cap import eqn_7_14b__A
from .eqn_7_14b__Q_condensor_heat_duty_cap import eqn_7_14b__Q_condensor_heat_duty
from .eqn_7_14b__U_cap import eqn_7_14b__U
from .eqn_7_14b__del_T_1_cap import eqn_7_14b__del_T_1
from .eqn_7_14b__del_T_2_cap import eqn_7_14b__del_T_2

class Precondensors:
    eqn_7_14b__A = staticmethod(eqn_7_14b__A)
    eqn_7_14b__Q_condensor_heat_duty = staticmethod(eqn_7_14b__Q_condensor_heat_duty)
    eqn_7_14b__U = staticmethod(eqn_7_14b__U)
    eqn_7_14b__del_T_1 = staticmethod(eqn_7_14b__del_T_1)
    eqn_7_14b__del_T_2 = staticmethod(eqn_7_14b__del_T_2)

    @kwasak_static
    def eqn_7_14b(A=None, Q_condensor_heat_duty=None, U=None, del_T_1=None, del_T_2=None, **kwargs):
        return
