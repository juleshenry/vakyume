from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_14b__A import eqn_7_14b__A
from .eqn_7_14b__Q_condensor_heat_duty import eqn_7_14b__Q_condensor_heat_duty
from .eqn_7_14b__U import eqn_7_14b__U
from .eqn_7_14b__del_T_1 import eqn_7_14b__del_T_1
from .eqn_7_14b__del_T_2 import eqn_7_14b__del_T_2


class Precondensors:
    eqn_7_14b__A = eqn_7_14b__A
    eqn_7_14b__Q_condensor_heat_duty = eqn_7_14b__Q_condensor_heat_duty
    eqn_7_14b__U = eqn_7_14b__U
    eqn_7_14b__del_T_1 = eqn_7_14b__del_T_1
    eqn_7_14b__del_T_2 = eqn_7_14b__del_T_2

    @kwasak
    def eqn_7_14b(
        self, A=None, Q_condensor_heat_duty=None, U=None, del_T_1=None, del_T_2=None
    ):
        return
