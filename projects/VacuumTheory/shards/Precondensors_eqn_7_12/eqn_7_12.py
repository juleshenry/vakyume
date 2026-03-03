from cmath import log, sqrt, exp
from math import e, pi
from sympy import I, Piecewise, LambertW, Eq, symbols, solve, powsimp
from scipy.optimize import newton
import numpy as np
from vakyume.config import UnsolvedException

from vakyume.kwasak import kwasak
from .eqn_7_12__A_cap import eqn_7_12__A
from .eqn_7_12__Q_cap_condensor_heat_duty import eqn_7_12__Q_condensor_heat_duty
from .eqn_7_12__U_cap import eqn_7_12__U
from .eqn_7_12__del_T_cap import eqn_7_12__del_T


class Precondensors:
    eqn_7_12__A = eqn_7_12__A
    eqn_7_12__Q_condensor_heat_duty = eqn_7_12__Q_condensor_heat_duty
    eqn_7_12__U = eqn_7_12__U
    eqn_7_12__del_T = eqn_7_12__del_T

    @kwasak
    def eqn_7_12(self, A=None, Q_condensor_heat_duty=None, U=None, del_T=None):
        return
